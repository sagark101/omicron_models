import argparse
import joey_utils as ut
import pandas as pd
import pyrosetta as pr

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-wt", "--wt_pdb", required=True, type=str,
        default='start_proteases/HCV.pdb', help="Pick starting PDB")
    parser.add_argument("-om", "--om_pdb", required=True, type=str,
        default='start_proteases/HCV.pdb', help="Pick starting PDB")
    parser.add_argument("-od", "--out_dir", type=str, 
        help="Name an output directory for the consolidated csv and decoys")
    parser.add_argument("-name", "--name", required=True, type=str)
    parser.add_argument("-suf", "--name_suffix", type=str,
        help="Add a suffix name to the output csv.")
    parser.add_argument("-site", "--design_focus", required=True, type=int, 
        nargs='+')
    parser.add_argument("-wd", "--wide_design", action="store_true")
    parser.add_argument('-n', "--n_decoys", type=int, default=100, 
        help="How many decoys do you want? (Default: 100)")
    args = parser.parse_args()
    return args


def determine_target_site_pose_residues(pose, target_res_list):
    # Determine pose chains
    pose_chains = ut.get_pose_chain_list(pose)

    # Get pose-numbered sites
    pose_res = []
    for chain in ['A', 'B']:
        if chain in pose_chains:
            for residue in target_res_list:
                pose_res.append(ut.pdb2pose(pose, residue, chain))

    return ut.index_selector(ut.join_list(pose_res, ','))


def get_hbonding_residues(site_selection):
    # Make main chain selections
    binding_selection = ut.chain_selector('H,L,X')

    # Make HBond selector
    hbonds = ut.hbond_selector(site_selection, include_sc=True, 
        include_bb_bb=False)
    bind_bonds = ut.selector_intersection(hbonds, binding_selection)

    return bind_bonds


def select_design_shell(wt_pose, wide_design=False):
    # Determine target site residues
    target_res = determine_target_site_pose_residues(wt_pose, args.design_focus)
    print('Target residues:', ut.selector_to_list(wt_pose, target_res, False))

    # Design larger shell around the hbonding residues
    if wide_design:
        # Determine shell around the H-bonding residues
        wide_shell = ut.close_interface_with_selection_selector(target_res)

        # Limit the shell to the antibody loops by excluding sheets and 
        # non-antibody chains
        ab_selector = ut.chain_selector('H,L')
        non_sheet = ut.secstruct_selector('HL')
        return ut.selector_intersection(wide_shell, ab_selector, non_sheet)

    # Design just H-bonding residues
    else:
        # Determine H-bonding residues
        hbond_res = get_hbonding_residues(target_res)
        print('Hbond residues:', ut.selector_to_list(wt_pose, hbond_res, False))

        # Convert H-bond residues selection to index selector, since Omicron 
        # model may not have WT's bonds
        hbond_res_list = ut.selector_to_list(wt_pose, hbond_res)
        return ut.index_selector(ut.join_list(hbond_res_list, ','))


def determine_sequence_changes(ref_pose, design_pose, chain):
    # Determine chain index
    pose_chains = ut.get_pose_chain_list(ref_pose)
    chain_ind = pose_chains.index(chain) + 1

    # Determine residue range of chain
    first_res = ref_pose.chain_begin(chain_ind)
    last_res = ref_pose.chain_end(chain_ind) + 1

    # Check for changes
    changes = []
    for res in range(first_res, last_res):
        wt = ut.find_res_aa(ref_pose, res)
        des = ut.find_res_aa(design_pose, res)
        if des != wt:
            pdb_number = ut.pose2pdb(ref_pose, res, include_chain=False)
            changes.append('{}{}{}'.format(wt, str(pdb_number), des))

    return changes


def main(args):
    # Loading poses
    wt_pose = ut.load_pose(args.wt_pdb)
    om_pose = ut.load_pose(args.om_pdb, res_type_cst=1)

    # Select design shell
    design_shell = select_design_shell(wt_pose, args.wide_design)
    print('Design shell:', ut.selector_to_list(wt_pose, design_shell, False))

    # Select repack shell
    repack_shell = ut.close_interface_with_selection_selector(design_shell)
    print('Repack shell:', ut.selector_to_list(wt_pose, repack_shell, False))

    # Set up mover
    sf = ut.get_sf(constrain=1.0)
    tf = ut.make_task_factory(design_shell, repack_shell)
    mm = ut.make_move_map(True, True, True)
    fr =  ut.fast_relax_mover(sf, tf, mm)
    #print(tf.create_task_and_apply_taskoperations(om_pose))

    # Set output names
    decoy_name = ut.output_file_name(args.name, path=args.out_dir, 
        suffix=args.name_suffix)
    print('Decoys:', decoy_name)
    table_name = ut.output_file_name(decoy_name, suffix='designs', 
        extension='csv')
    print('Table:', table_name)

    # Initialize table
    cols = ['Decoy', 'Score', 'H substitutions', 'L substitutions']
    df = pd.DataFrame(columns=cols)
    df.to_csv(table_name, index=False)

    # Run design protocol
    jd = ut.pyjobdistributor(decoy_name, args.n_decoys, sf)
    current_decoy = 1
    while not jd.job_complete:
        print('Decoy {}/{}\n'.format(current_decoy, args.n_decoys))
        pp = pr.Pose(om_pose)
        fr.apply(pp)
        current_decoy += 1

        # Tabulate changes
        decoy_change_dict = {}
        decoy_change_dict['Decoy'] = jd.current_name[:-4]
        decoy_change_dict['Score'] = sf(pp)
        decoy_change_dict['H substitutions'] = ut.join_list(
            determine_sequence_changes(om_pose, pp, 'H'), ';')
        decoy_change_dict['L substitutions'] = ut.join_list(
            determine_sequence_changes(om_pose, pp, 'L'), ';')
        df = pd.DataFrame(decoy_change_dict, columns=cols, index=[0])
        df.to_csv(table_name, index=False, mode='a', header=False)

        jd.output_decoy(pp)


if __name__ == '__main__':
    # Take user inputs
    args = parse_args()

    # Initializing PyRosetta
    ut.pyrosetta_init(verbose=False, extra_opts=['ignore_zero_occupancy false'])

    # Execute main
    main(args)