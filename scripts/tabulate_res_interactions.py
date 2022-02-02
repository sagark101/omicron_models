import argparse
import joey_utils as ut
import pandas as pd
import pyrosetta as pr
from os.path import basename

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", required=True, type=str)
    parser.add_argument("-c", "--csv", required=True, type=str)
    parser.add_argument("-o", "--outpath", required=True, type=str)
    args = parser.parse_args()
    return args

pr.init('-mute all')

abpdb = pd.read_csv('SARS-CoV-2Spike-AntibodyComplexes_12-01-2021_2.csv')

sabdab = pd.read_csv('sabdab_20211201_processed.csv')
sabdab['PDB'] = sabdab['PDB'].str.upper()
pdb_list = list(abpdb['ID'].unique())
sabdab.drop(sabdab[~(sabdab['PDB'].isin(pdb_list))].index, inplace=True)

sbd_summary = pd.read_csv('sabdab_20211201_summary_file.csv')
sbd_summary['pdb'] = sbd_summary['pdb'].str.upper()

args = parse_args()

cols = ['ID', 'Title', 'Released', 'Class of Spike-Binding Protein', 
    'Protein 1', 'Protein 2', 'Protein 3', 'Protein 4', 
    'Domain Target on Spike', 'Notes', 'antibody_number', 'antigen_chain', 
    'heavy_chain', 'light_chain', 'spike_residues', 'in_ntd', 'in_rbd', 
    'in_ace2bind', 'multichain', 'missing_antigen_res']

pdb = basename(args.pdb).replace('.pdb', '').upper()
print(pdb)

df = abpdb[abpdb['ID'] == pdb.upper()]

pdb_chain_associations = sbd_summary[sbd_summary['pdb'] == pdb]
pdb_chain_associations.sort_values(by='antigen_chain', inplace=True)

for antibody_number in range(len(pdb_chain_associations)):
    df = abpdb[abpdb['ID'] == pdb]
    df.loc[df['ID'] == pdb, 'antibody_number'] = antibody_number
    
    antigen_chain = pdb_chain_associations['antigen_chain'].iloc[antibody_number]
    antigen_chain = antigen_chain.split(' | ')
    antigen_chain = ut.join_list(sorted(antigen_chain), ',')
    multichain = int(len(antigen_chain) > 1)
    heavy_chain = pdb_chain_associations['Hchain'].iloc[antibody_number]
    light_chain = pdb_chain_associations['Lchain'].iloc[antibody_number]
    df.loc[df['ID'] == pdb, 'antigen_chain'] = antigen_chain
    df.loc[df['ID'] == pdb, 'multichain'] = multichain
    df.loc[df['ID'] == pdb, 'heavy_chain'] = heavy_chain
    df.loc[df['ID'] == pdb, 'light_chain'] = light_chain

    cdr_table = sabdab[(sabdab['PDB'] == pdb) & (sabdab['Fab'].str.contains(heavy_chain))]

    cdr_seqs = {}
    for cdrl in ['H1','H2','H3','L1','L2','L3']:
        cdr_seqs[cdrl] = cdr_table[cdrl].to_string(index=False)

    pose = ut.load_pose(args.pdb)
    pose_tab = ut.tabulate_pose_residues(pose)

    cdr_pdb_res = {}
    cdr_pose_res = {}
    for cdrl, seq in cdr_seqs.items():
        if seq == 'NaN':
            continue
        which_chain = {'H': heavy_chain, 'L': light_chain}[cdrl[0]]
        try:
            loop_table = ut.find_sequence_in_pose_table(pose_tab, seq)
        except Exception as e:
            print(e)
            pose_chain = pose_tab[pose_tab['pdb_chain'] == which_chain]['pose_chain'].unique()[0]
            aln = ut.global_align(seq, pose.chain_sequence(2), True)[0]
            aligned_indices = [n for n, i in enumerate(aln[0]) if i != '-']
            closest_seq = aln[1][min(aligned_indices):max(aligned_indices) + 1]
            loop_table = ut.find_sequence_in_pose_table(pose_tab, closest_seq)   
        loop_table = loop_table[loop_table['pdb_chain'] == which_chain]
        cdr_pdb_res[cdrl] = ut.join_list(loop_table['pdb_number'], ',')
        cdr_pose_res[cdrl] = ut.join_list(loop_table['pose_number'], ',') 

    cdr_residues = ut.join_list(cdr_pose_res.values(), ',')
    cdr_selection = ut.index_selector(cdr_residues)

    if multichain:
        spike_selection = ut.selector_union(
            *[ut.chain_selector(i) for i in antigen_chain.split(',')])
    else:
        spike_selection = ut.chain_selector(antigen_chain)

    contacted_spike_selection = ut.one_side_interface_selector(spike_selection, cdr_selection)
    spike_residues = [i[0] for i in ut.selector_to_list(pose, contacted_spike_selection, False)]
    spike_residues = list(set(spike_residues))
    df['spike_residues'] = ut.join_list(sorted(spike_residues), ';')

    ntd_range = list(range(14,306))
    rbd_range = list(range(319,542))
    ace2bind_range = sorted(list(range(438,507)))
    in_ntd, in_rbd, in_ace2bind = 0, 0, 0
    for res in spike_residues:
        if res in ntd_range:
            in_ntd = 1
        if res in rbd_range:
            in_rbd = 1
        if res in ace2bind_range:
            in_ace2bind = 1
    df.loc[df['ID'] == pdb, 'in_ntd'] = in_ntd
    df.loc[df['ID'] == pdb, 'in_rbd'] = in_rbd
    df.loc[df['ID'] == pdb, 'in_ace2bind'] = in_ace2bind

    antigen_tab = pose_tab[(pose_tab['pdb_chain'].isin(list(antigen_chain)))]
    if in_ntd == 1 and in_rbd == 1:
        antigen_residues_to_keep = \
            antigen_tab[antigen_tab['pdb_number'].isin(range(14, 542))]
    elif in_ntd == 1:
        antigen_residues_to_keep = \
            antigen_tab[antigen_tab['pdb_number'].isin(ntd_range)]
    else:
        antigen_residues_to_keep = \
            antigen_tab[antigen_tab['pdb_number'].isin(rbd_range)]
    antigen_res = list(antigen_residues_to_keep['pose_number'])
    antigen_to_keep_selection = ut.index_selector(antigen_res)
    antigen_pdb_numbers = list(antigen_residues_to_keep['pdb_number'])
    missing_antigen_res = [i for i in range(min(antigen_pdb_numbers), 
        max(antigen_pdb_numbers) + 1) if i not in antigen_pdb_numbers]
    df.loc[df['ID'] == pdb, 'missing_antigen_res'] = ut.join_list(missing_antigen_res, ';')

    antibody_residues_to_keep = pose_tab[(pose_tab['pdb_number'] < 120) &
        (pose_tab['pdb_chain'].isin([heavy_chain, light_chain]))]
    antibody_res = list(antibody_residues_to_keep['pose_number'])
    antibody_to_keep_selection = ut.index_selector(antibody_res)

    keep_chain_selection = ut.selector_union(antigen_to_keep_selection, antibody_to_keep_selection)
    trimmed_pose = ut.delete_selection_mover(pose, ut.not_selector(keep_chain_selection))

    trimmed_name = ut.output_file_name(pdb, args.outpath, extension='pdb', 
        suffix=antibody_number + 1)
    trimmed_pose.dump_pdb(trimmed_name)

    trimmed_pose_pdb = ut.tabulate_pdb_atom_lines(trimmed_name)
    if 'A' not in list(antigen_chain) or heavy_chain != 'H' or light_chain != 'L':
        trimmed_pose_pdb.loc[trimmed_pose_pdb['chain_id'] == heavy_chain,'chain_id'] = 'H'
        trimmed_pose_pdb.loc[trimmed_pose_pdb['chain_id'] == light_chain,'chain_id'] = 'L'
        for n, chain in enumerate(antigen_chain.split(',')):
            trimmed_pose_pdb.loc[trimmed_pose_pdb['chain_id'] == chain,'chain_id'] = \
                {0:'A', 1:'B'}[n]
        ut.atom_table_to_pdb(trimmed_pose_pdb, trimmed_name)

    df = df[cols]
    df.to_csv(args.csv, mode='a', header=False, index=False)
    print(df.to_dict())

    contiguity = ut.check_pdb_atom_table_contiguity(trimmed_pose_pdb)
    print(contiguity)
    if any(list(contiguity['insertions'].values())):
        insert_name = ut.output_file_name(trimmed_name, args.outpath, extension='txt', 
        suffix='insertions')
        with open(insert_name, 'a') as a:
            a.write('{}_{}\n'.format(pdb, antibody_number + 1))
