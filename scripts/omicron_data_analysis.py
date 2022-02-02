import joey_utils as ut
import argparse
import pandas as pd
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from os.path import basename, join, normpath
import pyrosetta as pr

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-od", "--out_dir", type=str, required=True)
	parser.add_argument("-t", "--threshold", type=float, default=1.4)
	parser.add_argument("-p", "--plot", type=bool, default=True)
	args = parser.parse_args()
	return args


class model_pdbs():
	def __init__(self, od):
		self.name = basename(normpath(od))
		print(self.name)
		self.pdb, self.antibody = self.name.split('_')
		self.rpk_cons_wt = None
		self.rpk_cons_om = None
		self.rpk_free_wt = None
		self.rpk_free_om = None
		self.af2_cons_wt = None
		self.af2_cons_om = None
		self.af2_free_wt = None
		self.af2_free_om = None
		self.missing = []

		# Determine PDB names
		rpk_cons_wt = ut.output_file_name('{}_rpk_cons_wt.pdb'.format(self.name), od)
		rpk_cons_om = ut.output_file_name('{}_rpk_cons_om.pdb'.format(self.name), od)
		rpk_free_wt = ut.output_file_name('{}_rpk_free_wt.pdb'.format(self.name), od)
		rpk_free_om = ut.output_file_name('{}_rpk_free_om.pdb'.format(self.name), od)
		af2_cons_wt = ut.output_file_name('{}_af2_cons_wt.pdb'.format(self.name), od)
		af2_cons_om = ut.output_file_name('{}_af2_cons_om.pdb'.format(self.name), od)
		af2_free_wt = ut.output_file_name('{}_af2_free_wt.pdb'.format(self.name), od)
		af2_free_om = ut.output_file_name('{}_af2_free_om.pdb'.format(self.name), od)

		# Get PDB list
		all_pdbs = glob(join(od, '*.pdb'))

		# Check present PDBs
		if rpk_cons_wt in all_pdbs:
			self.rpk_cons_wt = rpk_cons_wt
		else:
			self.missing.append('rpk_cons_wt')
		if rpk_cons_om in all_pdbs:
			self.rpk_cons_om = rpk_cons_om
		else:
			self.missing.append('rpk_cons_om')
		if rpk_free_wt in all_pdbs:
			self.rpk_free_wt = rpk_free_wt
		else:
			self.missing.append('rpk_free_wt')
		if rpk_free_om in all_pdbs:
			self.rpk_free_om = rpk_free_om
		else:
			self.missing.append('rpk_free_om')
		if af2_cons_wt in all_pdbs:
			self.af2_cons_wt = af2_cons_wt
		else:
			self.missing.append('af2_cons_wt')
		if af2_cons_om in all_pdbs:
			self.af2_cons_om = af2_cons_om
		else:
			self.missing.append('af2_cons_om')
		if af2_free_wt in all_pdbs:
			self.af2_free_wt = af2_free_wt
		else:
			self.missing.append('af2_free_wt')
		if af2_free_om in all_pdbs:
			self.af2_free_om = af2_free_om
		else:
			self.missing.append('af2_free_om')


class pdb_struc_info():
	def __init__(self, model_pdbs_obj, attr, variant=True):
		self.pdb = model_pdbs_obj.pdb
		self.antibody = model_pdbs_obj.antibody
		method, constraint, variant = attr.split('_')
		self.method = method
		self.constraint = constraint
		if variant:
			self.variant = variant


def missing_models_table(pdbs, od):
	missing_models = []
	for attr, val in pdbs.__dict__.items():
		# Skip non-PDB attributes
		if '_' not in attr:
			continue

		# Check if attribute is none, add missing to list
		if val == None:
			missing_mod_info = pdb_struc_info(pdbs, attr).__dict__
			missing_models.append(missing_mod_info)

	# If any are missing, make a missing models csv
	if len(missing_models) > 0:
		df = pd.DataFrame(missing_models)
		df.to_csv(ut.output_file_name(
			'{}_missing_models'.format(pdbs.name), od, 'csv'), index=False)


def e_table_colunns():
	etab_cols = ['pdb', 'antibody', 'method', 'constrained', 'site', 
		'wt', 'om', 'mutated', 'wt_tot', 'om_tot', 'd_tot', 
		'wt_int', 'om_int', 'd_int']
	return etab_cols


def main_table_columns():
	cols = ['pdb', 'antibody', 'method', 'constrained', 
		'wt_total', 'om_total', 'd_total', 
		'wt_interface', 'om_interface', 'd_interface', 
		'tot_stab_mut_close', 'tot_stab_mut_distal', 
		'tot_stab_other_close', 'tot_stab_other_distal', 
		'tot_destab_mut_close', 'tot_destab_other_close', 
		'tot_destab_mut_distal', 'tot_destab_other_distal', 
		'int_destab_mut', 'int_destab_other', 
		'int_stab_mut', 'int_stab_other']
	return cols


def per_residue_energy_table(model_info, wt_pose, om_pose, sf, antigen, antibody):
	# Calculate per-res total energies
	wt_per_res_e = np.array(ut.calc_single_res_Es(wt_pose, sf, antigen))
	om_per_res_e = np.array(ut.calc_single_res_Es(om_pose, sf, antigen))
	delta_per_res_e = om_per_res_e - wt_per_res_e

	# Calculate per-res interface energies	
	wt_per_res_inter_e = np.array(
		ut.calc_single_res_intEs(wt_pose, sf, antigen, antibody))
	om_per_res_inter_e = np.array(
		ut.calc_single_res_intEs(om_pose, sf, antigen, antibody))
	delta_per_res_inter_e = om_per_res_inter_e - wt_per_res_inter_e

	# Get list of all antigen residues, correcting AF2 numbering
	pose_res_list = ut.selector_to_list(wt_pose, antigen)
	reference_res_list = [i[0] for i in 
		ut.selector_to_list(wt_pose, antigen, False)]
	if reference_res_list[0] == 1:
		reference_res_list = [i + 319 for i in reference_res_list]

	res_energy_dicts_list = []
	for res in range(len(pose_res_list)):
		res_e = {k:v for k,v in model_info.items()}
		res_e['site'] = reference_res_list[res]
		res_e['wt'] = ut.find_res_aa(wt_pose, pose_res_list[res])
		res_e['om'] = ut.find_res_aa(om_pose, pose_res_list[res])
		res_e['mutated'] = int(res_e['om'] != res_e['wt'])
		res_e['wt_tot'] = wt_per_res_e[res]
		res_e['om_tot'] = om_per_res_e[res]
		res_e['d_tot'] = delta_per_res_e[res]
		res_e['wt_int'] = wt_per_res_inter_e[res]
		res_e['om_int'] = om_per_res_inter_e[res]
		res_e['d_int'] = delta_per_res_inter_e[res]
		res_energy_dicts_list.append(res_e)

	return pd.DataFrame(res_energy_dicts_list)


def make_per_res_total_energy_lineplot(pdbs, od, method, constrained, 
	res_e_tab):
	# Get X tick labels, separating different chains by 1000
	x = []
	for res in res_e_tab['site']:
		if res not in x:
			x.append(res)
		else:
			x.append(res + 1000)
	
	# Create line plots
	fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)
	axes[0].plot(x, list(res_e_tab['wt_tot']), linewidth=0.3, color='blue')
	axes[1].plot(x, list(res_e_tab['om_tot']), linewidth=0.3, color='red')
	axes[2].plot(x, list(res_e_tab['d_tot']), linewidth=0.3, color='black')
	axes[0].set_ylabel('WT Total Energy', fontsize=8)
	axes[1].set_ylabel('Om Total Energy', fontsize=8)
	axes[2].set_ylabel('Energy Total Difference', fontsize=8)

	# Adjust figure size
	fig.set_size_inches(7.5, 5)
	fig.set_dpi(400)
	bottom, top, left, right = [0.1, 0.95, 0.1, 0.95]
	vcenter = bottom + (top - bottom)/2
	vlower = bottom + (vcenter - bottom)/3
	vupper = vcenter + 2*(top - vcenter)/3
	hcenter = left + (right - left)/2
	fig.tight_layout()
	fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)

	# Add labels
	fig.text(0.001, vcenter, 'Residue interaction energy', fontsize=10, 
		rotation='vertical', ha='left', va='center')
	fig.text(hcenter, 0.001, 'Residue Number', fontsize=10, ha='center', 
		va='bottom')
	fig.text(hcenter, 0.99, 'Per-Residue Interaction Energies', fontsize=11, 
		ha='center', va='top')
	
	# Save figure
	plt_name = ut.output_file_name('{}_{}_{}_total_res_energies'.format(
		pdbs.name, method, constrained), path=od, extension='png')
	plt.savefig(plt_name)
	plt.close()


def make_per_res_interface_energy_lineplot(pdbs, od, method, constrained, 
	res_e_tab):
	# Get X tick labels, separating different chains by 1000
	x = []
	for res in res_e_tab['site']:
		if res not in x:
			x.append(res)
		else:
			x.append(res + 1000)

	# Create line plots
	fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)
	axes[0].plot(x, list(res_e_tab['wt_int']), linewidth=0.3, color='blue')
	axes[1].plot(x, list(res_e_tab['om_int']), linewidth=0.3, color='red')
	axes[2].plot(x, list(res_e_tab['d_int']), linewidth=0.3, color='black')
	axes[0].set_ylabel('WT Inter-chain Energy', fontsize=8)
	axes[1].set_ylabel('Om Inter-chain Energy', fontsize=8)
	axes[2].set_ylabel('Inter-chain Energy Difference', fontsize=8)

	# Adjust figure size
	fig.set_size_inches(7.5, 5)
	fig.set_dpi(400)
	bottom, top, left, right = [0.1, 0.95, 0.1, 0.95]
	vcenter = bottom + (top - bottom)/2
	vlower = bottom + (vcenter - bottom)/3
	vupper = vcenter + 2*(top - vcenter)/3
	hcenter = left + (right - left)/2
	fig.tight_layout()
	fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)

	# Add labels
	fig.text(0.001, vcenter, 'Residue interaction energy', fontsize=10, 
		rotation='vertical', ha='left', va='center')
	fig.text(hcenter, 0.001, 'Residue Number', fontsize=10, ha='center', 
		va='bottom')
	fig.text(hcenter, 0.99, 'Per-Residue Interaction Energies', fontsize=11,
		ha='center', va='top')

	# Save figure
	plt_name = ut.output_file_name('{}_{}_{}_interface_res_energies'.format(
		pdbs.name, method, constrained), path=od, extension='png')
	plt.savefig(plt_name)
	plt.close()


def calculat_aggregate_energies(model_info, single_res_energy_table):
	# Calculate total energies
	model_info['wt_total'] = single_res_energy_table['wt_tot'].sum()
	model_info['om_total'] = single_res_energy_table['om_tot'].sum()
	model_info['d_total'] = model_info['om_total'] - model_info['wt_total']

	# Calculate total interface energies
	model_info['wt_interface'] = single_res_energy_table['wt_int'].sum()
	model_info['om_interface'] = single_res_energy_table['om_int'].sum()
	model_info['d_interface'] = model_info['om_interface'] - \
		model_info['wt_interface']

	return model_info


def get_antigen_interface_list(wt_pose, om_pose, antigen_interface):
	select_res_list_wt = [i[0] for i in ut.selector_to_list(wt_pose, antigen_interface, False)]
	select_res_list_om = [i[0] for i in ut.selector_to_list(om_pose, antigen_interface, False)]
	return list(set(select_res_list_wt).union(set(select_res_list_om)))


def list_mutated_sites():
	mut_sites = [67, 69, 70, 95, 142, 143, 144, 145, 211, 212, 214, 
		339, 371, 373, 375, 417, 440, 446, 477, 478, 484, 493, 496, 
		498, 501, 505, 547, 614, 655, 679, 681, 764, 796, 856, 954, 
		969, 981]
	return mut_sites


def identify_influential_residues(model_info, single_res_energy_table, 
	ant_int_res, threshold=1.4):
	# Shorthand for name
	se_df = single_res_energy_table

	# Identify mutation sites
	mut_sites = list_mutated_sites()

	# Get lists of imfluencial residues
	tot_stab = se_df[(se_df['d_tot'] < -1 * threshold)]
	tot_stab_mut_close =  tot_stab[(tot_stab['site'].isin(mut_sites)) & 
		(tot_stab['site'].isin(ant_int_res))]['site'].unique()
	tot_stab_mut_distal = tot_stab[(tot_stab['site'].isin(mut_sites)) & 
		~(tot_stab['site'].isin(ant_int_res))]['site'].unique()
	tot_stab_other_close =  tot_stab[~(tot_stab['site'].isin(mut_sites)) & 
		(tot_stab['site'].isin(ant_int_res))]['site'].unique()
	tot_stab_other_distal = tot_stab[~(tot_stab['site'].isin(mut_sites)) & 
		~(tot_stab['site'].isin(ant_int_res))]['site'].unique()

	tot_destab = se_df[(se_df['d_tot'] > threshold)]
	tot_destab_mut_close =  tot_destab[(tot_destab['site'].isin(mut_sites)) & 
		(tot_destab['site'].isin(ant_int_res))]['site'].unique()
	tot_destab_mut_distal = tot_destab[(tot_destab['site'].isin(mut_sites)) & 
		~(tot_destab['site'].isin(ant_int_res))]['site'].unique()
	tot_destab_other_close =  tot_destab[~(tot_destab['site'].isin(mut_sites)) & 
		(tot_destab['site'].isin(ant_int_res))]['site'].unique()
	tot_destab_other_distal = tot_destab[~(tot_destab['site'].isin(mut_sites)) & 
		~(tot_destab['site'].isin(ant_int_res))]['site'].unique()

	int_stab = se_df[(se_df['d_int'] < -1 * threshold)]
	int_stab_mut =  int_stab[
		(int_stab['site'].isin(mut_sites))]['site'].unique()
	int_stab_other =  int_stab[
		~(int_stab['site'].isin(mut_sites))]['site'].unique()

	int_destab = se_df[(se_df['d_int'] > threshold)]
	int_destab_mut =  int_destab[
		(int_destab['site'].isin(mut_sites))]['site'].unique()
	int_destab_other =  int_destab[
		~(int_destab['site'].isin(mut_sites))]['site'].unique()

	# Add the influential residues to the energy comparisons dict
	model_info['tot_stab_mut_distal'] = ut.join_list(tot_stab_mut_distal, ';') 
	model_info['tot_stab_other_close'] = ut.join_list(tot_stab_other_close, ';') 
	model_info['tot_stab_other_distal'] = \
		ut.join_list(tot_stab_other_distal, ';') 
	model_info['tot_destab_mut_close'] = ut.join_list(tot_destab_mut_close, ';') 
	model_info['tot_destab_other_close'] = \
		ut.join_list(tot_destab_other_close, ';') 
	model_info['tot_destab_mut_distal'] = \
		ut.join_list(tot_destab_mut_distal, ';') 
	model_info['tot_destab_other_distal'] = \
		ut.join_list(tot_destab_other_distal, ';') 
	model_info['int_destab_mut'] = ut.join_list(int_destab_mut, ';') 
	model_info['int_destab_other'] = ut.join_list(int_destab_other, ';') 
	model_info['int_stab_mut'] = ut.join_list(int_stab_mut, ';') 
	model_info['int_stab_other'] = ut.join_list(int_stab_other, ';')  

	return model_info


def significant_interface_site_barplots(pdbs, od, site_groups):
	for grp, ms in site_groups.items():
		if len(ms) == 0:
			continue
		x = np.arange(len(ms))  # the label locations
		width = 0.15  # the width of the bars

		fig, ax = plt.subplots()
		rects1 = ax.bar(x - 3*width/2, [i[0] for i in ms.values()], width, label='RRMC')
		rects2 = ax.bar(x - width/2, [i[1] for i in ms.values()], width, label='RRMF')
		rects4 = ax.bar(x + width/2, [i[2] for i in ms.values()], width, label='AFRC')
		rects3 = ax.bar(x + 3*width/2, [i[3] for i in ms.values()], width, label='AFRF')
		
		rect = patches.Rectangle((min(x) - 2 * width, -1.4), 
			max(x) - min(x) + 4 * width, 2.8, facecolor=(.25, .25, .25, 0.75))
		ax.add_patch(rect)

		ax.set_ylabel('Interface Energy Change', fontsize=12)
		ax.set_xticks(x, ms.keys())
		ax.tick_params(axis='both', which='major', labelsize=12)
		ax.legend(fontsize=12)     
		plt.xticks(rotation=45, ha='center')

		fig.set_size_inches(7.5, 5)
		fig.set_dpi(350)
		bottom, top, left, right = [0.1, 0.95, 0.1, 0.95]
		vcenter = bottom + (top - bottom)/2
		vlower = bottom + (vcenter - bottom)/3
		vupper = vcenter + 2*(top - vcenter)/3
		hcenter = left + (right - left)/2
		fig.tight_layout()
		fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)

		plt_name = ut.output_file_name('{}_{}_interface_res_energies'.format(
			pdbs.name, grp), path=od, extension='png')
		plt.savefig(plt_name)
		plt.close()


def consensus_site_evaluation(pdbs, od, res_e_tab, plot=True, threshold=1.4):
	# Identify influential interfce sites
	major_inter_sites = res_e_tab[(res_e_tab['d_int'] > threshold) | 
		(res_e_tab['d_int'] < -1 * threshold)]['site'].unique()
	
	# Collect 4-member list of influential interface energies 
	# (for each method/constraint)
	method_singles = {}
	for dis in sorted(major_inter_sites):
		single_es = []
		rdf = res_e_tab[res_e_tab['site'] == dis]
		rc = rdf[(rdf['method'] == 'rpk') & (rdf['constrained'] == 'cons')]
		if len(rc) == 1:
			single_es.append(float(rc['d_int']))
		elif len(rc) == 0:
			single_es.append(0)
		else:
			single_es.append(float(rc['d_int'].sum()))
		rf = rdf[(rdf['method'] == 'rpk') & (rdf['constrained'] == 'free')]
		if len(rf) == 1:
			single_es.append(float(rf['d_int']))
		elif len(rf) == 0:
			single_es.append(0)
		else:
			single_es.append(float(rf['d_int'].sum()))
		ac = rdf[(rdf['method'] == 'af2') & (rdf['constrained'] == 'cons')]
		if len(ac) == 1:
			single_es.append(float(ac['d_int']))
		elif len(ac) == 0:
			single_es.append(0)
		else:
			single_es.append(float(ac['d_int'].sum()))
		af = rdf[(rdf['method'] == 'af2') & (rdf['constrained'] == 'free')]
		if len(af) == 1:
			single_es.append(float(af['d_int']))
		elif len(af) == 0:
			single_es.append(0)
		else:
			single_es.append(float(af['d_int'].sum()))
							 
		method_singles[dis] = single_es
	
	# Sort sites as mutant or not	
	mut_sites = list_mutated_sites()			 
	site_groups = {}
	site_groups['mutated'] = {k: v for k,v in method_singles.items() if k in mut_sites}
	site_groups['non-mutated'] = {k: v for k,v in method_singles.items() if k not in mut_sites}

	# Plotting significant sites as bars
	if plot:
		significant_interface_site_barplots(pdbs, od, site_groups)
	
	# Create consensus reports 
	res_reports = {}
	for res in mut_sites:
		if res not in site_groups['mutated'].keys():
			res_reports[res] = ''
		else:
			es = site_groups['mutated'][res]
			num_plus = sum([1 for i in es if i >= threshold])
			num_minus = sum([1 for i in es if i <= -1 * threshold])
			if num_plus > num_minus:
				ext = '+' * (num_plus - num_minus)
			elif num_plus < num_minus:
				ext = '-' * (-num_plus + num_minus)
			else:
				ext = '*'
			res_reports[res] = ext

	# Create table
	consensus_dict = {'pdb': pdbs.pdb, 'antibody': pdbs.antibody}
	consensus_dict.update(res_reports)
	consensus_df = pd.DataFrame(consensus_dict, index=[0])
	consensus_df.to_csv(ut.output_file_name(
		'{}_consensus_table'.format(pdbs.name), 
		path=od, extension='csv'), index=False)


def make_plots_and_tables(pdbs, od, threshold=1.4, plot=True):
	# Make score function and selections
	sf = ut.get_sf()
	antigen_selection = ut.chain_selector('A,B')
	antibody_selection = ut.not_selector(antigen_selection)
	antigen_interface = ut.one_side_close_and_interface_selector(
		antigen_selection, antibody_selection)

	# Initialize tables
	single_res_es_table = pd.DataFrame(columns=e_table_colunns())
	overall_summary_table = pd.DataFrame(columns=main_table_columns())

	methods = ['rpk', 'af2']
	constraints = ['cons', 'free']
	for method in methods:
		for constrained in constraints:
			# Get PDB names
			wt_base = '{}_{}_{}_wt'.format(pdbs.name, method, constrained)
			wt_pdb = ut.output_file_name(wt_base, path=od, extension='pdb')
			om_base = '{}_{}_{}_om'.format(pdbs.name, method, constrained)
			om_pdb = ut.output_file_name(om_base, path=od, extension='pdb')

			# Can't do comparison if a model is missing 
			if wt_base[7:] in pdbs.missing or om_base[7:] in pdbs.missing:
				continue

			# Get poses
			wt_pose = ut.load_pose(wt_pdb)
			om_pose = ut.load_pose(om_pdb)

			# Initialize energy collection
			energy_comparisons = {'pdb': pdbs.pdb, 'antibody': pdbs.antibody, 
				'method': method, 'constrained': constrained}

			# Get per-residue table, export (overwritten at each as checkpoint)
			se_df = per_residue_energy_table(energy_comparisons, wt_pose, 
				om_pose, sf, antigen_selection, antibody_selection)
			single_res_es_table = single_res_es_table.append(se_df, 
				ignore_index=True)
			single_res_es_table.to_csv(ut.output_file_name(
				'{}_single_res_energies'.format(pdbs.name), 
				path=od, extension='csv'), index=False)

			# Make per-residue lineplots
			make_per_res_total_energy_lineplot(pdbs, od, method, 
				constrained, se_df)
			make_per_res_interface_energy_lineplot(pdbs, od, method, 
				constrained, se_df)

			# Calculate aggregate energies
			energy_comparisons = calculat_aggregate_energies(
				energy_comparisons, se_df)

			# Identify influential residues
			antigen_interface_residue_list = \
				get_antigen_interface_list(wt_pose, om_pose, antigen_interface)
			energy_comparisons = identify_influential_residues(
				energy_comparisons, se_df, antigen_interface_residue_list, 
				threshold=threshold)

			# Add collected info to main table
			overall_summary_table = \
				overall_summary_table.append(energy_comparisons, ignore_index=True)

			# Export summary table (overwritten at each step as a checkpoint)
			overall_summary_table.to_csv(ut.output_file_name(
				'{}_energy_and_changed_residue_identification'.format(pdbs.name),  
				path=od, extension='csv'), index=False)

	# Perform consensus evaluation if all models are present
	if len(pdbs.missing) == 0:
		consensus_dict = consensus_site_evaluation(pdbs, od, 
			single_res_es_table, plot=plot, threshold=threshold)


def main(args):
	# Get PDBs
	pdbs = model_pdbs(args.out_dir)

	# Output missing models
	missing_models_table(pdbs, args.out_dir)

	# Make plots and tables
	make_plots_and_tables(pdbs, args.out_dir, threshold=args.threshold, 
		plot=args.plot)
			

if __name__ == '__main__':
	# Take user inputs
	args = parse_args()

	# Initializing PyRosetta
	ros_opts = ut.pyrosetta_init(verbose=False, extra_opts=['ignore_zero_occupancy false'])

	# Execute main
	main(args)