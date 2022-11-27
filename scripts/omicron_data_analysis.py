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
	parser.add_argument("-ab", "--antibody_name", type=str, required=True)
	parser.add_argument("-rcw", "--rpk_cons_wt", type=str)
	parser.add_argument("-rco", "--rpk_cons_om", type=str)
	parser.add_argument("-rfw", "--rpk_free_wt", type=str)
	parser.add_argument("-rfo", "--rpk_free_om", type=str)
	parser.add_argument("-acw", "--af2_cons_wt", type=str)
	parser.add_argument("-aco", "--af2_cons_om", type=str)
	parser.add_argument("-afw", "--af2_free_wt", type=str)
	parser.add_argument("-afo", "--af2_free_om", type=str)
	parser.add_argument("-od", "--out_dir", type=str, required=True)
	parser.add_argument("-t", "--threshold", type=float, default=1.4)
	parser.add_argument("-p", "--plot", type=bool, default=True)
	args = parser.parse_args()
	return args


class model_pdbs():
	def __init__(self, ab, rcw, rco, rfw, rfo, acw, aco, afw, afo):
		self.name = ab
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

		# Check present PDBs
		if rcw:
			self.rpk_cons_wt = rcw
		else:
			self.missing.append('rpk_cons_wt')
		if rco:
			self.rpk_cons_om = rco
		else:
			self.missing.append('rpk_cons_om')
		if rfw:
			self.rpk_free_wt = rfw
		else:
			self.missing.append('rpk_free_wt')
		if rfo:
			self.rpk_free_om = rfo
		else:
			self.missing.append('rpk_free_om')
		if acw:
			self.af2_cons_wt = acw
		else:
			self.missing.append('af2_cons_wt')
		if aco:
			self.af2_cons_om = aco
		else:
			self.missing.append('af2_cons_om')
		if afw:
			self.af2_free_wt = afw
		else:
			self.missing.append('af2_free_wt')
		if afo:
			self.af2_free_om = afo
		else:
			self.missing.append('af2_free_om')


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


def create_antigen_aligned_table(wt_pose, om_pose, antigen):
	# Create sequence alignment table
	wt_seq = ut.get_pose_sequence(wt_pose, antigen)
	om_seq = ut.get_pose_sequence(om_pose, antigen)
	alignment = ut.global_align(wt_seq, om_seq)[0]
	align_table = ut.tabulate_sequence_alignment(alignment)

	# Change to PDB numbering
	align_table.loc[align_table['s1_position'].notna(), 's1_position'] = \
		[i[0] for i in ut.selector_to_list(wt_pose, antigen, False)]
	align_table.loc[align_table['s2_position'].notna(), 's2_position'] = \
		[i[0] for i in ut.selector_to_list(om_pose, antigen, False)]

	return align_table


def aligned_per_res_energy_table(wt_pose, om_pose, sf, antigen, antibody):
	# Get initial per-res energy arrays
	## Per-res total energies
	wt_per_res_e = np.array(ut.calc_single_res_Es(wt_pose, sf, antigen))
	om_per_res_e = np.array(ut.calc_single_res_Es(om_pose, sf, antigen))
	## Per-res interface energies
	wt_per_res_int_e = np.array(
		ut.calc_single_res_intEs(wt_pose, sf, antigen, antibody))
	om_per_res_int_e = np.array(
		ut.calc_single_res_intEs(om_pose, sf, antigen, antibody))

	# Create alignment table with scores
	align_table = create_antigen_aligned_table(wt_pose, om_pose, antigen)
	align_table['wt_tot'] = 0
	align_table['om_tot'] = 0
	align_table['wt_int'] = 0
	align_table['om_int'] = 0
	align_table.loc[align_table['s1_position'].notna(), 'wt_tot'] = \
		wt_per_res_e
	align_table.loc[align_table['s2_position'].notna(), 'om_tot'] = \
		om_per_res_e
	align_table.loc[align_table['s1_position'].notna(), 'wt_int'] = \
		wt_per_res_int_e
	align_table.loc[align_table['s2_position'].notna(), 'om_int'] = \
		om_per_res_int_e
	align_table['d_tot'] = align_table.apply(
		lambda r: r['om_tot'] - r['wt_tot'], axis='columns')
	align_table['d_int'] = align_table.apply(
		lambda r: r['om_int'] - r['wt_int'], axis='columns')

	return align_table


def per_residue_energy_table(
		model_info, wt_pose, om_pose, sf, antigen, antibody):
	# Create sequence-aligned table with total and interface energies
	align_table = aligned_per_res_energy_table(
		wt_pose, om_pose, sf, antigen, antibody)

	# Add general info columns
	for k, v in model_info.items():
		align_table[k] = v

	# Trim columns
	col_list = list(model_info) + ['s1_position', 's1_sequence', 
		's2_sequence', 'identity', 'wt_tot', 'om_tot', 'd_tot', 
		'wt_int', 'om_int', 'd_int']
	align_table = align_table[col_list]

	# Rename columns
	replacements = {'s1_position': 'site', 's1_sequence': 'wt',
		's2_sequence': 'om', 'identity': 'mutated'}
	align_table.columns = [replacements[i] 
		if i in replacements else i for i in col_list]

	# Make site column all integers
	align_table['site'] = [ut.str2int(i) for i in align_table['site']]

	# Invert identity column to mutated
	align_table['mutated'] = [ut.str2int(not i) for i in align_table['mutated']]

	return align_table


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
	select_res_list_wt = [i[0] for i in 
		ut.selector_to_list(wt_pose, antigen_interface, False)]
	select_res_list_om = [i[0] for i in 
		ut.selector_to_list(om_pose, antigen_interface, False)]
	return sorted(list(set(select_res_list_wt).union(
		set(select_res_list_om))))


def list_mutated_sites(single_res_energy_table):
	# Get list of mutated sites
	return sorted([int(i) for i in list(single_res_energy_table[
		(single_res_energy_table['mutated'] == 1) & 
		(single_res_energy_table['site'].notna())]['site'].unique())])


def identify_influential_residues(model_info, single_res_energy_table, 
	ant_int_res, threshold=1.4):
	# Shorthand for name
	se_df = single_res_energy_table

	# Identify mutation sites
	mut_sites = list_mutated_sites(se_df)

	# Get lists of influencial residues
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
	model_info['tot_stab_mut_distal'] = \
		ut.join_list([ut.str2int(i) for i in tot_stab_mut_distal], ';') 
	model_info['tot_stab_other_close'] = \
		ut.join_list([ut.str2int(i) for i in tot_stab_other_close], ';') 
	model_info['tot_stab_other_distal'] = \
		ut.join_list([ut.str2int(i) for i in tot_stab_other_distal], ';') 
	model_info['tot_destab_mut_close'] = \
		ut.join_list([ut.str2int(i) for i in tot_destab_mut_close], ';') 
	model_info['tot_destab_other_close'] = \
		ut.join_list([ut.str2int(i) for i in tot_destab_other_close], ';')
	model_info['tot_destab_mut_distal'] = \
		ut.join_list([ut.str2int(i) for i in tot_destab_mut_distal], ';') 
	model_info['tot_destab_other_distal'] = \
		ut.join_list([ut.str2int(i) for i in tot_destab_other_distal], ';')
	model_info['int_destab_mut'] = \
		ut.join_list([ut.str2int(i) for i in int_destab_mut], ';') 
	model_info['int_destab_other'] = \
		ut.join_list([ut.str2int(i) for i in int_destab_other], ';')
	model_info['int_stab_mut'] = \
		ut.join_list([ut.str2int(i) for i in int_stab_mut], ';')
	model_info['int_stab_other'] = \
		ut.join_list([ut.str2int(i) for i in int_stab_other], ';')

	return model_info


def collect_influential_interface_res_e(res_e_tab, threshold=1.4):
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

	return method_singles


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
		ax.set_xticks(x, [str(int(i)) for i in ms.keys()])
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
	# Collect 4-member list of influential interface energies 
	# (for each method/constraint)
	method_singles = collect_influential_interface_res_e(
		res_e_tab, threshold=threshold)
	
	# Sort sites as mutant or not	
	mut_sites = list_mutated_sites(res_e_tab)			 
	site_groups = {}
	site_groups['mutated'] = \
		{k: v for k,v in method_singles.items() if k in mut_sites}
	site_groups['non-mutated'] = \
		{k: v for k,v in method_singles.items() if k not in mut_sites}

	# Plotting significant sites as bars
	if plot:
		significant_interface_site_barplots(pdbs, od, site_groups)
	
	# End here if not all models are present--no consensus possible
	if len(pdbs.missing) != 0:
		return

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

	return consensus_df


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
			print('\t' + method, constrained)
			# Get PDB names
			wt_base = '{}_{}_wt'.format(method, constrained)
			om_base = '{}_{}_om'.format(method, constrained)
			## Can't do comparison if a model is missing 
			if wt_base in pdbs.missing or om_base in pdbs.missing:
				continue
			wt_pdb = getattr(pdbs, wt_base)
			om_pdb = getattr(pdbs, om_base)

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

	# Perform consensus evaluation and plotting
	print('\tconsensus')
	consensus_dict = consensus_site_evaluation(pdbs, od, 
		single_res_es_table, plot=plot, threshold=threshold)
	print(consensus_dict)


def main(args):
	# Get PDBs
	pdbs = model_pdbs(args.antibody_name, args.rpk_cons_wt, args.rpk_cons_om, 
		args.rpk_free_wt, args.rpk_free_om, args.af2_cons_wt, args.af2_cons_om, 
		args.af2_free_wt, args.af2_free_om)

	# Make plots and tables
	make_plots_and_tables(pdbs, args.out_dir, threshold=args.threshold, 
		plot=args.plot)
			

if __name__ == '__main__':
	# Take user inputs
	args = parse_args()

	# Initializing PyRosetta
	ut.pyrosetta_init(verbose=False, extra_opts=['ignore_zero_occupancy false'])

	# Execute main
	main(args)