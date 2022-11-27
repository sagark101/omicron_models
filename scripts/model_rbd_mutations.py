import argparse
import joey_utils as ut
import pyrosetta as pr

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--start_model", required=True, type=str)
	parser.add_argument("-name", "--name", required=True, type=str)
	parser.add_argument("-od", "--out_dir", type=str)
	parser.add_argument("-n", "--number_decoys", required=True, type=int)
	parser.add_argument("-v", "--variant", default='om', choices=[
		'de', 's6p', 'om', 'osu', 'ba2', 'ba_2_12', 'ba4', 'new'])
	parser.add_argument("-f", "--fasta_file", type=str, help="If modeling a \
		new variant, provide a fasta file with a single RBD sequence.")
	parser.add_argument("-m", "--apply_mutations", action="store_true")
	parser.add_argument("-cst", "--apply_constraints", action="store_true")
	args = parser.parse_args()

	if args.variant == 'new' and not args.fasta_file:
		raise argparse.ArgumentTypeError('If modeling a new variant, provide \
			a fasta file with a single RBD sequence.')
	if args.variant != 'new' and args.fasta_file:
		raise argparse.ArgumentTypeError('If modeling an existing variant, no \
			fasta file will be used for model generation.')

	return args


class rbd_seqs():
	def __init__(self):
		self.wt = {319: 'R', 320: 'V', 321: 'Q', 322: 'P', 323: 'T', 324: 'E', 
			325: 'S', 326: 'I', 327: 'V', 328: 'R', 329: 'F', 330: 'P', 
			331: 'N', 332: 'I', 333: 'T', 334: 'N', 335: 'L', 336: 'C', 
			337: 'P', 338: 'F', 339: 'G', 340: 'E', 341: 'V', 342: 'F', 
			343: 'N', 344: 'A', 345: 'T', 346: 'R', 347: 'F', 348: 'A', 
			349: 'S', 350: 'V', 351: 'Y', 352: 'A', 353: 'W', 354: 'N', 
			355: 'R', 356: 'K', 357: 'R', 358: 'I', 359: 'S', 360: 'N', 
			361: 'C', 362: 'V', 363: 'A', 364: 'D', 365: 'Y', 366: 'S', 
			367: 'V', 368: 'L', 369: 'Y', 370: 'N', 371: 'S', 372: 'A', 
			373: 'S', 374: 'F', 375: 'S', 376: 'T', 377: 'F', 378: 'K', 
			379: 'C', 380: 'Y', 381: 'G', 382: 'V', 383: 'S', 384: 'P', 
			385: 'T', 386: 'K', 387: 'L', 388: 'N', 389: 'D', 390: 'L', 
			391: 'C', 392: 'F', 393: 'T', 394: 'N', 395: 'V', 396: 'Y', 
			397: 'A', 398: 'D', 399: 'S', 400: 'F', 401: 'V', 402: 'I', 
			403: 'R', 404: 'G', 405: 'D', 406: 'E', 407: 'V', 408: 'R', 
			409: 'Q', 410: 'I', 411: 'A', 412: 'P', 413: 'G', 414: 'Q', 
			415: 'T', 416: 'G', 417: 'K', 418: 'I', 419: 'A', 420: 'D', 
			421: 'Y', 422: 'N', 423: 'Y', 424: 'K', 425: 'L', 426: 'P', 
			427: 'D', 428: 'D', 429: 'F', 430: 'T', 431: 'G', 432: 'C', 
			433: 'V', 434: 'I', 435: 'A', 436: 'W', 437: 'N', 438: 'S', 
			439: 'N', 440: 'N', 441: 'L', 442: 'D', 443: 'S', 444: 'K', 
			445: 'V', 446: 'G', 447: 'G', 448: 'N', 449: 'Y', 450: 'N', 
			451: 'Y', 452: 'L', 453: 'Y', 454: 'R', 455: 'L', 456: 'F', 
			457: 'R', 458: 'K', 459: 'S', 460: 'N', 461: 'L', 462: 'K', 
			463: 'P', 464: 'F', 465: 'E', 466: 'R', 467: 'D', 468: 'I', 
			469: 'S', 470: 'T', 471: 'E', 472: 'I', 473: 'Y', 474: 'Q', 
			475: 'A', 476: 'G', 477: 'S', 478: 'T', 479: 'P', 480: 'C', 
			481: 'N', 482: 'G', 483: 'V', 484: 'E', 485: 'G', 486: 'F', 
			487: 'N', 488: 'C', 489: 'Y', 490: 'F', 491: 'P', 492: 'L', 
			493: 'Q', 494: 'S', 495: 'Y', 496: 'G', 497: 'F', 498: 'Q', 
			499: 'P', 500: 'T', 501: 'N', 502: 'G', 503: 'V', 504: 'G', 
			505: 'Y', 506: 'Q', 507: 'P', 508: 'Y', 509: 'R', 510: 'V', 
			511: 'V', 512: 'V', 513: 'L', 514: 'S', 515: 'F', 516: 'E', 
			517: 'L', 518: 'L', 519: 'H', 520: 'A', 521: 'P', 522: 'A', 
			523: 'T', 524: 'V', 525: 'C', 526: 'G', 527: 'P', 528: 'K', 
			529: 'K', 530: 'S', 531: 'T', 532: 'N', 533: 'L', 534: 'V', 
			535: 'K', 536: 'N', 537: 'K', 538: 'C', 539: 'V', 540: 'N', 
			541: 'F'}

		self.delta = {452: 'R', 478: 'K'}

		self.s6p = {417: 'N', 484: 'K', 501: 'Y'}

		self.omicron = {339: 'D', 371: 'L', 373: 'P', 375: 'F', 417: 'N', 
			440: 'K',  446: 'S', 477: 'N', 478: 'K', 484: 'A', 493: 'R', 
			496: 'S', 498: 'R', 501: 'Y', 505: 'H'}

		self.osu_omicron = {339: 'D', 346: 'K', 371: 'L', 373: 'P', 375: 'F', 
			417: 'N', 440: 'K', 446: 'S', 477: 'N', 478: 'K', 484: 'A', 
			493: 'R', 496: 'S', 498: 'R', 501: 'Y', 505: 'H'}

		self.omicron_ba2 = {339: 'D', 371: 'F', 373: 'P', 375: 'F', 376: 'A', 
			405: 'N', 408: 'S', 417: 'N', 440: 'K', 477: 'N', 478: 'K', 
			484: 'A', 493: 'R', 498: 'R', 501: 'Y', 505: 'H'}

		self.ba_2_12 = {339: 'D', 371: 'F', 373: 'P', 375: 'F', 376: 'A', 
			405: 'N', 408: 'S', 417: 'N', 440: 'K', 452: 'Q', 477: 'N', 
			478: 'K', 484: 'A', 493: 'R', 498: 'R', 501: 'Y', 505: 'H'}

		self.ba_4 = {339: 'D', 371: 'F', 373: 'P', 375: 'F', 376: 'A', 
			405: 'N', 408: 'S', 417: 'N', 440: 'K', 452: 'R', 477: 'N', 
			478: 'K', 484: 'A', 486: 'V', 498: 'R', 501: 'Y', 505: 'H'}

	def get_variant(self, variant, mutated=True):
		available_variants = ['de', 's6p', 'om', 'osu', 'ba2', 'ba_2_12', 'ba4', 
			'new']
		variant = variant.lower()
		if variant not in available_variants:
			raise ValueError(
				'Variant oprions are: {}'.format(', '.join(available_variants)))

		if variant == 'de':
			variant_dict = self.delta
		if variant == 's6p':
			variant_dict = self.s6p
		if variant == 'om':
			variant_dict = self.omicron
		if variant == 'osu':
			variant_dict = self.osu_omicron
		if variant == 'ba2':
			variant_dict = self.omicron_ba2
		if variant == 'ba_2_12':
			variant_dict = self.ba_2_12
		if variant == 'ba4':
			variant_dict = self.ba_4
		if variant == 'new':
			variant_dict = self.new_variant

		if mutated:
			return variant_dict
		else: 
			return {k: v for k, v in self.wt.items() if k in variant_dict}

	def add_new_variant(self, fasta_file):
		# Read in fasta sequence
		new_variant_sequence = ut.parse_fastafile(fasta_file)
		if len(new_variant_sequence) > 1:
			raise RuntimeError("""Fastsa file includes multiple sequences. Use a
				single-sequence fasta file input.""")
		new_variant_sequence = new_variant_sequence[0]

		# Align sequence with WT and confirm equal sequence length
		ref_seq = ut.join_list(self.wt.values(), '')
		alignment = ut.global_align(ref_seq, new_variant_sequence.seq)[0]
		alignment_table = ut.tabulate_sequence_alignment(alignment)
		if len(alignment_table[alignment_table['aligned'] == 0]) > 0:
			raise RuntimeError("""Sequence includes insertions and/or deletions. 
				This method is only suitable for modeling sequences of equal 
				length.""")

		# Renumber alignment table and make a dict of substitutions
		alignment_table['s1_position'] = alignment_table['s1_position'] + 318
		mut_dict = {}
		for i, r in alignment_table[alignment_table['identity'] == 0].iterrows():
			mut_dict[r['s1_position']] = r['s2_sequence']

		# Add new variant sequences dict
		self.new_variant = mut_dict


def convert_variant_subs(ptab, sequence_changes, rbd_chain='A'):
	site_changes = {}
	for pdb_num, aa in sequence_changes.items():
		# Identify table row with given PDB number in RBD (chain A)
		row_with_pdb_num = ptab[(ptab['pdb_chain'] == rbd_chain) & 
			(ptab['pdb_number'] == pdb_num)]

		# Case where residue is absent from PDB -- skip
		if len(row_with_pdb_num) == 0:
			continue

		# Convert to pose numbering
		pose_num = int(row_with_pdb_num['pose_number'])
		site_changes[pose_num] = aa

	return site_changes


def main(args):
	print('\npdb:\n', args.start_model, '\n')
	# Loading pose and tabulating
	if args.apply_constraints:
		pose = ut.load_pose(args.start_model, coord_cst=True)
	else:
		pose = ut.load_pose(args.start_model)
	ptab = ut.tabulate_pose_residues(pose)

	# Get sequence changes
	rbd_seq_obj = rbd_seqs()
	if args.variant == 'new':
		rbd_seq_obj.add_new_variant(args.fasta_file)
	seq_changes = rbd_seq_obj.get_variant(args.variant, args.apply_mutations)

	# Convert substitutions dict to pose numbers for this PDB
	# Handles multi-chain models
	substitutions_dict = {}
	antigen_chain_selections = []
	for chain in ['A', 'B']:
		if chain in list(ptab['pdb_chain']):
			substitutions_dict.update(convert_variant_subs(ptab, seq_changes, 
				chain))
			antigen_chain_selections.append(ut.chain_selector(chain))
	print('substitutions_dict:', substitutions_dict, '\n')

	# Make residue selections
	mutation_sites = ut.index_selector(list(substitutions_dict.keys()))
	mutation_shell = ut.close_interface_with_selection_selector(mutation_sites)
	antigen_chains = ut.selector_union(*antigen_chain_selections)
	print('mutation_sites:\n', ut.selector_to_list(pose, mutation_sites, False), '\n')
	print('mutation_shell:\n', ut.selector_to_list(pose, mutation_shell, False), '\n')
	print('antigen_chains:\n', ut.selector_to_list(pose, antigen_chains, False), '\n')

	# Set up FastRelax with score function and task factory
	## Score function
	sf = ut.get_sf(constrain=1.0)
	## Task Factory
	tf = ut.make_task_factory(res_changes=substitutions_dict, 
		repack_selection=mutation_shell)
	## Move map with everything mobile
	mm = ut.make_move_map(True, True, True)
	## Mover
	fr = ut.fast_relax_mover(sf, tf, mm)

	# Set decoy name
	decoy_name = ut.output_file_name(args.name, path=args.out_dir)
	print('decoy_name:', decoy_name)

	# Relax the pose
	jd = ut.pyjobdistributor(decoy_name, args.number_decoys, sf)
	current_decoy = 1
	while not jd.job_complete:
		print('Decoy {}/{}\n'.format(current_decoy, args.number_decoys))
		pp = pr.Pose(pose)
		fr.apply(pp)
		ut.interaction_energy(pp, sf, antigen_chains, 
			ut.not_selector(antigen_chains), apply_sm=True)
		ddg = ut.ddg_of_selections(pp, sf, antigen_chains, 
			ut.not_selector(antigen_chains))
		jd.additional_decoy_info = ddg
		current_decoy += 1
		jd.output_decoy(pp) 


if __name__ == '__main__':
	# Take user inputs
	args = parse_args()

	# Initializing PyRosetta
	ut.pyrosetta_init(verbose=False)

	# Execute main
	main(args)