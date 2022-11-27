# omicron_models

## Requirements

* All scripts and other neccessary files for our model generation and analysis are in the scripts folder
* Scripts utilize PyRosetta (which can be licensed at https://www.pyrosetta.org/) and several other Python modules, included in the Anaconda environment pyr.yml
* PDB files can be downloaded from the Protein Data Bank (https://www.rcsb.org/) 
* Development and use was done in Ubuntu 20.04.3 LTS and CentOS Linux release 7.9.2009

## Instalation

* Installing Anaconda (https://www.anaconda.com/) into a Linux terminal and used to duplicate our environment, in which all scripts should run correctly
* Install time may be approximately 1.5 hours

## Usage for scripts for omicron-antibody model generation and analysis

### Generation of trimmed complexes from PDB

* Trimming is done using the tabulate_res_interactions.py script with the following arguments:
	* `-p`/`--pdb`: Input PDB file
	* `-c`/`--csv`: Output name (with path) for the summary CSV file
	* `-o`/`--outpath`: Output path for the trimmed PDB
* Required input data tables. The first two are both extracted from SAbDab (http://opig.stats.ox.ac.uk/webapps/covabdab/). These tables include chain and CDR identifying information. The third is a table we generated associating PDB IDs with the names of corresponding therapeutic entities.
	* sabdab_20211201_processed.csv
	* sabdab_20211201_summary_file.csv
	* SARS-CoV-2Spike-AntibodyComplexes_12-01-2021_2.csv
* Antibody PDB files can be downloaded directly from the Protein Data Bank (https://www.rcsb.org/)
	* Some ligands cause loading errors in Rosetta. This can be prevented by deleting all non-ATOM linse from the PDB file.
* Example commands
	* `csv="trimmed_pdb_structures/pdb_complexes.csv"`
	* `outpath="trimmed_pdb_structures"`
	* `inpdb="6xc4.pdb"`
	* `python tabulate_res_interactions.py -p $inpdb -c $csv -o $outpath`
* Output:
	* A trimmed PDB file for each complex in the input PDB (a trimeric spike with three antibodies will produce three trimmed outputs)
	* A table summarizing the complexes (their source chain IDs and sets of antigen residues contacted)
	* If there are insertions/deletions in the PDB, they will be reported in an insertions.txt file for each complex
* Approximate runtime: 5 minutes
* In cases where automatic trimming script failed, manual trimming was performed by referencing SAbDab tables where available to identify relevant chains for each complex, then deleting extraneous residues and renaming chains in PyMOL.

### Performing baseline energetic minimization (Relax)

* Initial minimization is necessary to alleviate minor positional conflicts between experimental structure and Rosetta's scoring function such that the Monte Carlo algorithm has a reasonable basis for accepting/rejecting changes. (For example, a positional change of 0.1 A might make a substantial difference for Van der Waals repulsion.)
* Relax is performed using relax_new_pdb.py with the following arguments:
	* `pdb_file1`: What PDB file do you want to relax?
	* `-od`/`--out_dir`: Name an output directory for decoys (Default: current directory)
	* `-name`/`--name`: What do you want to name the relaxed PDB? (Default appends 'relaxed' to original name.)
	* `-suf`/`--name_suffix`: Add a suffix name to the output csv.
	* `-n`/`--n_decoys`: How many decoys do you want? (Default: 100)
	* `-cst`/`--constraints`: If EnzDes constraints are to be applied in addition to the default coordinate constraints, specify the file
	* `-lig`/`--ligand`: If there is a ligand, specify the params file
	* `-sug`/`--sugars`: Include sugars/glycans in the model.
	* `-sym`/`--symmetry`: If the relax should be symmetric, specify the symdef file
	* `-cwt`/`--constraint_weight`:Specify the constraints weight for coordinates and enzdes (Default: 1.0)
	* `-nocons`/`--no_constraints`:Option to not apply coordinate constraints to the pose when relaxing.
* Example commands
	* `pdb=trimmed_pdb_structures/6XC4_1.pdb`
	* `od=relaxed_pdbs`
	* `python relax_new_pdb.py $pdb -od $od -n 10`
* Output:
	* PDB files with suffix numbers from 0-n
	* A .fasc file summarizing decoy energies for all models made by a single processor
* Approximate runtime: 4 hours per decoy, depends heavily on complex size
* We generated ten decoys and selected the one with the lowest total score as the representative

### Making RRMC and RRMF models

* This type of modeling was only performed for complexes involving the RBD portion of the spike, since the repack-minimize protocol does not model insertions/deletions accurately, and there are several in the Omicron NTD
* Model generation was performed using model_rbd_mutations.py with the following arguments:
	* `-s`/`--start_model`: Input PDB file
	* `-name`/`--name`: Output name (without path) for the computed structural model
	* `-od`/`--out_dir`: Output path for the modeled PDB
	* `-n`/`--number_decoys`: How many decoys do you want?
	* `-v`/`--variant`: Which coronaviral strain to model (default: Omicron)
	* `-f`/`--fasta_file`: If modeling a new variant, provide a fasta file with a single RBD sequence.
	* `-m`/`--apply_mutations`: Option to model variant strain (instead of repacking and minimizing the relevant regions of the reference strain)
	* `-cst`/`--apply_constraints`: Option to not apply coordinate restraints when relaxing
* Example commands: 
	* `pdb=best_relaxed_pdbs/6XC4_1_relaxed.pdb`
	* `bn="$(basename $pdb .pdb)"`
	* `od=computed_omicron_models`
	* Reference RRMC
	* `jobname="${bn:0:6}"_rpk_cons_wt`
	* `python model_rbd_mutations.py -s $pdb -name ${jobname} -od $od -n 10 -cst`
	* Omicron RRMC
	* `jobname="${bn:0:6}"_rpk_cons_om`
	* `python model_rbd_mutations.py -s $pdb -name ${jobname} -od $od -n 10 -cst -m`
	* Reference RRMF
	* `jobname="${bn:0:6}"_rpk_free_wt`
	* `python model_rbd_mutations.py -s $pdb -name ${jobname} -od $od -n 10`
	* Omicron RRMF
	* `jobname="${bn:0:6}"_rpk_free_om`
	* `python model_rbd_mutations.py -s $pdb -name ${jobname} -od $od -n 10 -m`
* Output:
	* PDB files with suffix numbers from 0-n
	* A .fasc file summarizing decoy energies for all models made by a single processor
* Approximate runtime: 1 hour per decoy, depends on complex size
* We generated ten decoys for each variant and method and selected the one with the lowest total score as the representative

### Making AFRC and  AFRF models

* Starting models were generated in PyMOL by superimposing the appropriate AlphaFold2-generated RBD/NTD structures onto the corresponding domains of trimmed complex models, renumbering and renaming chains as appropriate, then outputing PDBs where the experimental domain was deleted and the AF2 domain was used instead.
* Models were then minimized using relax_new_pdb.py (described above)
* Example commands: 
	* `wtpdb=af2_pdbs/6XC4_1_af2_wt.pdb`
	* `ompdb=af2_pdbs/6XC4_1_af2_om.pdb`
	* `bn="$(basename $pdb .pdb)"`
	* `od=computed_omicron_models`
	* Reference AFRC
	* `jobname="${bn:0:6}"_af2_cons_wt`
	* `python relax_new_pdb.py $wtpdb -name ${jobname} -od $od -n 10`
	* Omicron AFRC
	* `jobname="${bn:0:6}"_af2_cons_om`
	* `python relax_new_pdb.py $ompdb -name ${jobname} -od $od -n 10`
	* Reference AFRF
	* `jobname="${bn:0:6}"_af2_free_wt`
	* `python relax_new_pdb.py $wtpdb -name ${jobname} -od $od -n 10 -nocons`
	* Omicron AFRF
	* `jobname="${bn:0:6}"_af2_free_om`
	* `python relax_new_pdb.py $ompdb -name ${jobname} -od $od -n 10 -nocons`
* Output:
	* PDB files with suffix numbers from 0-n
	* A .fasc file summarizing decoy energies for all models made by a single processor
* Approximate runtime: 4 hours per decoy, depends on complex size
* We generated ten decoys for each variant and method and selected the one with the lowest total score as the representative

### Plotting

* All eight models for each complex (when available) were collected evaluated and plots were generated by comparison. Consensus analysis was only possible when all eight were available, therefore excluding complexes not involving the RBD
* Analysis was performed using omicron_data_analysis.py with the following arguments:
	* `-ab`/`--antibody_name`: TE name for outputting the figures and tables. Should be in the form of [4-letter PDB code]\_[model number]
	* `-od`/`--out_dir`: Directory for outputting the figures and tables
	* `-rcw`/`--rpk_cons_wt`: RRMC reference CSM
	* `-rco`/`--rpk_cons_om`: RRMC omicron CSM
	* `-rfw`/`--rpk_free_wt`: RRMF reference CSM
	* `-rfo`/`--rpk_free_om`: RRMF omicron CSM
	* `-acw`/`--af2_cons_wt`: AFRC reference CSM
	* `-aco`/`--af2_cons_om`: AFRC omicron CSM
	* `-afw`/`--af2_free_wt`: AFRF reference CSM
	* `-afo`/`--af2_free_om`: AFRF omicron CSM
	* `-t`/`--threshold`: Rosetta energy threshold to use for consensus analysis (default: 1.4)
	* `-p`/`--plot`: Option to generate plots (default: True)
* Example commands: 
	* `ab="6XC4_1"`
	* `od=analysis/$ab`
	* `acw=6XC4_1_af2_cons_wt.pdb`
	* `aco=6XC4_1_af2_cons_om.pdb`
	* `afw=6XC4_1_af2_free_wt.pdb`
	* `afo=6XC4_1_af2_free_om.pdb`
	* `rcw=6XC4_1_rpk_cons_wt.pdb`
	* `rco=6XC4_1_rpk_cons_om.pdb`
	* `rfw=6XC4_1_rpk_free_wt.pdb`
	* `rfo=6XC4_1_rpk_free_om.pdb`
	* `python omicron_data_analysis.py -ab $ab -od $od -acw $acw -aco $aco -afw $afw -afo $afo -rcw $rcw -rco $rco -rfw $rfw -rfo $rfo`
* Output:
	* energy_and_changed_residue_identification.csv summarizing important interface changes
	* For each simulation condition: 
		* total_res_energies.png plot indicating per-residue total scores 
		* interface_res_energies.png plot indicating per-residue interfacial scores
		* single_res_energies.csv table indicating all per-residue energy differences between reference and Omicron models
		* mutated_interface_res_energies.png plot showing interface energies at mutated sites with consensus-identified stabilization/destabilization
		* non-mutated_interface_res_energies.png showing interface energies at non-mutated sites with consensus-identified stabilization/destabilization 
	* If all 8 models are present:
		* consensus_table.csv table identifying consensus stabilization/destabilization predictions at each mutated site
* Approximate runtime: 5 minutes