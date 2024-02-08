
## scalepopgen: phylogeny with TreeMix

This sub-workflow constructs the maximum likelihood tree based on allele frequencies and optionally infer migration events.

## Description of the parameters:
#### The sub-workflow is invoked by setting the argument ```treemix``` to TRUE and here are further options:
```k_snps```: the number of SNPs representing the size of the windows \
```treemix_args```: additional parameters to run the TreeMix \
```n_bootstrap```: number of bootstraps to be performed \
```set_random_seed```: an option to set random number for seed \
```n_mig```: the number of migration edges \
```n_iter```: the number of iterations for each migration edge \
```rand_k_snps```:  an option to randomize number of SNPs in window size for each iteration 

## Overview of the processed carried out in this sub-workflow: 

**1.** generating a file with samples and populations based on sample map provided by the user \
**2.** converting vcf to treemix input file format using python script \
**3.** merging converted treemix files (separated by chromosomes) into one \
**4.** running the treemix analysis without migration edges \
**5.** running the treemix analysis with bootstraps \
**6.** generating consensus tree out of the bootstrapped trees generated in the previous step \
**7.** the treemix analysis with migration edges \
**8.** identify optimal number of migration edges using the procedure as implemented in [OptM](https://academic.oup.com/biomethods/article/6/1/bpab017/6371180) package

> **Note:** If your input files are in Plink format, it will convert them into VCF before the first step.

## Validation and test-run of the sub-workflow:
For workflow validation, we have downloaded publicly available samples (see map below) with whole genome sequences from NCBI database (Alberto et al., 2018; Grossen et al., 2020; Henkel et al., 2019). We included domestic goats (*Capra hircus*) represented by various breeds from Switzerland. In addition to them, we also included Alpine ibex (*C. ibex*) and Bezoar wild goat (*C. aegagrus*). Since we need an outgroup when performing some of the analyses, we also added Urial sheep (*Ovis vignei*). We will use variants from chromosome 28 and 29 of, all together, 85 animals.

![Sample_info](https://github.com/NPogo/scalepopgen_README/assets/131758840/70fca73a-be60-4ebb-bb14-224e8efa2683)
Geographic map of the samples used for the test-run

###### Alberto et al. (2018). Convergent genomic signatures of domestication in sheep and goats. *Nature communications*, https://doi.org/10.1038/s41467-018-03206-y
###### Grossen et al. (2020). Purging of highly deleterious mutations through severe bottlenecks in Alpine ibex. *Nature communications*, https://doi.org/10.1038/s41467-020-14803-1
###### Henkel et al. (2019). Selection signatures in goats reveal copy number variants underlying breed-defining coat color phenotypes. *PLoS genetics*, https://doi.org/10.1371/journal.pgen.1008536

### 1. Required input files
The input data should be in the **VCF** or **PLINK binary** format files. 

All VCF files need to be splitted by the chromosomes and indexed with tabix. We have listed the input files in the CSV sheet (the example below) with the necessary header row. The first information in each row of the input sheet is chromosome id, next is path to the zipped VCF file and the last is path to the indexed VCF file. In our case, we inserted the links to the cloud stored data. Please note that the chromosome ID must not contain any punctuation marks.
```
chrom,vcf,vcf_idx
chr28,https://data.cyverse.org/dav-anon/iplant/home/maulik88/28_filt_samples.vcf.gz,https://data.cyverse.org/dav-anon/iplant/home/maulik88/28_filt_samples.vcf.gz.tbi
chr29,https://data.cyverse.org/dav-anon/iplant/home/maulik88/29_filt_samples.vcf.gz,https://data.cyverse.org/dav-anon/iplant/home/maulik88/29_filt_samples.vcf.gz.tbi
```
In addition to the VCF input sheet, it is also necessary to prepare a sample map file of individuals and populations. Sample map has two tab-delimited columns: in the first column are individual IDs and in the second are population IDs as demonstrated on the example below. It is also important that the name of the file ends with ".map" and unlike input sheet, there is no header.
```
SRX5250055_SRR8442974	Appenzell
SRX5250057_SRR8442972	Appenzell
SRX5250124_SRR8442905	Appenzell
SRX5250148_SRR8442881	Appenzell
SRX5250150_SRR8442879	Appenzell
SRX5250151_SRR8442878	Appenzell
SRX5250153_SRR8442876	Appenzell
SRX5250155_SRR8442874	Appenzell
SRX5250156_SRR8442873	Appenzell
SRX5250157_SRR8442872	Appenzell
340330_T1	Bezoar
340331_T1	Bezoar
340334_T1	Bezoar
340340_T1	Bezoar
340345_T1	Bezoar
340347_T1	Bezoar
340426_T1	Bezoar
470100_T1	Bezoar
470104_T1	Bezoar
470106_T1	Bezoar
...
454948_T1	Urial
ERR454947urial	Urial
SRR12396950urial	Urial
```
### 2. Optional input files

This workflow has an option to draw a geographic map with samples' origin. Please, check the README of the filtering sub-workflow, where we described the layout of the input files.


### 3. Setting the parameters
To assist the user in creating the parameter file, there is a Command-Line Interface (CLI). Please, refer to the general README for CLI installation. 
Start the CLI with:
```
python scalepopgen_cli.py
```
As we would like to create a YAML file that we do not have yet, click enter on "No".
![CLI1](https://github.com/NPogo/scalepopgen_README/assets/131758840/f2ca31c1-4364-4677-a9ae-16d49a6ef3b4)

At the beginning, we have to specify some of the general parameters, which can be found in the first tab of CLI: \
![CLI2](https://github.com/NPogo/scalepopgen_README/assets/131758840/f365a667-f1da-412b-bcea-1685bacc58b7)![CLI3](https://github.com/NPogo/scalepopgen_README/assets/131758840/3d4b11dc-1cac-473c-a418-30cc8f159370)

#### Setting the general parameters:
```input```: path to the ".csv" input sheet for the VCF files or ".p.csv" for the PLINK binary files \
```outDir```: the name of the output folder \
```sample_map```: path to the file with listed individuals and populations as addition to VCF inputs \
```color_map```: path to the file with specified colors for each population \
```outprefix```: the prefix of the output files \
```max_chrom```: maximum number of chromosomes \
```allow_extra_chrom```: set to true if the input contains chromosome ID in the form of string \
```chrom_length_map```: path to the file with listed lengths (base pairs) for each chromosome \
```fasta```: path to the reference genome fasta file that will be used for converting in case of PLINK inputs \
```outgroup```: the population ID of the outgroup \
```window_size```: window size relevant for summary statistics, Tajima's D, Pi, Fst and SweepFinder2  \
```step_size```: step size relevant for Tajima's D, Pi and Fst \

After completion of general parameters, we can move to the tabs dedicated to removal of sites and samples. Here we specify options described at the beginning of this documentation. At the end, save the parameters as YAML file:
![CLI4](https://github.com/NPogo/scalepopgen_README/assets/131758840/9bedd791-cef2-431c-b83c-b1d679a2078d)![CLI5](https://github.com/NPogo/scalepopgen_README/assets/131758840/6b765086-df4d-4ad0-8f72-c28c9611c98f)

With created YAML file of parameters, we are ready to start the workflow. Choose any container profile, we prefer mamba, and set the maximum number of processes, 10 in our case, that can be executed in parallel by each executor. From within the **scalepopgen** folder, execute the following command:
```
nextflow run scalepopgen  -params-file parameters.yml -profile mamba -qs 10
```
You can check all the other command running options with the option help :
```
nextflow run scalepopgen -help
```
If the module analyses are processed successfully, the command line output is looking like this:
```N E X T F L O W  ~  version 23.04.1
Launching `scalepopgen.nf` [shrivelled_sinoussi] DSL2 - revision: 9f9aaad1d2
executor >  local (24)
[71/4cd550] process > GENERATE_POP_COLOR_MAP (generating pop color map)                          [100%] 1 of 1 ✔
[2e/17a834] process > RUN_TREEMIX:PREPARE_POP_FILE (preparing_pop_file)                          [100%] 1 of 1 ✔
[ec/8f62b9] process > RUN_TREEMIX:VCF_TO_TREEMIX_INPUT (convert_vcf_to_treemix_input_CHR29)      [100%] 2 of 2 ✔
[18/54cea1] process > RUN_TREEMIX:MERGE_TREEMIX_INPUTS (merging_treemix_inputs)                  [100%] 1 of 1 ✔
[91/42f1bc] process > RUN_TREEMIX:RUN_TREEMIX_DEFAULT (run_treemix_default_merged_treemix_input) [100%] 1 of 1 ✔
[d0/f49d21] process > RUN_TREEMIX:RUN_TREEMIX_WITH_BOOTSTRAP (run_treemix_13895)                 [100%] 10 of 10 ✔
[6a/14795b] process > RUN_TREEMIX:RUN_CONSENSE (run_phylip_consensus)                            [100%] 1 of 1 ✔
[59/cac6d4] process > RUN_TREEMIX:ADD_MIGRATION_EDGES (adding_edge_3_2_treemix)                  [100%] 6 of 6 ✔
[ea/d5d5da] process > RUN_TREEMIX:EST_OPT_MIGRATION_EDGE (estimate_optimal_mig_edge)             [100%] 1 of 1 ✔
Completed at: 11-Aug-2023 16:00:16
Duration    : 7m 19s
CPU hours   : 0.7
Succeeded   : 24
```


### 4. Description of the output files generated by this sub-workflow:

According to different options that this tool is offering to run the TreeMix, the results will be stored in separated folders. In the main output folder you will find:

![folders](../../images/treemix_dir.png)

->**/input_files/**: input files for the program TreeMix\
->**/out_tree_default_m0/**: output files from the TreeMix analysis without bootstrapping and migration events (step 4) \
->**/out_tree_bootstrap/**: output files from the TreeMix analysis with bootstrapping (step 5) \
->**/out_tree_mig/**: output files from the TreeMix analysis with migration events (steps 7 and 8) \
->**/consensus_trees/**: output consensus trees (step 6)

The optimal number of migration events is suggested in plots **OptM_results.pdf**, which are directly in the output folder.

For our dataset the program suggested one migration event between breeds Chamois colored and Booted goats. Beside that, the topology showed a monophyletic clade containing all goat breeds from Switzerland, well separated from the Bezoar wild goat and Alpine ibex. Urials were used as the outgroup.

![plot](../../images/treemix_plot.png)

## References
Please cite the following papers if you use this sub-workflow in your study:

[1] Pickrell JK, Pritchard JK (2012) Inference of Population Splits and Mixtures from Genome-Wide Allele Frequency Data. PLOS Genetics 8(11): e1002967. https://doi.org/10.1371/journal.pgen.1002967

[2] Robert R Fitak, OptM: estimating the optimal number of migration edges on population trees using Treemix, Biology Methods and Protocols, Volume 6, Issue 1, 2021, bpab017, https://doi.org/10.1093

[3] Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319 (2017). https://doi.org/10.1038/nbt.3820


## License

MIT





