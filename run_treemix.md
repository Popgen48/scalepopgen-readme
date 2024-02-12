
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

![treemix drawio](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/3d4d1e06-8fd1-4783-8568-c0ffe3ce9baa)

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
```sample_map```: path to the file with listed individuals and populations as addition to VCF inputs (extension ".map") \
```color_map```: path to the file with specified colors for each population (extension ".map") \
```outprefix```: the prefix of the output files \
```max_chrom```: maximum number of chromosomes \
```allow_extra_chrom```: set to true if the input contains chromosome ID in the form of string \
```chrom_length_map```: path to the file with listed lengths (base pairs) for each chromosome (extension ".map") \
```fasta```: path to the reference genome fasta file that will be used for converting in case of PLINK inputs \
```outgroup```: the population ID of the outgroup \
```window_size```: window size relevant for summary statistics, Tajima's D, Pi, Fst and SweepFinder2  \
```step_size```: step size relevant for Tajima's D, Pi and Fst \

After completion of general parameters, we can move to the tabs dedicated to removal of sites and samples. Here we specify options described at the beginning of this documentation. At the end, save the parameters as YAML file:

![CLI4](https://github.com/NPogo/scalepopgen_README/assets/131758840/9bedd791-cef2-431c-b83c-b1d679a2078d)![CLI5](https://github.com/NPogo/scalepopgen_README/assets/131758840/6b765086-df4d-4ad0-8f72-c28c9611c98f)

With created YAML file of parameters, we are ready to start the workflow. Choose any container profile, we prefer mamba, and set the maximum number of processes, 10 in our case, that can be executed in parallel by each executor. From within the **scalepopgen** folder, execute the following command:
```
nextflow run scalepopgen -params-file Treemix.yml -profile mamba -qs 10
```
You can check all the other command running options with the option help :
```
nextflow run scalepopgen -help
```
If the module analyses are processed successfully, the command line output is looking like this:
```N E X T F L O W  ~  version 23.04.1
...
[84/3b2040] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:GAWK_GENERATE_COLORS (generating colors for plotting)     [100%] 1 of 1 ✔
[7c/a13ca5] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_VCF (CHR29)                                        [100%] 2 of 2 ✔
[08/88c905] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_MERGE_BED (merging_bed)                            [100%] 1 of 1 ✔
[58/086fbd] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:GAWK_UPDATE_CHROM_IDS (updating_chrom_ids)                [100%] 1 of 1 ✔
[2d/79ce64] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK_MAKE_BED (TRIAL)                                    [100%] 1 of 1 ✔
[d2/0fe0a4] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:PYTHON_CONVERT_VCF2TREEMIX (convert_vcf_to... [100%] 2 of 2 ✔
[68/7920b6] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:GAWK_MERGE_TREEMIX_INPUTS (merging_treemix... [100%] 1 of 1 ✔
[42/254376] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:TREEMIX_RUN_M0 (1_0_1234)                     [100%] 1 of 1 ✔
[6e/691b52] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:RSCRIPT_PLOT_TREE_M0 (treemix_tree)           [100%] 1 of 1 ✔
[7c/d16930] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:IMAGEMAGIK_RUN_M0 (Treemix_default)           [100%] 1 of 1 ✔
[56/c3b5c9] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:TREEMIX_RUN_M0_BOOTSTRAP (1_0_6411)           [100%] 10 of 10 ✔
[23/9354f5] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:RSCRIPT_PLOT_TREE_M0_BOOTSTRAP (treemix_tree) [100%] 10 of 10 ✔
[2a/7bea5b] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:PHYLIP_CONSENSE (consense)                    [100%] 1 of 1 ✔
[89/3ece08] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:TREEMIX_ADD_MIG (2_2_151)                     [100%] 20 of 20 ✔
[de/ec4a2c] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:RSCRIPT_PLOT_TREE_ADD_MIG (treemix_tree)      [100%] 20 of 20 ✔
[e5/5f5dab] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:RSCRIPT_OPTM (estimate_optimal_mig_edge)      [100%] 1 of 1 ✔
[35/739f09] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:IMAGEMAGIK_OPTM (OptM_results)                [100%] 1 of 1 ✔
[e4/e6be75] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:IMAGEMAGIK_RUN_ADD_MIG (TRIAL.5.2)            [100%] 20 of 20 ✔
[9c/7f6031] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_TREEMIX:IMAGEMAGIK_CONVERT_APPEND (Treemix_m2)        [100%] 2 of 2 ✔
[5a/4713d3] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:MULTIQC_GENETIC_STRUCTURE (1)                             [100%] 1 of 1 ✔
-[popgen48/scalepopgen] Pipeline completed successfully-
Completed at: 12-Feb-2024 18:00:09
Duration    : 4m 28s
CPU hours   : 1.2
Succeeded   : 98
```

### 4. Description of the output files generated by this sub-workflow:

The results are stored in different folders, the main is called **./treemix**:

![image](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/3a548f05-b870-4be6-a8f2-a1ed123a4512)

> **Note:** The output also contains a folder **./pipeline_info**, where are execution reports and used parameters.

In the directory named **./multiqc** is a link to all the plots produced by this sub-workflow. Just by opening it in any browser, we can quickly check the phylogeny.

For our dataset the program OptM suggested one migration event between Sannen breed and the clade of Appenzell and Toggenburg goats. Beside that, the topology showed a monophyletic clade containing all goat breeds from Switzerland together with the Bezoar wild goat that is separated from the Alpine ibex. Urials were used as the outgroup.

![treemix](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/ecc37564-2c00-4a52-ad35-acf57d2021f7)

## References
Please cite the following papers if you use this sub-workflow in your study:

[1] Pickrell JK, Pritchard JK (2012) Inference of Population Splits and Mixtures from Genome-Wide Allele Frequency Data. PLOS Genetics 8(11): e1002967. https://doi.org/10.1371/journal.pgen.1002967

[2] Robert R Fitak, OptM: estimating the optimal number of migration edges on population trees using Treemix, Biology Methods and Protocols, Volume 6, Issue 1, 2021, bpab017, https://doi.org/10.1093

[3] Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England), 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

[4] Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319 (2017). https://doi.org/10.1038/nbt.3820

## License

MIT





