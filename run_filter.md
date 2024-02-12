## scalepopgen: filtering and basic statistics
The workflow incorporates [plink2](https://www.cog-genomics.org/plink/2.0/) (v 2.00a3.7) and [VCFtools](https://vcftools.github.io/index.html) (v 0.1.16) to carry out sample and site filtering. The filtering process depends on the format of the input files. If the inputs are **PLINK binary** files, then sample and site filtering will be done with plink2, while for the **VCF** inputs, we used plink2 and VCFtools. 

Options for sample filtering:
 - removing samples according to the threshold of calculated kinship coefficient  
 - removing samples according to the amount of missing genotypes
 - removing user-defined samples

Options for site filtering:
 - filtering based on minor allele frequency threshold
 - filtering according to the Hardy-Weinberg equilibrium exact test
 - filtering SNPs according to the amount of missing information 
 - removing user-defined SNPs
 - filtering SNPs based on min and max average depths
 - filtering SNPs based on quality score threshold
> **Note:** For the VCF input all listed options of site filtering are available, while in the case of PLINK binary files the workflow considers only the first four options.
> 
## Description of the parameters:
#### When the argument ```apply_indi_filters``` is set to TRUE, there are further options:
```king_cutoff```: samples with relationship coefficient greater than this will be removed \
```mind```: samples with missing genotypes greater than the threshold selected here will be removed \
```rem_indi```:  the name of text file containing individuals to be removed (extension ".map")

#### If the argument ```apply_snp_filters``` is set to TRUE, further options can be specified:
```remove_snps```: the name of text file containing names of SNPs that will be removed (extension ".map") \
```maf```: SNPs with minor allele frequencies less than the threshold specified here will be removed \
```min_meanDP```: SNPs with average depth (across the samples) less than the threshold specified here will be removed (only for the VCF inputs) \
```max_meanDP```: SNPs with average depth (across the samples) greater than the threshold specified here will be removed (only for the VCF inputs) \
```hwe```: SNPs with HWE p-values of less than the threshold specified here will be removed \
```max_missing```: SNPs, for which the proportion of missing genotypes exceeded this threshold, will be removed \
```minQ```: SNPs with base quality less than specified here will be removed (only for the VCF inputs) \
```indiv_summary```: if set to true it will calculate sample-based summary statistics, after individual- and site-based filtering

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
This module enables to remove samples of you choice with option ```rem_indi```. You have to prepare a space/tab-delimited text file with population IDs in the first column and sample IDs in the second column as stated for option ```--remove``` in [PLINK documentation](https://www.cog-genomics.org/plink/2.0/filter#sample):
```
Grigia SRX7715567_SRR11076301
Grigia SRX7715568_SRR11076300
Grigia SRX7715579_SRR11076289
```
Similarly, you can also prepare a file with the SNPs that should be removed with option ```remove_snps```. The appearance of the file depends on the input formats. For the VCF input, we used option ```--exclude-positions``` of [VCFtools](https://vcftools.sourceforge.net/man_latest.html), which expect a text file with a tab-separated chromosome and position of the SNP in each line:
```
NC_030835.1	1460
NC_030835.1	1498
NC_030835.1	2140
NC_030835.1	2158
NC_030835.1	2173
NC_030835.1	44663424
NC_030835.1	44663455
NC_030835.1	44664291
NC_030835.1	44664392
NC_030835.1	44664464
NC_030836.1	1510
NC_030836.1	1525
NC_030836.1	1558
NC_030836.1	1593
NC_030836.1	1595
NC_030836.1	51331074
NC_030836.1	51331096
NC_030836.1	51331185
NC_030836.1	51331237
NC_030836.1	51331522
```

There is also an option to draw a geographic map with sample origin. For that we need to provide two files. In the first one we write down population ID in the first column and comma separated latitude and longitude in the second column:
```Bezoar	32.662864436650814,51.64853259116807
Urial	34.66031157,53.49391737
AlpineIbex	46.48952713,9.832698605
ChamoisColored	46.620927266181674,7.345747305114329
Appenzell	47.33229709563813,9.401363933224248
Booted	47.426361052956736,9.384330852599533
Peacock	46.321661051197026,8.804738507288173
Toggenburg	47.358160245764715,9.01070577172017
Grigia	46.24935612558498,8.700996940189137
Saanen	46.9570926960748,8.205509946726016
```
In the second file, we specify the hex codes of colors (2. column) that will represent each population (1. column):
```AlpineIbex	#008000
Appenzell	#ff5733
Booted	#0000FF
ChamoisColored	#d6b919
Grigia	#aee716
Peacock	#16e7cc
Saanen	#75baf3
Urial	#A52A2A
Toggenburg	#da4eed
Bezoar	#FFA500
```
The last file is not obligatory as the tool can choose random colors, while the first one with coordinates is necessary for the map plotting.

Please, pay attention that files for options ```rem_indi```, ```remove_snps``` and ```color_map``` have the extension ".map".

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
nextflow run scalepopgen -params-file Filtering.yml -profile mamba -qs 10
```
You can check all the other command running options with the option help :
```
nextflow run scalepopgen -help
```
After the filtering is done successfully, the command line output is looking like this:
```
N E X T F L O W  ~  version 23.04.1
...
[f3/e1cc4b] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:REMOVE_SAMPLE_LIST (combining_indiv_summary)   [100%] 1 of 1 ✔
[4d/406b14] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:VCFTOOLS_CONCAT (concate_vcf)                  [100%] 1 of 1 ✔
[6e/e6cd95] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:PLINK2_FILTER_SAMPLES (filter_indi_TRIAL)      [100%] 1 of 1 ✔
[b2/c88b3b] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:GAWK_EXTRACT_SAMPLEID (combining_indiv_summ... [100%] 1 of 1 ✔
[5e/7c4ed3] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:VCFTOOLS_KEEP (keep_indi_CHR29)                [100%] 2 of 2 ✔
[1c/2e0503] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:GAWK_PREPARE_NEW_MAP (preparing_new_map)       [100%] 1 of 1 ✔
[fd/a30a70] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:VCFTOOLS_FILTER_SITES (filter_sites_CHR29)     [100%] 2 of 2 ✔
[8a/0ea7b6] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:FILTER_VCF:TABIX (CHR29)                                  [100%] 2 of 2 ✔
[0d/310a0a] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:GAWK_GENERATE_COLORS (generating colors for plotting)     [100%] 1 of 1 ✔
-[popgen48/scalepopgen] Pipeline completed successfully-
Completed at: 12-Feb-2024 13:58:45
Duration    : 8m 10s
CPU hours   : 1.3
Succeeded   : 10
```

### 4. Description of the outputs generated by this sub-workflow:
The scheme of output files differs according to the input formats. In the case of VCF input it will look like this:

![image](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/08f1c5f1-43a3-4799-b707-fe0272deef8b)

The workflow will start with sample-based filtering. After that it will take the sample-filtered files and perform SNP-based filtering. We can find the final filtered files in the folder **./snp_filtering/vcftools/sites_filtes** and their indexed files in **./tabix**. Based on the remaining samples this tool will also prepare the updated sample map (_* *_new_sample_pop.map*_) in the folder **./gawk/prepare_new_map**.

> **Note:** The output also contains a folder **./pipeline_info**, where are execution reports and used parameters.

### 5. Generating the geographical map without running the workflow
For generating the geographic map of samples (without re-running the workflow), one can either use **geo_map_file.txt** generated by the workflow or prepare a tab-delimited file in the similar format. Run the python script with the following command: 
```
python3 plot_sample_info.py geo_map_file.txt plot_sample_on_map.yml tiles_info.yml
```
The python script is located in the bin folder of scalepopgen. The yaml files are located in the folder "/parameters/plots/". The parameters of the yaml files are described below:
```
tile: the name of the tiles to be used for plotting. Default: "world_gray_canvas". Note that the attributes of this tile should be present in tiles_info.yml 
marker_prop: in case you want to plot the area of the circle proportional to the sample size, specify the value here 
const_radius : radius of the circle to be plotted 
zoom_level : zoom level of the map to be plotted 
overlap_cordi: if the coordinates of the samples are overlappin, the value should be boolean, True or False
shift_cordi: if the coordinates are overlapping, shift the coordinates by this value. 
show_label: whether or not the label should be shown in the map, the value should be boolean, True or False 
label_loc : if the label to be displayed near the circle set to be true, specify the location from the circle 
label_size: size of the label to be plotted 
show_legend: whether or not to plot the legend 
display_sample_size: whether or not to display the sample size along with the pop label in the legend, the value should be boolean 
display_popup: whether or not to pop-up the information, the value should be boolean 
output_prefix: prefix of the output html file
```
## References
Please cite the following papers if you use this sub-workflow in your study:

[1] Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R. E., Lunter, G., Marth, G. T., Sherry, S. T., McVean, G., Durbin, R., & 1000 Genomes Project Analysis Group (2011). The variant call format and VCFtools. _Bioinformatics (Oxford, England)_, _27_(15), 2156–2158. https://doi.org/10.1093/bioinformatics/btr330

[2] Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., Maller, J., Sklar, P., de Bakker, P. I., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. _American journal of human genetics_, _81_(3), 559–575. https://doi.org/10.1086/519795

[3] Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319 (2017). https://doi.org/10.1038/nbt.3820

## License

MIT

