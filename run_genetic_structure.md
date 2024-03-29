﻿## scalepopgen: explore genetic structure
This sub-workflow is used for analyses with [Admixture](https://dalexander.github.io/admixture/) , [SmartPCA ](https://github.com/chrchang/eigensoft/tree/master/POPGEN), and clustering based on [IBS and Fst distances](https://www.cog-genomics.org/plink/2.0/). 

## Description of the parameters:
#### Setting the argument ```genetic_structure``` to TRUE will run this sub-workflow that has further options:
```rem_indi_structure```: the path to the file with listed samples that are removed in **all analyses** of this sub-workflow \
```ld_filt```: an option to use LD-based filtering according to parameters specified below \
```ld_window_size```: a window size in variant count or kilo bases \
```ld_step_size```: number of variants to shift the window at the end of each step\
```r2_threshold```: squared correlation threshold; at each step, only pairs of variants with r2 greater than the threshold are recognized \
```smartpca```: an option to run SmartPCA \
```smartpca_param```: the path to the file with additional parameters for SmartPCA\
```pca_plot_yml```: the path to the yml file containing the parameters to plot interactive PCA results\
```marker_map```: the path to the text file with specified marker shapes for each population (extension ".map") \
```chrom_map```: the path to the text file with specified new chromosome IDs (extension ".map") \
```admixture```: an option to run Admixture \
```start_k```: starting number of clusters for Admixture analysis \
```end_k```: maximal number of clusters for Admixture analysis \
```admixture_args```: additional arguments for Admixture analysis (like "--cv=") \
```admixture_colors```: the path to the text file with custom colors for Admixture plot \
```admixture_plot_pop_order```: the path to the text file with population IDs in order they should be plotted on Admixture plot \
```admixture_plot_yml```: the path to the yml file containing the parameters to plot interactive Admixture results \
```pairwise_global_fst```: an option to calculate pairwise Fst distances between populations \
```fst_plot_yml```: the path to the yml file containing parameters for plotting interactive NJ tree with Fst distances \
```ibs_dist```: an option to calculate IBS distances between individuals \
```ibs_plot_yml```: the path to the yml file containing parameters for plotting interactive NJ tree with IBS distances

## Overview of the processed carried out in this sub-workflow: 

![genetic_structure drawio](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/65c3e9ab-2cb2-4ca6-bd67-54eec0c38db7)

## Validation and test-run of the sub-workflow:
For workflow validation, we have downloaded publicly available samples (see map below) with whole genome sequences from NCBI database (Alberto et al., 2018; Grossen et al., 2020; Henkel et al., 2019). We included domestic goats (*Capra hircus*) represented by various breeds from Switzerland. In addition to them, we also included Alpine ibex (*C. ibex*) and Bezoar wild goat (*C. aegagrus*). Since we need an outgroup when performing some of the analyses, we also added Urial sheep (*Ovis vignei*). We will use variants from chromosome 28 and 29 of, all together, 85 animals.

![Sample_info](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/1a50f0ef-6673-41a2-877b-f7f9e9031dc1)
Geographic map of the samples used for the test-run

###### Alberto et al. (2018). Convergent genomic signatures of domestication in sheep and goats. *Nature communications*, https://doi.org/10.1038/s41467-018-03206-y
###### Grossen et al. (2020). Purging of highly deleterious mutations through severe bottlenecks in Alpine ibex. *Nature communications*, https://doi.org/10.1038/s41467-020-14803-1
###### Henkel et al. (2019). Selection signatures in goats reveal copy number variants underlying breed-defining coat color phenotypes. *PLoS genetics*, https://doi.org/10.1371/journal.pgen.1008536

### 1. Required input files
The input data should be in the **VCF** or **PLINK binary** format files. 

All VCF files need to be splitted by the chromosomes and indexed with tabix. We have listed the input files in the CSV sheet (the example below) with the necessary header row. The first information in each row of the input sheet is chromosome id, next is path to the zipped VCF file and the last is path to the indexed VCF file. In our case, we inserted the links to the cloud stored data. Please note that the chromosome ID must not contain any punctuation marks.
```
chrom,vcf,vcf_idx
NC_030835.1,https://data.cyverse.org/dav-anon/iplant/home/maulik88/28_filt_samples.vcf.gz,https://data.cyverse.org/dav-anon/iplant/home/maulik88/28_filt_samples.vcf.gz.tbi
NC_030836.1,https://data.cyverse.org/dav-anon/iplant/home/maulik88/29_filt_samples.vcf.gz,https://data.cyverse.org/dav-anon/iplant/home/maulik88/29_filt_samples.vcf.gz.tbi
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
This sub-workflow provide the option (```rem_indi_structure```) to remove desired samples from all the analyzes within it. For example, during the filtering sub-workflow, the samples of Grigia goat breed were removed. For the analyses of genetic structure, we would like to exclude samples that belong to the outgroup. A space/tab-delimited text file with population IDs in the first column and sample IDs in the second column should be provided:
```
Urial 454948_T1
Urial ERR454947urial
Urial SRR12396950urial
```
The desired colors for each "*K*" of Admixture analysis (```admixture_colors```) need to be specified in a text file with one column:
```
#0000FF
#d6b919
#16e7cc
#008000
#ff5733
#75baf3
#da4eed
#FFA500
```
Additionally, the Admixture plot can be made with a certain order of populations (```admixture_plot_pop_order```). For that, we need to prepare a text file with ordered population IDs in one column :
```
Appenzell
Booted
ChamoisColored
Peacock
Saanen
Toggenburg
Bezoar
AlpineIbex
```
Similarly, mark shapes for each population (```marker_map```) can also be provided for the PCA plotting. Available shapes are listed in **./extra/markershapes.txt**. Again, population IDs are in the first column and specified shapes in the second:
```
AlpineIbex	square_default
Appenzell	circle_default
Booted	circle_default
ChamoisColored	circle_default
Peacock	circle_default
Saanen	circle_default
Toggenburg	circle_default
Bezoar	triangle_default
```

In the case of PCA, you can also provide your own file with optional parameters (```smartpca_param```). To make one, please use the [instructions of the EIGENSOFT software](https://github.com/chrchang/eigensoft/tree/master/POPGEN).

### 3. Setting the parameters
To assist the user in creating the parameter file, there is a Command-Line Interface (CLI). Please, refer to the general README for CLI installation. 
Start the CLI with:
```
python scalepopgen_cli.py
```
As we would like to create a YAML file that we do not have yet, click enter on "No".

![CLI1](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/8e1dbf6f-d6a8-4f74-99d4-f98e96a39640)

At the beginning, we have to specify some of the general parameters, which can be found in the first tab of CLI:

![CLI2](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/41c162a6-25e1-4db3-b951-32e10c1fda1b)
![CLI3](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/60110de7-19af-4d33-9cd3-f6202180c0c9)

#### Setting the general parameters:
```input```: path to the ".csv" input sheet for the VCF files or ".p.csv" for the PLINK binary files \
```outDir```: the name of the output folder \
```sample_map```: path to the file with listed individuals and populations as addition to VCF inputs (extension ".map") \
```color_map```: path to the file with specified colors for each population (extension ".map") \
```outprefix```: the prefix of the output files \
```max_chrom```: maximum number of chromosomes \
```allow_extra_chrom```: set to true if the input contains chromosome ID in the form of string \
```chrom_length_map```: path to the file with listed lengths (base pairs) for each chromosome (extension ".map") \
```chrom_id_map```: path to the file with listed chromosome IDs in the form of string (first column) and number (second column) (extension ".map") \
```fasta```: path to the reference genome fasta file that will be used for converting in case of PLINK inputs \
```outgroup```: the population ID of the outgroup \
```window_size```: window size relevant for summary statistics, Tajima's D, Pi, Fst and SweepFinder2  \
```step_size```: step size relevant for Tajima's D, Pi and Fst \
```indiv_summary```: if set to true it will calculate sample-based summary statistics, after individual- and site-based filtering

After completion of general parameters, we can move to the tab dedicated to explore genetic structure. Here we specify options described at the beginning of this documentation. At the end, save the parameters as YAML file:

![CLI4](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/4bc056b6-e1d0-458e-bf11-5f6b6d57a79e)
![CLI5](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/c817da5f-d23d-4617-af16-e5400980a56a)


With created YAML file of parameters, we are ready to start the workflow. Choose any container profile, we prefer mamba, and set the maximum number of processes, 10 in our case, that can be executed in parallel by each executor. From within the **scalepopgen** folder, execute the following command:
```
nextflow run scalepopgen  -params-file Gen_structure.yml -profile mamba -qs 10
```
You can check all the other command running options with the option help :
```
nextflow run scalepopgen -help
```
If the module analyses are processed successfully, the command line output is looking like this:
```
N E X T F L O W  ~  version 23.04.1
...
[f0/85291b] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:GAWK_GENERATE_COLORS (generating colors for plotting)     [100%] 1 of 1 ✔
[94/37a987] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_VCF (CHR28)                                        [100%] 2 of 2 ✔
[ef/9f70e4] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_MERGE_BED (merging_bed)                            [100%] 1 of 1 ✔
[e8/cbe397] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_REMOVE_CUSTOM_INDI (remove_indi_pca_TRIAL_rem_i... [100%] 1 of 1 ✔
[91/3447c7] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK2_INDEP_PAIRWISE (ld_filtering_TRIAL_rem_indi)       [100%] 1 of 1 ✔
[ea/92afc6] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:GAWK_UPDATE_CHROM_IDS (updating_chrom_ids)                [100%] 1 of 1 ✔
[ab/045070] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:PLINK_MAKE_BED (TRIAL_rem_indi_ld_filtered)               [100%] 1 of 1 ✔
[89/5fdcc8] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:PLINK2_EXPORT_PED (merging_bed)                   [100%] 1 of 1 ✔
[a7/03be2a] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:PYTHON_CREATE_EIGENSTRAT_PAR (create smartpca ... [100%] 1 of 1 ✔
[43/742896] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:GAWK_MODIFY_PHENO_COL_PED (preparing_new_map)     [100%] 1 of 1 ✔
[1b/f16b34] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:EIGENSOFT_CONVERTF (TRIAL_rem_indi_ld_filtered... [100%] 1 of 1 ✔
[ac/b1bd83] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:PYTHON_CREATE_SMARTPCA_PAR (create smartpca par)  [100%] 1 of 1 ✔
[95/ef31ea] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:EIGENSOFT_SMARTPCA (TRIAL_rem_indi_ld_filtered... [100%] 1 of 1 ✔
[eb/3b8ebc] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:RUN_PCA:PYTHON_PLOT_PCA (plot_interactive_pca)            [100%] 1 of 1 ✔
[92/713c8f] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:CALC_FST:GAWK_MAKE_CLUSTER_FILE (making_cluster_file)     [100%] 1 of 1 ✔
[85/6bd339] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:CALC_FST:PLINK2_CALC_PAIRWISE_FST (pairwise_fst_null)     [100%] 1 of 1 ✔
[18/a9ebe6] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:CALC_FST:PYTHON_PLOT_PAIRWISE_FST (plot_pairwise_fst_n... [100%] 1 of 1 ✔
[93/0dbf80] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:CALC_1_MIN_IBS_DIST:PLINK_CALC_1_MIN_IBS_DIST (1_min_i... [100%] 1 of 1 ✔
[a7/2f845d] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:CALC_1_MIN_IBS_DIST:PYTHON_PLOT_1_MIN_IBS_DIST (plot_1... [100%] 1 of 1 ✔
[e2/cb686d] process > POPGEN48_SCALEPOPGEN:SCALEPOPGEN:MULTIQC_GENETIC_STRUCTURE (1)                             [100%] 1 of 1 ✔
-[popgen48/scalepopgen] Pipeline completed successfully-
Completed at: 12-Feb-2024 17:01:57
Duration    : 1m 21s
CPU hours   : (a few seconds)
Succeeded   : 21
```

### 4. Description of the output files generated by this sub-workflow:
The results are stored in the different folders: 

![image](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/15d0d09c-2cd8-4cc4-802b-f8b78e9bc559)

The workflow produce separace folder according to different analyses:  \
 -the eigenvalues and eigenvectors of the PCA are stored in the folder **./pca/eigensoft/smartpca/**  \
 -the q-matrixes of admixture analysis are in the folder **./admixture/ADMIXTURE/**  \
 -the calculated IBS-1 distances can be found in the folder **./ibs_clustering/plink/calc_1_mins_ibs_dist/**  \
 -the calculated pairwise FST distances are in the folder **./fst_clustering/plink2/**

> **Note:** The output also contains a folder **./pipeline_info**, where are execution reports and used parameters.

In the directory named **./multiqc/** is a link to all the plots produced by this sub-workflow. Just by opening it in any browser, we can quickly check the results and gain insight into the genetic structure of the samples. As an example, let's take a look at the PCA results. Our samples cluster into three groups. In the first one we can found all breeds of domestic goats from Switzerland. In another cluster are Alpine ibexes and in the third one are Bezoar wild goats.

![PCA](https://github.com/Popgen48/scalepopgen-readme/assets/131758840/30d6d0d7-812b-4fb3-add0-dd3cc98a92b2)
Principal Component Analysis

### 5. Generating the interactive plots without running the workflow
For generating the interactive PCA plot only (without re-running the workflow), one can use the following command: 
```
python3 plot_interactive_pca.py <eigenvect_file> <eigenval_file> pop_markershape_col.txt pca.yml <output_prefix>
```
The python script (**/bin/plot_interactive_pca.py**) and the yaml (**/extra/plots/pca.yml**) file are located in the scalepopgen's folder. The **<eigenvect_file>** and the **<eigenval_file>** are located in the respective output folder of smartpca (**./pca/eigensoft/smartpca/**). The **pop_markershape_col.txt** is located in the output folder **./pca/python/plot/pca/** or one can also create this tab-delimited file with this format: the first column is pop_id, the second column is shape_id, the third column is hex color code. Refer to "./extra/markershapes.txt" to see the list of shapes implemented in this bokeh-dependent python script. The parameters of the yaml files are described below:
```
 plot_width: plot-width size in pixel
 plot_height: plot-height size in pixel
 pc_x_to_plot: which pc to plot on the x-axis
 pc_y_to_plot: which pc to plot on the y-axis
 fill_alpha: fill-color intensity 
 line_alpha: line-color intensity
 marker_size : size of the markers to be plotted; input should be boolean; True or False
 show_sample_label: whether or not to show the sample label for each dot during hovering. 
```
> **Note:** For plotting large number of samples and populations, increase the plot_width and plot_height size, reduce the marker_size and set show_sample_label to false.

Next, the admixture plot can be also generated with the following command: 
```
python3 plot_interactive_q_mat.py -q <Q_matrix_file> -f <plink_fam_file> -y admixture.yml -c color.txt -o <output_prefix> -s plot_pop_order.txt
```
The python script located in the bin folder of scalepopgen. The q-matrixes are located in the respective output folder of admixture (**./admixture/ADMIXTURE/**) and should be listed in one column of a text file **<Q_matrix_file>**. The file <plink_fam_file> is located in folder **./plink/make_bed/**. Text files **color.txt** and **plot_pop_order.txt** can be created by the user. File **color.txt** has hex color codes listed in one column:
```
#FF0000
#00FF00
#0000FF
#FFFF00
#FF00FF
#00FFFF
#FFA500
#800080
#008000
#800000
```
Similarly, in file **plot_pop_order.txt** population IDs are listed in one column according to the order they should be plotted:
```
Appenzell
Booted
ChamoisColored
Peacock
Saanen
Toggenburg
Bezoar
AlpineIbex
```
The yaml file **admixture.yml** is located in scalepopgen folder **/extra/plots/** and contains parameters described below:
```
 width: plot-width size in pixel
 height: plot-height size in pixel
 bar_width: the width of each sample bar
 sample_label_orientation: degrees of anlge at which the sample labels should be written
 pop_label_orientation: degrees of anlge at which the population labels should be written
 space_pop_group: the width of space between populations
 legend_font_size: font size of the legend
 num_legend_per_col: number of different K per column
 label_font_size: font size of the labels
 fil_alpha: fill-color intensity 
```

For generating the IBS-dist interactive NJ trees, one can use the following command:
```
python3 make_ibs_dist_nj_tree.py -r <outgroup> -i <square_mat_mdist_file> -m <mdist.id_file> -c pop_sc_color.map -y ibs_nj.yml -o <output_prefix>
```
The python script is located in the bin folder of scalepopgen. The ```<outgroup>``` refers to the population to be used for rooting the tree, **<square_mat_mdist_file>** refers to the square matrix of 1-ibs distance between pairwise samples (**./ibs_clustering/plink/calc_1_mins_ibs_dist/**), **<mdist.id_file>** refers to id file generated along with the square matrix (**./ibs_clustering/plink/calc_1_mins_ibs_dist/**, **pop_sc_color.map* is the tab-delimited file containing the first column as pop_id, second column as sample size and the third column as hex color code. Also generated by the workflow and saved in the output folder as **pop_sc_color.map**. The yml file is located in **/extra/plots/** folder. The parameter of the yml files are described below:
```
 width: plot-width size in pixel
 height: plot-height size in pixel
 layout: the tree layout, valid options: 'c','r','d', for circular, right and down layout of the tree
 tip_label_align: whether or not to align the tip labels; input should be boolean; True or False
 tip_label_font_size: the font size of the tip labels, default: "12px"
 edge_widths: the width of edges, default: 1
 node_sizes: the size of nodes, default:6
 node_hover: whether or not to show details info of the node while hovering; input should be boolean; True or False
```
The python script can also be run directly using the file containing tree in newick format. For more details run:```python3 make_ibs_dist_nj_tree.py -h```. As the python script is dependent on Toytree, for more details of these parameters, refer to [Toytree](https://toytree.readthedocs.io/en/latest/) documentation.

For generating the Fst-based NJ tree, one can use the following command:
```
python3 make_fst_dist_nj_tree.py -i <fst_summary_file_generated_by_plink2> -r <outgroup> -o <output_prefix> -y fst_nj.yml -c pop_sc_color.map
``` 
The python script is located in the bin folder of scalepopgen. The yml file is located in **/extra/plots/** folder. The parameter of the yml files are described below:
```
 width: plot-width size in pixel
 height: plot-height size in pixel
 layout: the tree layout, valid options: 'c','r','d', for circular, right and down layout of the tree
 tip_label_align: whether or not to align the tip labels; input should be boolean; True or False
 tip_label_font_size: the font size of the tip labels, default: "12px"
 node_sizes: the size of nodes, default:6
 node_hover: whether or not to show details info of the node while hovering; input should be boolean; True or False
```

The ```<outgroup>``` refers to the population to be used for rooting the tree, **<fst_summary_file_generated_by_plink2>** refers to the table with calculated FST distances between populations located in **./fst_clustering/plink2/**, **pop_sc_color.map* is the tab-delimited file containing the first column as pop_id, second column as sample size and the third column as hex color code. Also generated by the workflow and saved in the output folder as **pop_sc_color.map**.

## References
Please cite the following papers if you use this sub-workflow in your study:

[1] Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome research, 19(9), 1655-1664. https://doi.org/10.1101/gr.094052.109

[2] Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. PLoS genetics, 2(12), e190. https://doi.org/10.1371/journal.pgen.0020190

[3] Price, A. L., Patterson, N. J., Plenge, R. M., Weinblatt, M. E., Shadick, N. A., & Reich, D. (2006). Principal components analysis corrects for stratification in genome-wide association studies. Nature genetics, 38(8), 904-909. https://doi.org/10.1038/ng1847

[4] Francis, R. M. (2017). pophelper: an R package and web app to analyse and visualize population structure. Mol Ecol Resour, 17: 27–32. doi:10.1111/1755-0998.12509

[5] Huerta-Cepas, J. et al.,(2016). ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data, Molecular Biology and Evolution, Volume 33, Issue 6, June 2016, Pages 1635–1638, https://doi.org/10.1093/molbev/msw046.

[6] Eaton, Deren. (2019). Toytree: A minimalist tree visualization and manipulation library for Python. Methods in Ecology and Evolution. 11. 10.1111/2041-210X.13313. 

[7] Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England), 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

[8] Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319 (2017). https://doi.org/10.1038/nbt.3820

## License

MIT

