# Amplicon Sequence Cleaning Workflow
by Kathryn Atherton

The purpose of this project is to document the amplicon sequence cleaning process and decisions that I use in my research so that others can reproduce this process on my data or on their own datasets. 

My workflow starts after sequences have been clustered and annotated by DADA2 (I use the [BU16S pipeline](https://github.com/Boston-University-Microbiome-Initiative/BU16s) written by Dr. Michael Silverstein). It guides a user through testing for batch effect and sequence depth influence, identifying and dropping outliers, removing contaminants from samples based on negative control sequence data, and normalizing data appropriately for various downstream analyses, including rarefaction, CLR and z-score transformations, and calculating Aitchison distance. The workflow also allows for annotating fungal and bacterial functional guilds.

## Requirements
To run this pipeline, you will need the following R packages (in parentheses are the versions I used; I ran this pipeline in R version 4.3.1):
- compositions (v.)
- dplyr (v.1.1.4)
- geosphere (v.)
- ggplot2 (v.3.5.1)
- gridExtra (v.2.3)
- optparse (v.1.7.5)
- phyloseq (v.1.44.0)
- purrr (v.1.0.4)
- readr (v.2.1.4)
- readxl (v.1.4.2)
- robCompositions (v.2.4.1)
- sva (v.)
- tidyr (v.)
- tidyverse (v.)
- vegan (v.2.6-4)
- vroom (v1.6.3)

## Workflow
### Step 1: Run DADA2
**Filename: `01_DADA2/01_bu16s_example_inputs.sh`**
- Summary: This script will run BU16S, a pipeline that runs DADA2 to perform quality control on your samples and then assign taxonomy. 
- Required inputs:
  - project name + run ID (example: M-BUDS_Run1)
  - path to directory where the sequencing run fastq files are
  - output directory where results should be saved [default is the present working directory]
  - file naming format of the forward/reverse read files (example: Tufts names the forward file `_L001_R1_001.fastq.gz` and the reverse file `_L001_R2_001.fastq.gz`)
  - forward and reverse primer sequences (leave blank if primers were removed by sequencing facility)
  - path to taxonomic annotation database (defaults to BU Microbiome Initiative SILVA version 132)
  - advanced DADA2 parameters, if necessary (i.e. primer end, truncation length, cutadapt arguments, paired sequencing, etc).
- Outputs:
  - ASV table: rows are ASVs, columns are samples, elements are counts of ASVs in each sample
  - Taxonomy table: matches ASV ID to taxonomic information
  - Other intermediate files, including representative sequences file (see BU16S documentation for more information)

See `01_DADA2/00_how_to_run_dada2.txt` for more information on how to use the example inputs script to run BU16S. If you are working on the BU SCC, you will need to load the BU16S module (`module load bu16s`) and then run the script with the following command: `qsub -P <SCC project name> -N <job name> $SCC_BU16S_DIR/bu16s.qsub <inputs file name>.sh`


## Step 2: Format Raw DADA2 ASV tables
**Filename: `02_Format_Raw_DADA2_ASV_Tables/02_run_format_raw_dada2_asv_tables.sh`**
- Summary: This script will help you to format the DADA2 outputs for downstream analysis. 
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Should the metadata file be edited? (flag: -e or --edit): do you want to add the DADA2 sequence count to the metadata file? options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - File with paths to ASV tables output by DADA2 (flag: -v or --asvtable): the file should be a .txt where each line is the full path to each DADA2 output ASV table. See `example_asv_table_paths_file.txt`
    - **Make this file yourself, but the file should point to the ASV tables created by BU16S in Step 1!**
  - File with paths to taxonomy tables output by DADA2 (flag: -t or --taxonomy): the file should be a .txt where each line is the full path to each DADA2 output taxonomy table.
    - **Make this file yourself, but the file should point to the taxonomy tables created by BU16S in Step 1! See `example_taxonomy_paths_file.txt`**
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.). The metadata file should also have a column called is_control with values TRUE (in rows for negative controls) or FALSE (in rows for samples).
  - File with negative control naming patterns, matched to the sample types that the negative controls were amplified and sequenced with (flag: -c or --negcontrols): the file should be a .csv and have a column called sample_type with the sample type that matches the sample_type names in the metadata file and another column called negative_control_name_pattern that contains strings that define the naming patterns of the negative controls. You can use the whole negative control sample names, or just a pattern that will allow a function like grep to find the negative control samples.
    - **Make this file yourself! See `example_negative_control_name_file.csv`**
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
- Outputs:
  - Raw (i.e. not filtered) formatted ASV tables: one for all sample types and one for each individual sample type. Saved as .csv
  - Raw (i.e. not filtered) formatted taxonomy tables: for all sample types and one for each individual sample type. Saved as .csv
  - Raw (i.e. not filtered) phyloseq objects: one for each individual sample type. The phyloseq objects contain the ASV table, taxonomy, and metadata for each sample. Saved as .RDS                          

To run this script, adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 02_run_format_raw_dada2_asv_tables.sh>`


## Step 3: Filter Samples
### Step 3a: Identify outlier samples
**Filename: `03_Filter_Samples/03_a_identify_outlier_samples.sh`**
- Summary: This script will help you to identify any outlier samples by making a scatterplot with the sample names as the points. 
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.). The metadata file should also have a column called is_control with values TRUE (in rows for negative controls) or FALSE (in rows for samples).
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
  - File with metadata column names for NMDS plotting (flag: -v or --variables): the file should be a .txt and each row should have a column name within the metadata file. The variables included in this file will be used to color NMDS plots of the samples to visualize the data structure as determined by these variables. Ideally, the variables will be factors, but continuous variables are also okay.
    - **Make this file yourself! See example file: `example_nmds_variables_file.txt`**
- Outputs:
  - Figure: NMDS of the samples where the points are the sample names to ID outliers. Use your best judgement to remove any outliers -- I typically remove any samples whose names I can clearly read that are far from any other samples. You'll know it when you see it.
  - Figure: multipanel NMDS of the data structure before removing outliers, coloring the samples by different variables.
  - Figure: histogram of the sample sequencing depths before dropping any samples.

To run this script, comment out lines 45-98 (put a "#" in front of each line) and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 03_run_filter_samples.sh>`

### Step 3b: Evaluate drop thresholds
**Filename: `03_Filter_Samples/03_b_evaluate_drop_thresholds.R`**
- Summary: This script helps you to evaluate multiple low sequence count thresholds before you officially decide what to use moving forward. It will output the results from multiple drop threshold scenarios so that you can compare them side by side. 
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.). The metadata file should also have a column called is_control with values TRUE (in rows for negative controls) or FALSE (in rows for samples).
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
  - File with sequence read depth thresholds to test (flag: -t or --thresholds): the file should be a .csv with three columns: sample_type (which contains the sample types matching the sample_type column in the metadata file, not including negative controls), threshold1, and threshold2 (numbers used to test the threshold at which to drop samples if they have a post-DADA2 read count less than the threshold -- I typically evaluate thresholds between 5000 and 10000).
    - **Make this file yourself. See example file: `example_filter_test_thresholds_file.csv`** 
  - File with the names of outlier samples to drop (flag: -o or --outliers): the file should be a .txt and each row should have a sample name which you identified in Step 3a.
    - **Make this file yourself based on your findings from step 3a. See example file: `example_outliers_file.txt`**
  - File with metadata column names for NMDS plotting (flag: -v or --variables): the file should be a .txt and each row should have a column name within the metadata file. The variables included in this file will be used to color NMDS plots of the samples to visualize the data structure as determined by these variables. Ideally, the variables will be factors, but continuous variables are also okay.
    - **Make this file yourself. See example file: `example_nmds_variables_file.txt`**
- Outputs:
  - Figure: NMDS of the samples where the points are the sample names to ID outliers, after the previously identified outliers were removed. Use this to see if you missed any outliers.
  - Figure: multipanel NMDS of the data structure after removing outliers, coloring the samples by different variables.
  - Figure: multipanel NMDS of the data structure after removing outliers and samples with a sequence count < threshold1, coloring the samples by different variables.
  - Figure: multipanel NMDS of the data structure after removing outliers and samples with a sequence count < threshold2, coloring the samples by different variables.
    - The main thing you're looking for here is that there isn't a significant restructuring of your data before/after dropping samples. Based on these NMDS figures, decide which threshold to go with. The bash output will list the number of samples and names of samples dropped for each threshold. Use that to determine if too many samples total or if too many samples from one treatment were dropped by either threshold in order to make your decision.
   
To run this script, comment out lines 24-43 and 77-98 (put a "#" in front of each line) and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 03_run_filter_samples.sh>`

### Step 3c: Drop outliers and low read count samples
**Filename: `03_Filter_Samples/03_c_drop_outliers_and_low_read_count_samples.R`**
- Summary: This script will finally remove the outlier samples and filter out any samples with a low read count.
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Edit metadata file (flag: -e or --edit): Do you want to make edits to the metadata file; options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.). The metadata file should also have a column called is_control with values TRUE (in rows for negative controls) or FALSE (in rows for samples).
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
  - File with the threshold to drop samples with low reads (flag: -t or --thresholds): the file should be a .csv with three columns: sample_type (which contains the sample types matching the sample_type column in the metadata file, not including negative controls), and threshold (number used to drop samples if they have a post-DADA2 read count less than the threshold).
    - **Make this file yourself, based on your findings from Step 3b. See example file: `example_filter_thresholds_file.csv`** 
  - File with the names of outlier samples to drop (flag: -o or --outliers): the file should be a .txt and each row should have a sample name which you identified in Step 3a.
    - **Make this file yourself based on your findings from Step 3a. See example file: `example_outliers_file.txt`**
  - File with metadata column names for NMDS plotting (flag: -v or --variables): the file should be a .txt and each row should have a column name within the metadata file. The variables included in this file will be used to color NMDS plots of the samples to visualize the data structure as determined by these variables. Ideally, the variables will be factors, but continuous variables are also okay.
    - **Make this file yourself. See example file: `example_nmds_variables_file.txt`**
- Outputs:
  - Figure: histogram of the sample sequencing depths after dropping all outliers and low read count samples.
  - File: ASV table without dropped samples (csv and phyloseq RDS versions)
  - File: ASV table without dropped samples, including negative controls (phyloseq RDS)
  - If -e flag == Y: File: metadata table with new column "sequences_dropped". Containes values "No" for samples kept for downstream analysis, "Outlier" for samples identified as outliers in steps 3a, or "Low read count" for samples dropped because they have a read count < threshold.

To run this script, comment out lines 24-73 (put a "#" in front of each line) and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 03_run_filter_samples.sh>`

## Step 4: Decontaminate ASV Tables
**Filename: `04_Decontam_ASV_Tables/04_decontaminate_samples.sh`**
- Summary: This script will remove any contaminating ASVs. If you have negative control sample(s) for a sample type and sequencing batch, it will use the prevalence of ASVs in the negative controls vs. the true samples to identify contaminants. Otherwise, it will use the frequency of ASVs vs. the DNA concentration in the sample before you pooled your samples for sequencing to identify contaminants. 
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Edit metadata file (flag: -e or --edit): Do you want to make edits to the metadata file; options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.) and a column name called sequencing_batch, which defines the sequencing batches. The metadata file should also have a column called is_control with values TRUE (in rows for negative controls) or FALSE (in rows for samples).
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
- Outputs:
  - Figure: multipanel figure of scatterplots for each sequencing run plotting the pre-decontam sample read count vs. the rank of the read count for the sequencing run. Samples are colored by their control status.
  - Figure (prevalence method): If you have negative controls for sequencing runs, multipanel figure of scatterplots for each sequencing run plotting the prevalence of ASVs in negative controls vs. in true samples, colored by their contaminant status.
  - Figure (frequency method): If you don't have negative controls for sequencing runs, multipanel figure of scatterplots, showing the frequency of ASVs vs the DNA concentrations of the samples they appear in. This figure only shows plots for ASVs determined to be contaminants.
  - Figure: histogram of sequencing depth after decontam. 
  - File: ASV table with contaminating ASVs removed and negative control samples dropped (csv and phyloseq RDS versions)
  - File: Table of the taxonomy of the ASVs removed as contaminants (csv, one per sequencing run)
  - If -e flag == Y: File: metadata table with new column "decontam_seq_count". Contains the summed read counts for samples after decontam.
 
To run this script, and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 04_run_decontaminate_samples.sh>`

## Step 5: Normalize Data
**Filename: `05_Normalize_Data/05_normalize_data.sh`**
- Summary: This script will check if batch correction across your determined variable is necessary, and then perform it if necessary. Then, it will rarefy, CLR transform, Z-score transform, and/or calculate the Aitchison Distance Matrix for your data.
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Edit metadata file (flag: -e or --edit): Do you want to make edits to the metadata file; options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.).
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
  - Batch correction (flag: -b or --batchcorrect): Do you want to check for a batch effect and correct for a batch effect, if the check recommends it?; options: Y or N [default = N]
  - Rarefaction (flag: -r or --rarefy): Do you want to rarefy your data?; options: Y or N [default = N]
  - Rarefaction threshold (flag -t or --threshold): the minimum threshold to rarefy to. Remove if you are not rarefying, or if you want to use the default of the minimum of the sample read sums.
  - CLR transformation (flag -l or --clr): Do you want to CLR transform your data?; options: Y or N [default = N]
  - Z-Score transformation (flag -z or --zscore): Do you want to Z-Score transform your data?; options: Y or N [default = N]
  - Aitchison Distance calculation (flag -d or --aitchisondistance): Do you want to calculate the Aitchison Distance Matrix for your data?; options: Y or N [default = N]
  - Batch Correction Variables (flag -v or --variables): path to text file with variables you want to color by in your batch correction plots. One variable per line. Only include if you use -b "Y".
  - Variable to Batch Correct by (flag --batchvariable): Column name of the variable that you want to check for a batch effect with. Only include if you use -b "Y".
    - - **Make this file yourself! See example file: `example_batch_correction_variables_file.txt`**
  - Batch Correction Covariates (flag --covariates): Path to text file with variables to maintain the signature of in batch correction. Only include if you use -b "Y".
    - **Make this file yourself! See example file: `example_batch_correction_covariates_file.txt`**
- Outputs:
  - If -b == "Y":
    - File: A batch-corrected ASV table (csv).
    - If batch correction was done, any downstream outputs use the batch-corrected data.
  - If -r == "Y":
    - File: A rarefied ASV table for each sample type (csv and phyloseq object)
    - Figure: Multipanel figure of scatter plots of the post-rarefied vs. pre-rarefied count where points are ASVs, colored by Phylum and by Sample. One per sample type.
    - Figure: Rarefaction curves for each sample type.
  - If -c == "Y":
    - File: A CLR-transformed ASV table, one per sample type plus one of all combined sample types (csv).
  - If -z == "Y":
    - File: A z-score transformed ASV table, one per sample type plus one of all combined sample types (csv).
  - If -d == "Y":
    - File: An Aitchison distance matrix, one per sample type plus one of all combined sample types (csv).
   
To run this script, and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 05_run_normalize_data.sh>`


## Step 6: Calculate Diversity and Function
**Filename: `06_Calculate_Diversity_and_Function/01_blast_16s_pathogens.sh`**
- Summary: This script takes the representative sequences for each ASV created by DADA2 and BLASTs them against the Multiple Bacterial Pathogen Database. Default parameters keep hits that are at least 95% similar to 100% coverage of the query sequence.
- Required inputs:
  - Script directory (variable: script_dir): path to the function script directory where the 00_functions.R script is held
  - Working directory (variable: project_dir): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
- Outputs:
  - File: Table of the top pathogen hits from the Multiple Bacterial Pathogen Database.
 
To run this script, and adjust the script and project directories and run `qsub -P <SCC project name> -N <job name> 01_blast_16s_pathogens.sh>`

**Filename: `06_Calculate_Diversity_and_Function/06_calculate_diversity_and_function.R`**
- Summary: This script will calculate diversity metrics and the percent relative abundances of functional guilds in each sample. For bacterial data, functional guilds are compared against an in-house database (00_Databases/werbin_bacterial_functional_groups_20200220.csv). For fungal data, functional guilds are compared against the FungalTraits database (00_Databases/fungal_traits_database.xlsx). 
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Edit metadata file (flag: -e or --edit): Do you want to make edits to the metadata file; options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - Metadata file (flag: -m or --metadata): the file should be a .csv.
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
  - Rarefied data folder path (flag: -r or --raredata): Path to the folder where all rarefied data phyloseq objects for all sample types are stored. 
  - Taxonomy file (flag: -t or --taxonomy): Path to the taxonomy file created in step 2.
  - Pathogen BLAST results (flag --pathogen16s): Path to the MBPD BLAST results. Remove if you are not working with 16S data.
- Outputs:
  - File: The percent relative abundance of the functional guilds in each sample. Columns are functional guilds and rows are samples.
  - If amplicon == "16S": File: The percent relative abundance of plant, animal, and zoonotic pathogens. Columns are pathogen types and rows are samples.
  - If -e == "Y": A new metadata table with Shannon, Simpson, Inverse Simpson, and Fisher diversity metrics and the percent relative abundance of functional guilds (and, if amplicon == "16S", percent relative abundance of bacterial pathogen types) concatenated to the end of the file.

To run this script, and adjust the inputs for each respective flag within the script and then run `qsub -P <SCC project name> -N <job name> 06_run_calculate_diversity_and_function.sh>`

