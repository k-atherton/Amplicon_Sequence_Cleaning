# Amplicon Sequence Cleaning Workflow
by Kathryn Atherton

The purpose of this project is to document the amplicon sequence cleaning process and decisions that I use in my research so that others can reproduce this process on my data or on their own datasets. 

My workflow starts after sequences have been clustered and annotated by DADA2 (I use the [BU16S pipeline](https://github.com/Boston-University-Microbiome-Initiative/BU16s) written by Dr. Michael Silverstein). It guides a user through testing for batch effect and sequence depth influence, identifying and dropping outliers, removing contaminants from samples based on negative control sequence data, and normalizing data appropriately for various downstream analyses, including rarefaction, CLR and z-score transformations, and calculating Aitchison distance. The workflow also allows for annotating fungal and bacterial functional guilds.

## Requirements
To run this pipeline, you will need the following R packages:
- dplyr (v.1.1.4)
- optparse (v.1.7.5)
- phyloseq (v.1.44.0)
- purrr (v.1.0.4)
- readr (v.2.1.4)
- readxl (v.1.4.2)
- vroom (v1.6.3)

## Workflow
### STEP 1: Run DADA2
**Filename: `01_DADA2/01_bu16s_example_inputs.sh`**
- Required inputs:
  - project name + run ID (example: M-BUDS_Run1)
  - path to directory where the sequencing run fastq files are
  - output directory where results should be saved [default is the present working directory]
  - file naming format of the forward/reverse read files (example: Tufts names the forward file `\_L001\_R1\_001.fastq.gz` and the reverse file `\_L001\_R2\_001.fastq.gz`)
  - forward and reverse primer sequences (leave blank if primers were removed by sequencing facility)
  - path to taxonomic annotation database (defaults to BU Microbiome Initiative SILVA version 132)
  - advanced DADA2 parameters, if necessary (i.e. primer end, truncation length, cutadapt arguments, paired sequencing, etc).
- Outputs:
  - ASV table: rows are ASVs, columns are samples, elements are counts of ASVs in each sample
  - Taxonomy table: matches ASV ID to taxonomic information
  - Other intermediate files, including representative sequences file (see BU16S documentation for more information)

See `01_DADA2/00_how_to_run_dada2.txt` for more information on how to use the example inputs script to run BU16S. If you are working on the BU SCC, you will need to load the BU16S module (module load bu16s) and then run the script with the following command: `qsub -P <SCC project name> -N <job name> $SCC_BU16S_DIR/bu16s.qsub <inputs file name>.sh`

## Step 2: Format Raw DADA2 ASV tables
**Filename: `02_Format_Raw_DADA2_ASV_Tables/02_run_format_raw_dada2_asv_tables.sh`**
- Required inputs:
  - Amplicon type (flag: -a or --amplicon): amplicon type of the dataset; options: 16S or ITS [default = 16S]
  - Last name (flag: -n or --name): your last name for output file naming scheme [default = atherton]
  - Should the metadata file be edited? (flag: -e or --edit): do you want to add the DADA2 sequence count to the metadata file? options: Y or N [default = N]
  - Working directory (flag: -p or --pwd): Directory for saving the outputs of the script [default is the result of the function getwd() in R]
  - File with paths to ASV tables output by DADA2 (flag: -v or --asvtable): the file should be a .txt where each line is the full path to each DADA2 output ASV table. See `example_asv_table_paths_file.txt`
  - File with paths to taxonomy tables output by DADA2 (flag: -t or --taxonomy): the file should be a .txt where each line is the full path to each DADA2 output taxonomy table. See `example_taxonomy_paths_file.txt`
  - Metadata file (flag: -m or --metadata): the file should be a .csv and have a column name called sample_type, which defines the sample types (e.g. "leaf", "negative control", "soil", etc.).
  - File with negative control naming patterns, matched to the sample types that the negative controls were amplified and sequenced with (flag: -c or --negcontrols): the file should be a .csv and have a column called sample_type with the sample type that matches the sample_type names in the metadata file and another column called negative_control_name_pattern that contains strings that define the naming patterns of the negative controls. You can use the whole negative control sample names, or just a pattern that will allow a function like grep to find the negative control samples. See `example_negative_control_name_file.csv`
  - Script directory (flag: -s or --scriptdir): path to the function script directory where the 00_functions.R script is held
- Outputs:
  - Raw (i.e. not filtered) formatted ASV tables: one for all sample types and one for each individual sample type. Saved as .csv
  - Raw (i.e. not filtered) formatted taxonomy tables: for all sample types and one for each individual sample type. Saved as .csv
  - Raw (i.e. not filtered) phyloseq objects: one for each individual sample type. The phyloseq objects contain the ASV table, taxonomy, and metadata for each sample. Saved as .RDS                          

To run this script, adjust the inputs for each respective flag within the script and then run qsub -P <SCC project name> -N <job name> 02_run_format_raw_dada2_asv_tables.sh>
