#!/bin/bash -l

# Set SCC project
#$ -P microbiome

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Give job a name
#$ -N memento_diversity

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o Bash_Outputs/MEMENTO_atherton_calculate_diversity_function_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir="/projectnb/talbot-lab-data/Katies_data/Amplicon_Sequence_Cleaning"
project_dir="/projectnb/talbot-lab-data/Katies_data/MEMENTO"

echo "Running on 16S data"
module load R/4.3.1 
Rscript $script_dir/06_Calculate_Diversity_and_Function/06_calculate_diversity_and_function.R \
	-a "16S" \
	-n "atherton" \
	-p $project_dir \
	-m "$project_dir/01_Collect_Data/01_Sample_Metadata/MEMENTO_atherton_16S_sample_metadata_dada2_dropsamples_decontam_20250510.csv" \
	-s $script_dir \
	-e "Y" \
	-r "$project_dir/02_Clean_Data/05_Transform_Data/16S/Rarefied_ASV_Tables/" \
	-t "$project_dir/02_Clean_Data/02_DADA2_ASV_Tables/16S/atherton_16S_taxonomy_allsampletypes_raw_20250507.csv" \
	--pathogen16s "/projectnb/talbot-lab-data/Katies_data/M-BUDS/03_Data_Analysis/02_Identify_Pathogens/atherton_16S_mbpd_blast_match_table_20240925.txt"

echo "Running on ITS data"
Rscript $script_dir/06_Calculate_Diversity_and_Function/06_calculate_diversity_and_function.R \
	-a "ITS" \
	-n "atherton" \
	-p $project_dir \
	-m "$project_dir/01_Collect_Data/01_Sample_Metadata/MEMENTO_atherton_ITS_sample_metadata_dada2_dropsamples_decontam_20250510.csv" \
	-s $script_dir \
	-e "Y" \
	-r "$project_dir/02_Clean_Data/05_Transform_Data/ITS/Rarefied_ASV_Tables/" \
	-t "$project_dir/02_Clean_Data/02_DADA2_ASV_Tables/ITS/atherton_ITS_taxonomy_allsampletypes_raw_20250507.csv"

echo "Done."
