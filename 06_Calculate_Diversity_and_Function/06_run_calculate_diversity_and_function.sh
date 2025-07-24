#!/bin/bash -l

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o calculate_diversity_function_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir="<path to Amplicon_Sequence_Cleaning directory>"
project_dir="<path to project directory>"

echo "Running on 16S data"
module load R/4.3.1 

# Input information: 
# -a: amplicon sequence type [options: 16S or ITS]
# -n: your last name, used for naming output files so we know who ran the script
# -p: path to where outputs should be saved
# -m: path to metadata file, see documentation
# -s: path to script directory
# -e: do you want to edit the metadata file to add the DADA2 read count? [options: Y or N]
# -r: path to directory with rarefied data, see documentation
# -t: path to taxonomy table created by Step 2, see documentation
# --pathogen16s: path to 01_blast_16s_pathogens.sh output, only needed for 16S data; see documentation

Rscript $script_dir/06_Calculate_Diversity_and_Function/06_calculate_diversity_and_function.R \
	-a "16S" \
	-n "atherton" \
	-p $project_dir \
	-m <path to metadata file> \
	-s $script_dir \
	-e "Y" \
	-r <path to directory with rarefied data> \
	-t <path to taxonomy table> \
	--pathogen16s <path to 01_blast_16s_pathogens.sh output>

echo "Running on ITS data"
Rscript $script_dir/06_Calculate_Diversity_and_Function/06_calculate_diversity_and_function.R \
	-a "ITS" \
	-n "atherton" \
	-p $project_dir \
	-m <path to metadata file> \
	-s $script_dir \
	-e "Y" \
	-r <path to directory with rarefied data> \
	-t <path to taxonomy table> \

echo "Done."
