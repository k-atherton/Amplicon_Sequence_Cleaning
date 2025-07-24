#!/bin/bash -l

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o normalize_data_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir=<path to Amplicon_Sequence_Cleaning directory>
project_dir=<path to where you want your outputs saved>

echo "Running on 16S data"
module load R/4.3.1 

# Input Information:
# -a: amplicon type [options: 16S or ITS]
# -n: your last name, used for naming output files so we know who ran the script
# -p: path to where outputs should be saved
# -m: path to metadata file, see documentation
# -s: path to script directory
# -b: do you want to check for and, if needed, run batch correction? [options: Y or N]
# -r: do you want to rarefy your data?  [options: Y or N]
# -t: threshold for rarefying data, remove if you are not rarefying, or if you want to use the default of the minimum of the sample read sums
# -l: do you want to clr-transform your data? [options: Y or N]
# -z: do you want to z-score transform your data? [options: Y or N]
# -d: do you want to calculate the aitchison distance matrix for your data? [options: Y or N]
# -v: path to text file with variables for batch correction plots, see documentation; only needed if you use -b "Y"
# --batchvariable: variable that data will be batch corrected by, only needed if you use -b "Y"
# --covariates: path to text file with variables to maintain the signature of in batch correction, see documentation; only needed if you use -b "Y"

Rscript $script_dir/05_Normalize_Data/05_normalize_data.R \
	-a "16S" \
	-n "atherton" \
	-p $project_dir \
	-m <path to metadata file> \
	-s $script_dir \
	-b "Y" \
	-r "Y" \
	-t <number> \
	-l "Y" \
	-z "Y" \
	-d "Y" \
	-v <path to text file with variables for batch correction plots> \
	--batchvariable <variable> \
	--covariates <path to text file with variables to maintain the signature of in batch correction> # see documentation, only needed if you use -b "Y"

echo "Running on ITS data"
Rscript $script_dir/05_Normalize_Data/05_normalize_data.R \
	-a "ITS" \
	-n "atherton" \
	-p $project_dir \
	-m <path to metadata file> \
	-s $script_dir \
	-b "Y" \
	-r "Y" \
	-t <number> \
	-l "Y" \
	-z "Y" \
	-d "Y" \
	-v <path to text file with variables for batch correction plots> \
	--batchvariable <variable> \
	--covariates <path to text file with variables to maintain the signature of in batch correction> # see documentation, only needed if you use -b "Y"

echo "Done."
