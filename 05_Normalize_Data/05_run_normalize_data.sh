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
Rscript $script_dir/05_Normalize_Data/05_normalize_data.R \
	-a "16S" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-b "Y" \ # do you want to run batch correction?
	-r "Y" \ # do you want to rarefy your data?
	-t <number> \ # threshold for rarefying data, remove if you are not rarefying, or if you want to use the default of the minimum of the sample read sums
	-l "Y" \ # do you want to clr-transform your data?
	-z "Y" \ # do you want to z-score transform your data?
	-d "Y" \ # do you want to calculate the aitchison distance matrix for your data?
	-v <path to text file with variables for batch correction plots> \ # see documentation, only needed if you use -b "Y"
	--batchvariable <variable> \ # variable that data will be batch corrected by, only needed if you use -b "Y"
	--covariates <path to text file with variables to maintain the signature of in batch correction # see documentation, only needed if you use -b "Y"

echo "Running on ITS data"
Rscript $script_dir/05_Normalize_Data/05_normalize_data.R \
	-a "ITS" \ \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-b "Y" \ # do you want to run batch correction?
	-r "Y" \ # do you want to rarefy your data?
	-t <number> \ # threshold for rarefying data, remove if you are not rarefying, or if you want to use the default of the minimum of the sample read sums
	-l "Y" \ # do you want to clr-transform your data?
	-z "Y" \ # do you want to z-score transform your data?
	-d "Y" \ # do you want to calculate the aitchison distance matrix for your data?
	-v <path to text file with variables for batch correction plots> \ # see documentation, only needed if you use -b "Y"
	--batchvariable <variable> \ # variable that data will be batch corrected by, only needed if you use -b "Y"
	--covariates <path to text file with variables to maintain the signature of in batch correction # see documentation, only needed if you use -b "Y"


echo "Done."
