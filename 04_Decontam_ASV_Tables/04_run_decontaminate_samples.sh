#!/bin/bash -l

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o decontaminate_sequences_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir=<path to Amplicon_Sequence_Cleaning directory>
project_dir=<path to where you want your outputs saved>

echo "Running on 16S data"
module load R/4.3.1 
Rscript $script_dir/04_Decontaminate_Samples/04_decontaminate_samples.R \
	-a "16S" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-e "Y" \ # do you want to edit the metadata file to add the DADA2 read count? [options: Y or N]
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir

echo "Running on ITS data"
Rscript $script_dir/04_Decontaminate_Samples/04_decontaminate_samples.R \
	-a "ITS" \ \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-e "Y" \ # do you want to edit the metadata file to add the DADA2 read count? [options: Y or N]
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir

echo "Done."
