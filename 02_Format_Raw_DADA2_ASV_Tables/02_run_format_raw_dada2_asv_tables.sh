#!/bin/bash -l

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o format_raw_dada2_asv_tables_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir=<path to Amplicon_Sequence_Cleaning directory>
project_dir=<path to where you want your outputs saved>

# Input information: 
# -a: amplicon type [options: 16S or ITS]
# -n: your last name, used for naming output files so we know who ran the script
# -e: do you want to edit the metadata file to add the DADA2 read count? [options: Y or N]
# -p: path to where outputs should be saved 
# -v: path to file with list of ASV table paths, see documentation
# -t: path to file with list of taxonomy table paths, see documentation
# -m: path to metadata file, see documentation
# -c: path to file with negative control naming patterns, see documentation
# -s: path to the scripts directory

echo "Running on 16S data"
module load R/4.3.1 
Rscript $script_dir/02_Format_Raw_DADA2_ASV_Tables/02_format_raw_dada2_asv_tables.R \
	-a "16S" \
	-n "atherton" \
	-e "Y" \
	-p $project_dir \
	-v <path to file with list of ASV table paths> \
	-t <path to file with list of taxonomy table paths> \
	-m <path to metadata file> \
	-c <path to file with negative control naming patterns> \
	-s $script_dir

echo "Running on ITS data"
Rscript $script_dir/02_Format_Raw_DADA2_ASV_Tables/02_format_raw_dada2_asv_tables.R \
	-a "ITS" \
	-n "atherton" \
	-e "Y" \
	-p $project_dir \
	-v <path to file with list of ASV table paths> \
	-t <path to file with list of taxonomy table paths> \
	-m <path to metadata file> \
	-c <path to file with negative control naming patterns> \
	-s $script_dir
 
echo "Done."
