#!/bin/bash -l

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o identify_and_drop_outlier_low_read_samples.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir=<path to Amplicon_Sequence_Cleaning directory>
project_dir=<path to where you want your outputs saved>
module load R/4.3.1 

# I'd recommend commenting out steps 3b and 3c the first time you run this!

echo "IDENTIFY OUTLIER SAMPLES"
echo "Running on 16S data"
Rscript $script_dir/03_Filter_Samples/03_a_identify_outlier_samples.R \
	-a "16S" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-v <path to file with list of variables to color data by for NMDS plots> \ # see documentation

echo "Running on ITS data"
Rscript $script_dir/03_Filter_Samples/03_a_identify_outlier_samples.R \
	-a "ITS" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-v <path to file with list of variables to color data by for NMDS plots> \ # see documentation

echo "=========================================================="

# Once you're satisfied that step 3a has run properly, uncomment just step 3b, and 
# (optionally) comment out step 3a. Rerun step 3b until you are satisfied that you have \
# identified all outliers and identified a proper drop threshold.

echo "EVALUATE DROP THRESHOLDS"
echo "Running on 16S data"

Rscript $script_dir/03_Filter_Samples/03_b_evaluate_drop_threshold.R \
	-a "16S" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-t <path to file with filter test thresholds table> \ # see documentation
	-o <path to file with outlier sample names list> \ # see documentation
	-v <path to file with variable names list for coloring NMDS plots by> \ # see documentation

echo "Running on ITS data"
Rscript $script_dir/03_Filter_Samples/03_b_evaluate_drop_threshold.R \
	-a "ITS" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-t <path to file with filter test thresholds table> \ # see documentation
	-o <path to file with outlier sample names list> \ # see documentation
	-v <path to file with variable names list for coloring NMDS plots by> \ # see documentation

echo "=========================================================="

# Once step 3b is done, comment it out and then uncomment step 3c. 

echo "DROP OUTLIERS AND LOW READ SAMPLES"
echo "Running on 16S data"
Rscript $script_dir/03_Filter_Samples/03_c_drop_outliers_and_low_read_count_samples.R \
	-a "16S" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-e "Y" \ # Do you want to edit the metadata table to reflect which samples have been dropped and for what reason? [options: Y or N]
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-t <path to file with filter thresholds table> \ # see documentation
	-o <path to file with outlier sample names list> \ # see documentation
	
echo "Running on ITS data"
Rscript $script_dir/03_Filter_Samples/03_c_drop_outliers_and_low_read_count_samples.R \
	-a "ITS" \ # amplicon type [options: 16S or ITS]
	-n "atherton" \ # your last name, used for naming output files so we know who ran the script
	-e "Y" \ # Do you want to edit the metadata table to reflect which samples have been dropped and for what reason? [options: Y or N]
	-p $project_dir \ # path to where outputs should be saved
	-m <path to metadata file> \ # see documentation
	-s $script_dir \
	-t <path to file with filter thresholds table> \ # see documentation
	-o <path to file with outlier sample names list> \ # see documentation

echo "Done."
