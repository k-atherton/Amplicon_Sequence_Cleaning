#!/bin/bash -l

# Set SCC project
#$ -P microbiome

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Give job a name
#$ -N m-buds_blast_16s_pathogens

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o Bash_Outputs/M-BUDS_ED_atherton_blast_bacteria_pathogens_output.txt

# Request more cores and memory
#$ -pe omp 8
#$ -l mem_per_core=16G


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

script_dir="/projectnb/talbot-lab-data/Katies_data/Amplicon_Sequence_Cleaning"
project_dir="/projectnb/talbot-lab-data/Katies_data/M-BUDS_ED"

# To make the query file, you will need to concatenate your DADA2 representative sequences
# files from your 16S data into one file.
# Ensure that your DADA2 individual run folders for your 16S data all start with 16S in 
# order to run this line for only your bacterial data. Otherwise, adjust the path to the
# representative sequences files. 
echo "Concatenating multiple runs' 16S representative sequences."
for i in $project_dir/02_Clean_Data/01_DADA2/16S_*/intermediate/dada2/*rep_seqs.fasta/dna-sequences.fasta;
	do cat "$i" >> $project_dir/02_Clean_Data/06_Calculate_Diversity_and_Function/16s_representative_sequences_with_duplicates.fasta;
done

# Then, remove the duplicate representative sequences from across runs
echo "Removing duplicated representative sequences."
module load seqkit
seqkit rmdup -s < $project_dir/02_Clean_Data/06_Calculate_Diversity_and_Function/16s_representative_sequences_with_duplicates.fasta > $project_dir/02_Clean_Data/06_Calculate_Diversity_and_Function/16s_representative_sequences_unique.fasta

# Unzip the MBPD pathogen database. 
echo "Unzipping the MBPD database."
unzip $project_dir/00_Databases/pathogen.fasta -d $project_dir/00_Databases/

# BLAST against the Multiple Bacterial Pathogen Database
echo "BLASTing against the pathogen database."
module load blast+
blastn -query $project_dir/02_Clean_Data/06_Calculate_Diversity_and_Function/16s_representative_sequences_unique.fasta\
	-db $project_dir/00_Databases/pathogen.fasta \
	-perc_identity 95 \
	-qcov_hsp_perc 100 \
	-outfmt 7 \
	-out $project_dir/02_Clean_Data/06_Calculate_Diversity_and_Function/16s_pathogens_blast_output.txt

# Format the output into a table.
echo "Formatting BLAST output."
printf "asv_id\tpathogen_id\t\perc_id\talignment_length\tmismatches\tgap_opens\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score\n" >> 16s_pathogens_results.txt
grep -v "#" >> 16s_pathogens_results.txt

echo "Done."
