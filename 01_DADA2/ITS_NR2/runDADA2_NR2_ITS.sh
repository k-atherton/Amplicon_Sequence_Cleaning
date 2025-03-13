# python /share/pkg.7/bu16s/1.0/install/create_inputs.py --project ST_ITS_NR2_trunc --input_dir ../../Nahuel_Run2/fastq_files/ITS/ --output_dir ./ --fwd _L001_R1_001.fastq.gz --fprimer 
export PROJECTNAME=ST_ITS_NR2_trunc
export INPUTDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/Nahuel_Run2/fastq_files/ITS
export OUTPUTDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/dada2_output/ITS_NR2_trunc
export INTERMEDIATEDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/dada2_output/ITS_NR2_trunc/intermediate
export RUNPARAMETERS=/projectnb/talbot-lab-data/Katies_data/Street_Trees/dada2_output/ITS_NR2_trunc/.runparams
export FWD_FMT=_L001_R1_001.fastq.gz
export REV_FMT=None
export FWD_PRIMER=
export REV_PRIMER=GGACTACHVHHHTWTCTAAT
export PRIMER_END=5
export DADA2_TRUNC_LEN_F=0
export DADA2_TRUNC_LEN_R=0
export CUTADAPT_ARGS=""
export DADA2_ARGS=""
export PAIRED=False
export SCRIPTSDIR=/share/pkg.7/bu16s/1.0/install/scripts
export SILVA_SEQUENCES=/projectnb/microbiome/ref_db/UNITE/sh_refs_qiime_ver9_dynamic_25.07.2023.qza
export SILVA_TAXONOMY=/projectnb/microbiome/ref_db/UNITE/sh_taxonomy_qiime_ver9_dynamic_25.07.2023.qza
# Load modules and inputs
module purge
module load miniconda/4.7.5
module load qiime2/2020.2
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
conda activate $SCC_QIIME2_DIR
