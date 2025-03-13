# python /share/pkg.7/bu16s/1.0/install/create_inputs.py --project ST_16S_NR1_trunc --input_dir ../../Nahuel_Run1/fastq_files/16S/ --output_dir ./ --fwd _L001_R1_001.fastq.gz --rev _L001_R2_001.fastq.gz --fprimer  --rprimer  --dada2_args=--p-trunc-q 30
export PROJECTNAME=ST_16S_ReRun
export INPUTDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/ReRun/fastq_files/16S
export OUTPUTDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/dada2_output/16S_ReRun
export INTERMEDIATEDIR=/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/dada2_output/16S_ReRun/intermediate
export RUNPARAMETERS=/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/dada2_output/16S_ReRun/.runparams
export FWD_FMT=_L001_R1_001.fastq.gz
export REV_FMT=_L001_R2_001.fastq.gz
export FWD_PRIMER=
export REV_PRIMER=
export PRIMER_END=5
export DADA2_TRUNC_LEN_F=0
export DADA2_TRUNC_LEN_R=0
export CUTADAPT_ARGS=""
export DADA2_ARGS="--p-trunc-q 0"
export PAIRED=True
export SCRIPTSDIR=/share/pkg.7/bu16s/1.0/install/scripts
export SILVA_SEQUENCES=/projectnb/talbot-lab-data/Katies_data/Databases/dada2_silva/silva_138.1_refseqs.qza
export SILVA_TAXONOMY=/projectnb/talbot-lab-data/Katies_data/Databases/dada2_silva/silva_138.1_tax.qza
# Load modules and inputs
module purge
module load miniconda/4.7.5
module load qiime2/2020.2
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
conda activate $SCC_QIIME2_DIR
