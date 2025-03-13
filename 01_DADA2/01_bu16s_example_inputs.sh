export PROJECTNAME=<name_of_DADA2_run>
export INPUTDIR=<path_to_input_directory>
export OUTPUTDIR=$(pwd)
export INTERMEDIATEDIR=$(pwd)/intermediate
export RUNPARAMETERS=$(pwd)/.runparams
export FWD_FMT=_L001_R1_001.fastq.gz # this is the format that Tufts uses, change if your sequencing facility has a different format
export REV_FMT=_L001_R2_001.fastq.gz # this is the format that Tufts uses, change if your sequencing facility has a different format
export FWD_PRIMER=#if your primers were removed by your sequencing facility (e.g. Tufts), leave these blank
export REV_PRIMER=#if your primers were removed by your sequencing facility (e.g. Tufts), leave these blank
export PRIMER_END=5
export DADA2_TRUNC_LEN_F=0
export DADA2_TRUNC_LEN_R=0
export CUTADAPT_ARGS=""
export DADA2_ARGS=""
export PAIRED=True
export SCRIPTSDIR=/share/pkg.7/bu16s/1.0/install/scripts
export SILVA_SEQUENCES=/projectnb/microbiome/ref_db/silva_132_99_16S.qza # change if you need a different version of SILVA or if you are annotating ITS sequences
export SILVA_TAXONOMY=/projectnb/microbiome/ref_db/silva_132_99_majority_taxonomy.qza # change if you need a different version of SILVA or if you are annotating ITS sequences
# Load modules and inputs
module purge
module load miniconda/4.7.5
module load qiime2/2020.2
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
conda activate $SCC_QIIME2_DIR