### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(phyloseq)
library(ggplot2)
library(robCompositions)
library(gridExtra)
library(vegan)
library(readr)

### SCRIPT SETUP ##############################################################
print("SETTING UP SCRIPT:")
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/M-BUDS/"

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default="16S", 
              help="amplicon dataset to filter; options: 16S or ITS [default= %default]", 
              metavar="amplicon"),
  make_option(c("-n", "--name"), type="character", default="atherton", 
              help="last name for output file naming scheme [default= %default]", 
              metavar="your_name"),
  make_option(c("-p", "--pwd"), type="character", default=getwd(),
              help="the directory for saving the outputs of this script [default= %default]",
              metavar="working_directory"),
  make_option(c("-m", "--metadata"), type="character", default=NA,
              help="path to sample metadata file", metavar="metadata"),
  make_option(c("-s","--scriptdir"), type="character", 
              default="/projectnb/talbot-lab-data/Katies_data/Amplicon_Sequence_Cleaning/",
              help="path to function script directory, [default = %default]",
              metavar="script_directory"),
  make_option(c("-t", "--threshold"), type="numeric", default=NA, 
              help="a .csv file that has three columns: sample_type (which matches the sample types (not including negative controls) in the metadata file) and threshold1 and threshold2 (defines the minimum number of reads needed for a sample to be kept for downstream analysis; recommend ~5000-8000 and ~8000-10000)[default= %default]",
              metavar="threshold"),
  make_option(c("-o", "--outliers"), type="character", default=NA,
              help="text file with list of outlier samples to remove from any sample type",
              metavar="outliers")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

amplicon <- opt$amplicon
yourname <- opt$name
edit_metadata <- opt$edit
if(endsWith(opt$pwd, "/")){
  pwd <- opt$pwd
} else{
  pwd <- paste0(opt$pwd, "/")
}
script_dir <- opt$scriptdir
print(opt$metadata)
if (!is.na(opt$metadata)) {
  metadata_path <- opt$metadata
} else {
  stop("Path to metadata file must be provided. See script usage (--help)")
}
if (!is.na(opt$threshold)) {
  threshold_path <- opt$threshold
  thresholds <- read.csv(threshold_path)
} else {
  stop("Path to threshold file must be provided. See script usage (--help)")
}
if (!is.na(opt$threshold)) {
  outliers_path <- opt$outliers
  outliers <- readLines(outliers_path, warn = FALSE)
} else {
  stop("Path to outliers file must be provided. See script usage (--help)")
}

setwd(script_dir)
source("00_functions.R")

### READ IN PHYLOSEQ OBJECTS ##################################################
# Read in metadata
print("Reading in metadata file...")
metadata <- read_csv(metadata_path)

# get sample types
types <- unique(metadata$sample_type)
print("Sample types:")
types <- types[!grepl("Negative Control", types)]
print(types)
for(i in 1:length(types)){
  print(paste0("Processing ", types[i], " samples"))
  
  setwd(paste0(pwd,"02_Clean_Data/02_DADA2_ASV_Tables/",amplicon))
  print("Reading in most recent version of phyloseq object:")
  ps <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                     "_phyloseq_", types[i], 
                                     "_raw_withnegcontrols_"), ".RDS")
  
  print("Saving negative controls")
  ps_nc <- prune_samples(sample_data(ps)$is_control == TRUE, ps)
  
  print("Removing negative controls from filtering analysis")
  ps_samples <- prune_samples(sample_data(ps)$is_control == FALSE, ps)
  
  # filter outliers out of data and metadata
  ps_filtered <- prune_samples(!(sample_names(ps) %in% outliers), ps)
  sample_metadata <- as.data.frame(as.matrix(ps_filtered@sam_data))
  
  print("Visualizing removal of outliers:")
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
               "/Figures"))
  # evaluate whether you have removed all outliers; re-run until you feel you have 
  # removed all outliers in the above lines
  id_outliers_evaluate_seq_depth(ps_filtered, sample_metadata, types[i], yourname, 
                                 amplicon, date, "outliers_removed")
  
  print("Testing drop thresholds:")
  threshold1 <- thresholds$threshold1[which(thresholds$sample_type == types[i])]
  threshold2 <- thresholds$threshold2[which(thresholds$sample_type == types[i])]
  test_drop_threshold(ps_filtered, sample_metadata, types[i], yourname, amplicon, 
                      date, threshold1)
  test_drop_threshold(ps_filtered, sample_metadata, types[i], yourname, amplicon, 
                      date, threshold2)
}

print("Done with all sample types. Check the saved figures to see if you have any remaining outliers and/or if the one of the drop thresholds you tested is okay. I'd suggest testing 2+ drop thresholds!")