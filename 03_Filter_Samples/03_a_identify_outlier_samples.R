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
  make_option(c("-v", "--variables"), type="character", default=NA, 
              help="text file with list of varibles to color data by for NMDS plots",
              metavar="variables")
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
if (!is.na(opt$metadata)) {
  metadata_path <- opt$metadata
} else {
  stop("Path to metadata file must be provided. See script usage (--help)")
}
if (!is.na(opt$variables)) {
  variables_path <- opt$variables
  color_vars <- readLines(variables_path, warn = FALSE)
  if(!"seq_bin" %in% color_vars){
    color_vars <- c(color_vars, "seq_bin")
  }
} else {
  stop("Path to variables file must be provided. See script usage (--help)")
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
  
  print("Saving metadata in an object")
  sample_metadata <- as.data.frame(as.matrix(ps_samples@sam_data))
  
  # make sequencing depth sample_metadata
  sample_metadata$seq_count_data2 <- as.numeric(sample_metadata$seq_count_dada2)
  
  print("Binning sequencing depth")
  sample_metadata <- bin_seq_depth(sample_metadata)
  
  print("Visualizing raw data")
  ensure_directory_exists(paste0(pwd,
                                 "02_Clean_Data/03_Filter_Samples_ASV_Tables/", 
                                 amplicon, "/Figures"))
  ensure_directory_exists(paste0(pwd,
                                 "02_Clean_Data/03_Filter_Samples_ASV_Tables/", 
                                 amplicon, "/Figures/",types[i]))
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
               "/Figures"))
  plot_prefilter_seq_depth(sample_metadata, types[i], 8000, date)
  
  print("Looking for outliers:")
  ensure_directory_exists(paste0(pwd,
                                 "02_Clean_Data/03_Filter_Samples_ASV_Tables/", 
                                 amplicon, "/Figures/", types[i]))
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
               "/Figures/", types[i]))
  
  otu_raw <- otu_table(ps_samples)
  id_outliers_evaluate_seq_depth(otu_raw, sample_metadata, types[i], yourname, 
                                 amplicon, date, "no_outliers_removed", 
                                 color_vars)
}

print("Done with all sample types. Check the saved figures to identify outliers.")