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
              metavar="variables"),
  make_option(c("-e", "--edit"), type="character", default="N", 
              help="do you want to edit the metadata file? options: Y or N [default= %default]", 
              metavar="edit_metadata")
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