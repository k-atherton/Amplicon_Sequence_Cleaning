### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(phyloseq)
library(robCompositions)
library(vegan)
library(ggplot2)
library(gridExtra)
library(sva)
library(dplyr)
library(tidyr)
library(compositions)
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
  make_option(c("-b", "--batchcorrect"), type="character", default="N", 
              help="do you want to run batch correction? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-r", "--rarefy"), type="character", default="N", 
              help="do you want to rarefy your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=NA, 
              help="rarefaction threshold, default = minimum(sample sums)", 
              metavar="threshold"),
  make_option(c("-l", "--clr"), type="character", default="N", 
              help="do you want to clr-transform your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-z", "--zscore"), type="character", default="N", 
              help="do you want to z-score transform your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-d", "--aitchisondistance"), type="character", default="N", 
              help="do you want to calculate the aitchison distance matrix for your samples? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-v", "--variables"), type="character", default=NA, 
              help="text file with list of varibles to color data by for NMDS plots",
              metavar="variables"),
  make_option(c("--batchvariable"), type="character", default="sequencing_batch",
              help="variable to use for batch correction [default= %default]",
              metavar="character"),
  make_option(c("--covariates"), type="character", default=NA,
              help="path to text file containing variable(s) to maintain the signature of during batch correction [default= %default",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

amplicon <- opt$amplicon
yourname <- opt$name
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
bc <- opt$batchcorrect
rarefy_transform <- opt$rarefy
clr_transform <- opt$clr
zscore_transform <- opt$zscore
aitchison_transform <- opt$aitchisondistance
if(bc == "Y"){
  print("Batch correction will happen! Reading in necessary files.")
  if(!is.na(opt$variables)) {
    variables_path <- opt$variables
    color_vars <- readLines(variables_path, warn = FALSE)
  } else {
    stop("Path to variables file must be provided to run batch control. See script usage (--help)")
  }
  batch_var <- opt$batchvariable
  if(is.na(batch_var)){
    stop("A variable for batch control must be provided to run batch control. See script usage (--help)")
  }
  if (!is.na(opt$covariates)) {
    covars_path <- opt$covariates
    covars <- readLines(covars_path, warn = FALSE)
  } else{
    print("Warning: path to covariates file was not provided. Batch correction will not maintain signature of any variables. See script usage (--help)")
  }
}
if(rarefy_transform == "Y"){
  print("Rarefaction will happen! Reading in the rarefaction threshold.")
  threshold <- opt$threshold
}

setwd(script_dir)
source("00_functions.R")

### TRANSFORMATIONS ###########################################################
# Read in metadata
print("Reading in metadata file...")
metadata <- read_csv(metadata_path)

# Get sample types
types <- unique(metadata$sample_type[which(metadata$is_control == FALSE)])
print("Sample types:")
print(types)

if(bc == "Y"){
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon))

  print("Read in all data types phyloseq object:")
  setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
  ps_all <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                         "_alldatatypes_decontam_"), ".RDS")
  
  print("CHECKING FOR BATCH EFFECT:")
  check_batch_effect(ps_all, amplicon, yourname, date, color_vars)
  print("RUNNING BATCH CORRECTION: ")
  corrected <- batch_correct(ps_all, amplicon, yourname, date, batch_var, 
                             covars)
  setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon))
  print("Saving batch corrected ASV table.")
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon, "/Batch_Corrected_ASV_Tables"))
  write.csv(corrected, paste0("Batch_Corrected_ASV_Tables/", yourname, "_", 
                              amplicon, "_batch_corrected_ASV_table", date, 
                              ".csv"))
}

if(rarefy_transform == "Y"){
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon, "/Rarefied_ASV_Tables"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                                 "/Rarefied_ASV_Tables/Figures"))
  # rarefy separately for each sample type
  for(i in 1:length(types)){
    print(paste0("Reading in ", types[i], " data:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps <- read_in_file(getwd(), paste0(yourname, "_", amplicon, "_", types[i],
                                       "_decontam_"), ".RDS")
    print("RAREFYING DATA: ")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/Rarefied_ASV_Tables"))
    rarefy_data(ps, types[i], amplicon, yourname, date, threshold)
  }
}


# clr-transform separately for each sample type and all together
if(clr_transform == "Y"){
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon, "/CLR_Transformed_Tables"))
  # rarefy separately for each sample type
  for(i in 1:length(types)){
    print(paste0("Reading in ", types[i], " data:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps <- read_in_file(getwd(), paste0(yourname, "_", amplicon, "_", types[i],
                                       "_decontam_"), ".RDS")
    print("CLR-TRANSFORMING DATA: ")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/CLR_Transformed_Tables"))
    clr_data(ps, types[i], amplicon, yourname, date)
  }
  if(length(types) > 1){
    print(paste0("Reading in all data types:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps_all <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                           "_alldatatypes_decontam_"), ".RDS")
    print("CLR-TRANSFORMING DATA: ")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/CLR_Transformed_Tables"))
    clr_data(ps_all, "allsampletypes", amplicon, yourname, date)
  }
}

# zscore-transform separately for each sample type and all together
if(zscore_transform == "Y"){
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon, "/ZScore_Transformed_Tables"))
  for(i in 1:length(types)){
    print(paste0("Reading in ", types[i], " data:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps <- read_in_file(getwd(), paste0(yourname, "_", amplicon, "_", types[i],
                                       "_decontam_"), ".RDS")
    print("Z-Score-TRANSFORMING DATA: ")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/ZScore_Transformed_Tables"))
    zscore_data(ps, types[i], amplicon, yourname, date)
  }
  if(length(types) > 1){
    print(paste0("Reading in all data types:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps_all <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                           "_alldatatypes_decontam_"), ".RDS")
    print("Z-Score-TRANSFORMING DATA: ")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/ZScore_Transformed_Tables"))
    zscore_data(ps_all, "allsampletypes", amplicon, yourname, date)
  }
}

# calculate aitchison_distance sparately for each sample type and all together
# for the all together distance matrix, use the batch-corrected data
if(aitchison_transform == "Y"){
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/05_Transform_Data/", 
                                 amplicon, "/Aitchison_Distance_Tables"))
  for(i in 1:length(types)){
    print(paste0("Reading in ", types[i], " data:"))
    setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
    ps <- read_in_file(getwd(), paste0(yourname, "_", amplicon, "_", types[i],
                                       "_decontam_"), ".RDS")
    print("CALCULATING AITCHISON DISTANCE MATRIX FOR DATA:")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/Aitchison_Distance_Tables"))
    aitchison_data(ps, types[i], amplicon, yourname, date)
  }
  if(length(types) > 1 & bc == "Y"){
    print(paste0("Reading in all data types:"))
    print("CALCULATING AITCHISON DISTANCE MATRIX FOR DATA:")
    setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
                 "/Aitchison_Distance_Tables"))
    aitchison_data(corrected, "allsampletypes", amplicon, yourname, date)
  }
}