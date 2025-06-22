### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(geosphere)
library(vegan)
library(dplyr)
library(tidyverse)
library(phyloseq)

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
  make_option(c("-e", "--edit"), type = "character", default = "N",
              help="edit metadata file? options: Y or N [default = %default]",
              metavar="edit"),
  make_option(c("-r", "--raredata"), type="character", default=NA,
              help="path to folder where rarefied data for all sample types is stored", 
              metavar="raredata"),
  make_option(c("-t", "--taxonomy"), type="character", default=NA,
              help="path to taxonomy file", metavar="taxonomy"),
  make_option(c("--pathogen16s"), type = "character", default = NA,
              help="path to MBPD BLAST results", metavar = "pathogens")
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
edit_metadata <- opt$edit
if (!is.na(opt$raredata)) {
  rarefied_data_path <- opt$raredata
} else {
  stop("Path to rarefied data file must be provided. See script usage (--help)")
}
if (!is.na(opt$taxonomy)) {
  taxonomy_path <- opt$taxonomy
} else {
  stop("Path to taxonomy file must be provided. See script usage (--help)")
}
if(amplicon == "16S"){
  if(!is.na(opt$pathogens)) {
    pathogens_path <- opt$pathogens
  } else {
    stop("Path to MBPD BLAST results file must be provided. See script usage (--help)")
  } 
}

setwd(script_dir)
source("00_functions.R")

### READ IN DATA ##############################################################
# Read in metadata
print("Reading in metadata file...")
metadata <- read_csv(metadata_path)

# Collect sample types
sample_types <- unique(metadata$sample_type)
sample_types <- sample_types[which(sample_types != "Negative Control")]

# Read in rarefied sequence data
print("Reading in rarefied sequence data files...")
seq_data <- vector(mode = "list", length = length(sample_types))
for(i in 1:length(seq_data)){
  print(paste0("Reading in ", sample_types[i], "sequence data."))
  seq_data_ps <- read_in_file(rarefied_data_path, 
                              paste0(yourname, "_", amplicon, "_", 
                                     sample_types[i], "_rarefied_phyloseq_"), 
                              ".RDS")
  seq_data[[i]] <- as.data.frame(seq_data_ps@otu_table)
}

# Read in taxonomy file
print("Reading in taxonomy file...")
taxonomy <- read_csv(taxonomy_path)
colnames(taxonomy)[1] <- "ASV_ID"

### TRANSFORM DATA TO PERCENTAGES #############################################
print("Transforming data to percentages:")
perc_data <- vector(mode = "list", length = length(sample_types))
perc_data_guild <- vector(mode = "list", length = length(sample_types))
for(i in 1:length(perc_data_guild)){
  # find the sums of the columns
  print("Finding column sums:")
  sums <- colSums(seq_data[[i]]) 
  
  # divide the individual counts by the mean of the sums (because the data was 
  # rarefied, all sums are the same) and multiply by 100 in order to get the 
  # percent abundance (relative abundance)
  print("Transforming to relative abundance:")
  perc_data[[i]] <- seq_data[[i]]/mean(sums) * 100
  
  # aggregate by guild
  print("Aggregating counts by guild:")
  if(amplicon == "16S"){
    perc_data_guild[[i]] <- aggregate_guilds_16S(seq_data[[i]], perc_data[[i]], 
                                                 taxonomy)
  } else if(amplicon == "ITS"){
    perc_data_guild[[i]] <- aggregate_guilds_ITS(seq_data[[i]], perc_data[[i]], 
                                                 taxonomy)
  }
  
}

if(amplicon == "16S"){
  guild_abund <- vector(mode = "list", length = length(sample_types))
  for(i in 1:length(perc_data_guild)){
    print("Recording guild abundances:")
    guild_abund[[i]] <- record_16s_guild_abund(perc_data_guild[[i]])
  }
  if(length(perc_data_guild) > 1){
    result <- bind_rows(perc_data_guild, .id = "column_label")
  } else{
    result <- guild_abund[[1]]
  }
  result$sample_name <- rownames(result)
  
  print("Saving relative abundances to file.")
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function/16S"))
  write.csv(result, 
            paste0(pwd, 
                   "02_Clean_Data/06_Calculate_Diversity_and_Function/16S/",
                   yourname, "_", amplicon, "_guildrelativeabundances", date, 
                   ".csv"))
  
  if(edit_metadata == "Y"){
    print("Merge metadata with the guild relative abundances.")
    metadata <- merge(metadata, result, by = "sample_name", all.x = "TRUE")
    write.csv(metadata, 
              paste0(pwd, "01_Collect_Data/01_Sample_Metadata/", yourname, "_",
                     amplicon, 
                     "_sample_metadata_dada2_dropsamples_decontam_guildabund",
                     date, ".csv"))
  }
} else if(amplicon == "ITS"){
  guild_abund <- vector(mode = "list", length = length(sample_types))
  for(i in 1:length(perc_data_guild)){
    print("Recording guild abundances:")
    guild_abund[[i]] <- record_its_guild_abund(perc_data_guild[[i]], 
                                               perc_data[[i]], taxonomy)
  }
  if(length(guild_abund) > 1){
    result <- bind_rows(guild_abund, .id = "column_label")
  } else{
    result <- guild_abund[[1]]
  }
  result$sample_name <- rownames(result)
  
  print("Saving relative abundances to file.")
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function/ITS"))
  write.csv(result, 
            paste0(pwd, 
                   "02_Clean_Data/06_Calculate_Diversity_and_Function/ITS/",
                   yourname, "_", amplicon, "_guildrelativeabundances", date, 
                   ".csv"))
  
  if(edit_metadata == "Y"){
    print("Merge metadata with the guild relative abundances.")
    metadata <- merge(metadata, result, by = "sample_name", all.x = T)
    write.csv(metadata, 
              paste0(pwd, "01_Collect_Data/01_Sample_Metadata/", yourname, "_",
                     amplicon, 
                     "_sample_metadata_dada2_dropsamples_decontam_guildabund",
                     date, ".csv"))
  }
}

### PATHOGENS #################################################################
if(amplicon == "16S"){
  print("Reading in MBPD BLAST results.")
  pathogens <- vroom::vroom(pathogens_path)
  pathogens$row <- rownames(pathogens)
  
  print("Reading in pathogens database taxonomy file.")
  pathogens_tax <- vroom::vroom(paste0(script_dir, 
                                       "00_Databases/pathogen_tax.txt"))
  
  taxonomy <- process_16s_pathogens(pathogens_file, pathogens_tax, taxonomy)
  
  # Aggregate counts by pathogen type
  perc_data_pathogen <- vector(mode = "list", length = length(sample_types))
  for(i in 1:length(perc_data_guild)){
    # aggregate by guild
    print("Aggregating counts by pathogen type:")
    perc_data_pathogen[[i]] <- aggregate_16s_pathogens(seq_data[[i]], 
                                                       perc_data[[i]], taxonomy)
  }
  
  pathogen_abund <- vector(mode = "list", length = length(sample_types))
  for(i in 1:length(perc_data_pathogen)){
    print("Recording pathogen abundances:")
    pathogen_abund[[i]] <- record_16s_pathogen_abund(perc_data_pathogen[[i]])
  }
  if(length(perc_data_pathogen) > 1){
    result <- bind_rows(perc_data_pathogen, .id = "column_label")
  } else{
    result <- pathogen_abund[[1]]
  }
  result$sample_name <- rownames(result)
  
  print("Saving relative abundances to file.")
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/06_Calculate_Diversity_and_Function/16S"))
  write.csv(result, 
            paste0(pwd, 
                   "02_Clean_Data/06_Calculate_Diversity_and_Function/16S/",
                   yourname, "_", amplicon, "_pathogenrelativeabundances", date, 
                   ".csv"))
  
  if(edit_metadata == "Y"){
    print("Merge metadata with the pathogen relative abundances.")
    metadata <- merge(metadata, result, by = "sample_name", all.x = "TRUE")
    write.csv(metadata, 
              paste0(pwd, "01_Collect_Data/01_Sample_Metadata/", yourname, "_",
                     amplicon, 
                     "_sample_metadata_dada2_dropsamples_decontam_guildabund_pathogenabund",
                     date, ".csv"))
  }
  
}

### CALCULATE DIVERSITY #######################################################
diversity <- vector(mode = "list", length = length(sample_types))
print("Calculating alpha diversity of the samples.")
for(i in 1:length(sample_types)){
  print(paste0(sample_types[[i]], " sample diversity:"))
  # transpose the data
  seq_t <- t(seq_data[[i]])
  
  print("Calculating shannon diversity.")
  shannon <- diversity(seq_t, index="shannon",base=2)
  
  print("Calculating simpson diverstiy.")
  simpson <- diversity(seq_t, index="simpson")
  
  print("Calculating inverse simpson diversity.")
  invsimpson <- diversity(seq_t, index="invsimpson")
  
  print("Calculating fisher diversity.")
  fisher <- fisher.alpha(seq_t)
  
  # merging all diversity metrics together
  diversity[[i]] <- cbind(shannon, simpson, invsimpson, fisher)
}

# make one diversity dataframe 
if(length(sample_types) > 1){
  diversity <- bind_rows(diversity, .id = "column_label")
} else{
  diversity <- as.data.frame(diversity[[i]])
}

# make sample name a column for merging 
diversity$sample_name <- rownames(diversity)

print("Merging diversity metrics with metadata.")
metadata <- merge(metadata, diversity, by = "sample_name", all.x = T)

print("Saving metadata with diversity metrics.")
if(amplicon == "16S"){
  write.csv(metadata, 
            paste0(pwd, "01_Collect_Data/01_Sample_Metadata/", yourname, "_",
                   amplicon, 
                   "_sample_metadata_dada2_dropsamples_decontam_guildabund_pathogenabund_diversity",
                   date, ".csv"))
} else{
  write.csv(metadata, 
            paste0(pwd, "01_Collect_Data/01_Sample_Metadata/", yourname, "_",
                   amplicon, 
                   "_sample_metadata_dada2_dropsamples_decontam_guildabund_diversity",
                   date, ".csv"))
}
