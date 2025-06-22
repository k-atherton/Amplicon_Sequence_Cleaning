### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(phyloseq)
library(ggplot2)
library(decontam)
library(readr)
library(gridExtra)

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
  make_option(c("-e", "--edit"), type="character", default="N", 
              help="do you want to edit the metadata file? options: Y or N [default= %default]", 
              metavar="edit_metadata"),
  make_option(c("-p", "--pwd"), type="character", default=getwd(),
              help="the directory for saving the outputs of this script [default= %default]",
              metavar="working_directory"),
  make_option(c("-m", "--metadata"), type="character", default=NA,
              help="path to sample metadata file", metavar="metadata"),
  make_option(c("-s","--scriptdir"), type="character", 
              default="/projectnb/talbot-lab-data/Katies_data/Amplicon_Sequence_Cleaning/",
              help="path to function script directory, [default = %default]",
              metavar="script_directory")
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

setwd(script_dir)
source("00_functions.R")

### Separate the sequence data by batch #######################################
# Read in metadata
print("Reading in metadata file...")
metadata <- read_csv(metadata_path)

# Get sample types
types <- unique(metadata$sample_type[which(metadata$is_control == FALSE)])
print("Sample types:")
print(types)

# Loop through sample types
for(i in 1:length(types)){
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
  print(paste0("Decontaminating ", types[i], " samples."))
  print(paste0("Reading in ", types[i], 
               " phyloseq object with negative controls:"))
  
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
  ps <- read_in_file(getwd(), paste0("atherton_", amplicon, "_", types[i], 
                                     "_phyloseq_filteredsamples_withnc_"), 
                     ".RDS")
  
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", 
                                 amplicon))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", 
                                 amplicon, "/Figures/"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", 
                                 amplicon, "/Figures/", types[i]))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", 
                                 amplicon, "/Contaminant_Taxonomy/"))
  ensure_directory_exists(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", 
                                 amplicon, "/Contaminant_Taxonomy/", types[i]))
  setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon))
  
  print("Running Decontam Function:")
  decontam_type <- decontaminate_samples(ps, types[i], yourname, amplicon, date)
  saveRDS(decontam_type, paste0(yourname, "_", amplicon, "_", types[i], 
                                "_decontam", date, ".RDS"))
  if(i > 1){
    print("Merging with last sample type data:")
    decontam_data <- merge_phyloseq(decontam_data, decontam_type)
  } else{
    decontam_data <- decontam_type
  }
}

if(i > 1){
  saveRDS(decontam_data, paste0(yourname, "_", amplicon, "_alldatatypes_decontam", 
                                date, ".RDS")) 
}

if(edit_metadata == "Y") {
  print("Writing decontam sequence count to metadata:")
  # add decontam sequence count column
  metadata$seq_count_decontam <- NA
  
  # sum the post-dada2 processing sequence counts in each sample
  decontam_seq_count <- colSums(decontam_data@otu_table)
  
  # record sequence counts in metadata table
  for (i in 1:length(decontam_seq_count)) {
    sample_name <- names(decontam_seq_count)[i]
    metadata$seq_count_decontam[which(metadata$sample_name == sample_name)] <-
      decontam_seq_count[i]
  }
  
  # write this data to the file
  new_metadata_file_name <- sub("_\\d{8}\\.csv$", paste0("_decontam", date, 
                                                         ".csv"), metadata_path)
  write.csv(metadata, new_metadata_file_name, row.names = FALSE)
  
  print("Visualize how decontam impacted sequencing depth:")
  setwd(paste0(pwd, "02_Clean_Data/04_Decontaminate_Samples/", amplicon, 
               "/Figures"))
  for(i in 1:length(types)){
    plot_decontam_seq_depth(metadata, types[i], date) 
  }
}
