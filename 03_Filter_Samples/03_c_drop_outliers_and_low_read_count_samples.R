### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(phyloseq)
library(ggplot2)

### SCRIPT SETUP ##############################################################
print("SETTING UP SCRIPT:")
date <- format(Sys.Date(),"_%Y%m%d")

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default="16S", 
              help="amplicon dataset to filter; options: 16S or ITS [default= %default]", 
              metavar="amplicon"),
  make_option(c("-n", "--name"), type="character", default="atherton", 
              help="last name for output file naming scheme [default= %default]", 
              metavar="yourname"),
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
              metavar="script_directory"),
  make_option(c("-t", "--threshold"), type="character", default=NA, 
              help="a .csv file that has two columns: sample_type (which matches the sample types (not including negative controls) in the metadata file) and threshold (defines the minimum number of reads needed for a sample to be kept for downstream analysis) [default= %default]",
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
if (!is.na(opt$metadata)) {
  metadata_path <- opt$metadata
} else {
  stop("Path to metadata file must be provided. See script usage (--help)")
}
if (!is.na(opt$threshold)) {
  threshold_path <- opt$threshold
  thresholds <- vroom::vroom(threshold_path)
} else {
  stop("Path to threshold file must be provided. See script usage (--help)")
}
if (!is.na(opt$outliers)) {
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
if(edit_metadata == "Y"){
  metadata$sequences_dropped <- "No"
  for(i in 1:nrow(metadata)){
    if(metadata$sample_data[i] %in% outliers){
      metadata$sequences_dropped <- "Outlier"
    }
  }
}

# get sample types
types <- unique(metadata$sample_type[which(metadata$is_control == FALSE)])
print("Sample types:")
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
  print("Removing outliers")
  ps_filtered <- prune_samples(!(sample_names(ps_samples) %in% outliers), 
                               ps_samples)
  sample_metadata <- as.data.frame(as.matrix(ps_filtered@sam_data))
  
  # subset the filtered samples out of the phyloseq
  threshold <- thresholds$threshold[which(thresholds$sample_type == types[i])]
  ps_final <- subset_samples(ps_filtered, seq_count_dada2 > threshold)
  
  # write to file
  ensure_directory_exists(paste0(pwd,
                                 "02_Clean_Data/03_Filter_Samples_ASV_Tables/", 
                                 amplicon))
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
  write.csv(otu_table(ps_final), paste0(yourname, "_", amplicon, 
                                        "_", types[i], "_ASV_table_filteredsamples", 
                                        date, ".csv"))
  saveRDS(ps_final, paste0(yourname, "_", amplicon, "_", types[i],
                           "_phyloseq_filteredsamples", date, ".RDS"))
  
  # merge_with_negative_controls
  ps_final_with_nc <- merge_phyloseq(ps_final, ps_nc)
  saveRDS(ps_final_with_nc, paste0(yourname, "_", amplicon, "_", types[i],
                                   "_phyloseq_filteredsamples_withnc", date, 
                                   ".RDS"))
  
  # plot histogram of post-filtering sequencing depth
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
               "/Figures/", types[i]))
  plot_filter_seq_depth(sample_metadata, types[i], threshold, date)
  
  if(edit_metadata == "Y"){
    for(i in 1:nrow(metadata)){
      if(metadata$sample_type == types[i] & metadata$seq_count_dada2 < threshold){
        metadata$sequences_dropped <- "Low read count"
      }
    }
  }
}

if(edit_metadata == "Y"){
  setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
  metadata_file_name <- basename(metadata_path)
  new_metadata_file_name <- sub("_\\d{8}\\.csv$", paste0(date, ".csv"), 
                                metadata_file_name)
  write.csv(metadata, new_metadata_file_name)
}

print("Done with all sample types.")