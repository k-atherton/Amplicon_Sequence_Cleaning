### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(dplyr)
library(phyloseq)
library(readxl)
library(purrr)
library(readr)

### SCRIPT SETUP ##############################################################
print("SETTING UP SCRIPT:")
date <- format(Sys.Date(),"_%Y%m%d")

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default="16S", 
              help="amplicon dataset to filter; options: 16S or ITS [default= %default]", 
              metavar="character"),
  make_option(c("-n", "--name"), type="character", default="atherton", 
              help="last name for output file naming scheme [default= %default]", 
              metavar="character"),
  make_option(c("-e", "--edit"), type="character", default="N", 
              help="do you want to edit the metadata file? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-p", "--pwd"), type="character", default=getwd(),
              help="the directory for saving the outputs of this script [default= %default]",
              metavar="character"),
  make_option(c("-v", "--asvtable"), type="character", default=NA,
              help="file with paths to ASV tables output by DADA2",
              metavar="character"),
  make_option(c("-t", "--taxonomy"), type="character", default=NA,
              help="file with paths to taxonomy tables output by DADA2",
              metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NA,
              help="path to sample metadata file", metavar="character"),
  make_option(c("-c", "--negcontrols"), type="character", default=NA,
              help="path to negative controls naming scheme file", 
              metavar = "character"),
  make_option(c("-s","--scriptdir"), type="character", 
              default="/projectnb/talbot-lab-data/Katies_data/Amplicon_Sequence_Cleaning/",
              help="path to function script directory, [default = %default]",
              metavar="character")
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

if (!is.na(opt$asvtable)) {
  asv_table_paths <- opt$asvtable
} else {
  stop("File with ASV table paths must be provided. See script usage (--help)")
}
if (!is.na(opt$taxonomy)) {
  taxonomy_paths <- opt$taxonomy
} else {
  stop("File with taxonomy table paths must be provided. See script usage (--help)")
}
if (!is.na(opt$metadata)) {
  metadata_path <- opt$metadata
} else {
  stop("Path to metadata file must be provided. See script usage (--help)")
}
if (!is.na(opt$negcontrols)) {
  negcontrols_path <- opt$negcontrols
} else {
  stop("Path to metadata file must be provided. See script usage (--help)")
}

setwd(script_dir)
source("00_functions.R")
ensure_directory_exists(paste0(pwd,"02_Clean_Data"))
setwd(paste0(pwd, "02_Clean_Data"))

### READ IN AND FORMAT ASV TABLES #############################################
print("READING IN AND FORMATTING ASV TABLES:")
# load sequencing data
asv_tables <- load_files(asv_table_paths)

# change first column to "ASV" for merging
for(i in 1:length(asv_tables)){
  colnames(asv_tables[[i]])[1] <- "ASV"  
}

# create dataframe of all sequence data together
data <- purrr::reduce(asv_tables, full_join, by = "ASV")

# change NA to 0 in sequence data
data[is.na(data)] <- 0

# make ASV names rownames, remove ASV column
rownames(data) <- data$ASV
data <- data[,-1]

### READ IN AND FORMAT TAXONOMY TABLES ########################################
print("READING IN AND FORMATTING TAXONOMY TABLES:")
# load taxonomic data
taxonomy_tables <- load_files(taxonomy_paths)

# keep only the taxonomic information
for(i in 1:length(taxonomy_tables)){
  taxonomy_tables[[i]] <- taxonomy_tables[[i]][,c(1,2)]
}

# create dataframe of all taxonomy information together
tax <- purrr::reduce(taxonomy_tables, full_join, by = c("Feature.ID", "Taxon"))

# separate the taxonomy information into ranks
names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
           "Species")
split_taxa <- stringr::str_split(tax$Taxon, pattern = ";")
taxa_names <- lapply(split_taxa, function(x) x[1:length(names)])
tax_table <- taxa_names %>%
  unlist() %>% matrix(ncol = length(names), byrow = TRUE) %>%
  data.frame() %>%
  magrittr::set_colnames(names) %>%
  mutate(Genus = replace(Genus, 
                         Genus == "Escherichia-Shigella", "Escherichia"))

# add taxonomic ranks back to taxon table, 
# remove column with ranks as one string
tax <- cbind(tax, tax_table)
tax <- tax[,-2]

# make Feature ID the rownames, remove column
rownames(tax) <- tax[,1]
tax <- tax[,-1]
rm(list = c('asv_tables', 'taxonomy_tables', 'names','split_taxa','taxa_names',
            'tax_table'))

# add functional guild to taxonomy
if(amplicon == "16S"){
  functional_guilds <- read_in_file(paste0(script_dir, "00_Databases/"), 
                                    "werbin_bacteria_functional_groups", ".csv")
  tax <- add_16s_guild(tax, functional_guilds)
} else{
  functional_guilds <- read_excel(paste0(script_dir, 
                                         "00_Databases/fungal_traits_database.xlsx"))
  tax <- add_its_guild(tax, functional_guilds)
}

# write files to csv
ensure_directory_exists(paste0(pwd,"02_Clean_Data/02_DADA2_ASV_Tables"))
ensure_directory_exists(paste0(pwd,"02_Clean_Data/02_DADA2_ASV_Tables/", 
                               amplicon))
setwd(paste0(pwd, "02_Clean_Data/02_DADA2_ASV_Tables/", amplicon))
write.csv(data, paste0(yourname, "_", amplicon, "_ASV_table_allsampletypes_raw",
                       date,".csv"))
write.csv(tax,paste0(yourname, "_", amplicon, "_taxonomy_allsampletypes_raw",
                     date,".csv"))

### READ IN SAMPLE METADATA ###################################################
print("READING IN SAMPLE METADATA:")
ensure_directory_exists(paste0(pwd,"01_Collect_Data/01_Sample_Metadata/"))
setwd(paste0(pwd,"01_Collect_Data/01_Sample_Metadata/"))
# record raw dada2 read counts in metadata table
metadata <- read_csv(metadata_path)

# remove samples not in metadata
data <- data[ ,colnames(data) %in% metadata$sample_name]

if(edit_metadata %in% c("Y", "y")){
  ### EDIT SAMPLE METADATA FILE ###############################################
  print("EDITING SAMPLE METADATA FILE:")
  # add is.control column
  metadata$is_control <- FALSE
  # record negatvie controls as controls in is.control column
  metadata$is_control[which(metadata$sample_type == "Negative Control")] <- 
    TRUE
  
  # sum the post-dada2 processing sequence counts in each sample
  dada2_seq_count <- colSums(data)
  
  # record sequence counts in metadata table
  metadata$seq_count_dada2 <- NA
  for(i in 1:length(dada2_seq_count)){
    sample_name <- names(dada2_seq_count)[i]
    metadata$seq_count_dada2[which(metadata$sample_name == sample_name)] <- 
      dada2_seq_count[i]
  }
  
  # write this data to the file
  write.csv(metadata, paste0(getwd(), "/", yourname, "_", amplicon, 
                             "_sample_metadata", date, ".csv"),
            row.names = FALSE)
}

### SEPARATE ASV TABLE BY SAMPLE TYPE #########################################
print("SEPARATING ASV TABLE BY SAMPLE TYPE:")
setwd(paste0(pwd,"02_Clean_Data/02_DADA2_ASV_Tables/", amplicon))
# 
# # make phyloseq objects
# if(amplicon == "ITS"){
#   colnames(data) <- gsub("neg_control_1.1", "neg_control_1-1", 
#                          colnames(data))
#   colnames(data) <- gsub("neg_control_1.2", "neg_control_1-2", 
#                          colnames(data))
# }
tax <- tax_table(as.matrix(tax))
data <- otu_table(as.matrix(data), taxa_are_rows = TRUE)

ps <- phyloseq(data, tax)

meta <- as.data.frame(metadata)
meta <- meta[!is.na(meta$sample_name),]
rownames(meta) <- meta$sample_name
meta <- meta[,-1]

meta <- sample_data(meta)

ps_meta <- merge_phyloseq(ps, meta)

print("Getting negative control names.")
negative_control_names <- read_csv(negcontrols_path)
sample_types <- unique(negative_control_names$sample_type)

for(i in 1:length(sample_types)){
  type <- sample_types[i]
  print(paste0("Working with ", type, " samples."))
  print("Getting negative controls and merging with samples.")
  ncs <- negative_control_names$negative_control_name_pattern[which(negative_control_names$sample_type == type)]
  ps_w_nc <- subset_phyloseq(data, ps_meta, type, ncs)
  
  print("Making raw data tables.")
  raw <- as.data.frame(as.matrix(ps_w_nc@otu_table))
  raw_tax <- as.data.frame(as.matrix(ps_w_nc@tax_table))
  
  print("Writing to files.")
  write.csv(raw, paste0(getwd(), "/", yourname, "_", amplicon, "_ASV_table_", 
                        type, "_raw", date, ".csv"))
  write.csv(raw_tax, paste0(getwd(), "/", yourname, "_", amplicon, "_taxonomy_", 
                            type, "_raw", date, ".csv"))
  saveRDS(ps_w_nc, paste0(getwd(), "/", yourname, "_", amplicon, "_phyloseq_", 
                          type, "raw_withnegcontrols", date, ".RDS"))
}
