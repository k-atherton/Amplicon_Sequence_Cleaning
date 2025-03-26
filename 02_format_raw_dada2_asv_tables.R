### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
library(vroom)
library(plyr)
library(dplyr)
library(phyloseq)
library(readxl)

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
  make_option(c("-c", "-negcontrols"), type="character", default=NA,
              help="path to negative controls naming scheme file", 
              metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

amplicon <- opt$amplicon
yourname <- opt$name
edit_metadata <- opt$edit
pwd <- opt$pwd
script_dir <- getwd()

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

setwd(pwd)
source("00_functions.R")
ensure_directory_exists(paste0(pwd,"02_Clean_Data"))
setwd(paste0(pwd, "02_Clean_Data"))

### READ IN AND FORMAT ASV TABLES #############################################
print("READING IN AND FORMATTING ASV TABLES:")
# load sequencing data
asv_tables <- load_files(asv_table_paths)

# create dataframe of all sequence data together
data <- join_all(asv_tables, by = "ASV", type = "full")

# change NA to 0 in sequence data
data[is.na(data)] <- 0

# make ASV names rownames, remove ASV column
rownames(data) <- data$ASV
data <- data[,-1]

### READ IN AND FORMAT TAXONOMY TABLES ########################################
print("READING IN AND FORMATTING TAXONOMY TABLES:")
# load taxonomic data
taxonomy_tables <- load_files(taxonomy_paths)

# create dataframe of all taxonomy information together
tax <- join_all(taxonomy_tables, by = c("Feature ID", "Taxon"), type = "full")

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
rm(list = c('names','split_taxa','taxa_names','tax_table'))

# add functional guild to taxonomy
if(amplicon == "16S"){
  functional_guilds <- read_in_file(script_dir, 
                                    "werbin_bacteria_functional_groups", ".csv")
  tax <- add_16s_guild(tax, functional_guilds)
} else{
  functional_guilds <- read_excel(paste0(script_dir, 
                                         "fungal_traits_database.xlsx"))
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
setwd(paste0(pwd,"01_Collect_Data/01_Sample_Metadata"))
# record raw dada2 read counts in metadata table
metadata <- read.csv(metadata_path)

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
  for(i in 1:length(dada2_seq_count)){
    sample_name <- names(dada2_seq_count)[i]
    metadata$seq_count_dada2[which(metadata$sample_name == sample_name)] <- 
      dada2_seq_count[i]
  }
  
  # write this data to the file
  write.csv(metadata, paste0(getwd(), yourname, "_sample_metadata_", 
                             amplicon, date, ".csv"),
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

negative_control_names <- read.csv(negcontrols_path)

# make leaf dataset
leaf_ncs <- negative_control_names[which(negative_control_names$sample_type == "leaf"),]
ps_leaf_w_nc <- subset_phyloseq(data, ps_meta, "Leaf", leaf_ncs)

# make root dataset
root_ncs <- negative_control_names[which(negative_control_names$sample_type == "root"),]
ps_root_w_nc <- subset_phyloseq(data, ps_meta, "Root", root_ncs)

# make M soil dataset
msoil_ncs <- negative_control_names[which(negative_control_names$sample_type == "msoil"),]
ps_msoil_w_nc <- subset_phyloseq(data, ps_meta, "MSoil", msoil_ncs)

# make O soil dataset
osoil_ncs <- negative_control_names[which(negative_control_names$sample_type == "osoil"),]
ps_osoil_w_nc <- subset_phyloseq(data, ps_meta, "OSoil", osoil_ncs)

# save as data frames
leaf_raw <- as.data.frame(as.matrix(ps_leaf_w_nc@otu_table))
leaf_raw_tax <- as.data.frame(as.matrix(ps_leaf_w_nc@tax_table))
root_raw <- as.data.frame(as.matrix(ps_root_w_nc@otu_table))
root_raw_tax <- as.data.frame(as.matrix(ps_root_w_nc@tax_table))
msoil_raw <- as.data.frame(as.matrix(ps_msoil_w_nc@otu_table))
msoil_raw_tax <- as.data.frame(as.matrix(ps_msoil_w_nc@tax_table))
osoil_raw <- as.data.frame(as.matrix(ps_osoil_w_nc@otu_table))
osoil_raw_tax <- as.data.frame(as.matrix(ps_osoil_w_nc@tax_table))

# write to SCC
write.csv(leaf_raw, paste0(getwd(), "/", yourname, "_", amplicon, 
                           "_ASV_table_leaf_raw", date, ".csv"))
write.csv(root_raw, paste0(getwd(), "/", yourname, "_", amplicon, 
                           "_ASV_table_root_raw", date, ".csv"))
write.csv(msoil_raw, paste0(getwd(), "/", yourname, "_", amplicon, 
                            "_ASV_table_msoil_raw", date, ".csv"))
write.csv(osoil_raw, paste0(getwd(), "/", yourname, "_", amplicon, 
                            "_ASV_table_osoil_raw", date, ".csv"))
write.csv(leaf_raw_tax, paste0(getwd(), "/", yourname, "_", amplicon, 
                               "_taxonomy_leaf_raw", date, ".csv"))
write.csv(root_raw_tax, paste0(getwd(), "/", yourname, "_", amplicon, 
                               "_taxonomy_root_raw", date, ".csv"))
write.csv(msoil_raw_tax, paste0(getwd(), "/", yourname, "_", amplicon, 
                                "_taxonomy_msoil_raw", date, ".csv"))
write.csv(osoil_raw_tax, paste0(getwd(), "/", yourname, "_", amplicon, 
                                "_taxonomy_osoil_raw", date, ".csv"))

saveRDS(ps_leaf_w_nc, paste0(getwd(), "/", yourname, "_", amplicon, 
                             "_phyloseq_leaf_raw_withnegcontrols", date, 
                             ".RDS"))
saveRDS(ps_root_w_nc, paste0(getwd(), "/", yourname, "_", amplicon, 
                             "_phyloseq_root_raw_withnegcontrols", date, 
                             ".RDS"))
saveRDS(ps_msoil_w_nc, paste0(getwd(), "/", yourname, "_", amplicon, 
                              "_phyloseq_msoil_raw_withnegcontrols", date, 
                              ".RDS"))
saveRDS(ps_osoil_w_nc, paste0(getwd(), "/", yourname, "_", amplicon, 
                              "_phyloseq_osoil_raw_withnegcontrols", date, 
                              ".RDS"))
