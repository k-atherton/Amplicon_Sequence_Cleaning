### LOAD IN PACKAGES ##########################################################
library(vroom)
library(plyr)
library(dplyr)
library(phyloseq)

### SCRIPT SETUP ##############################################################
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"
amplicon <- "16S" # options: 16S or ITS
yourname <- "atherton" # user's last name for file storage purposes
edit_metadata <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")
setwd("02_Clean_Data")

### READ IN AND FORMAT ASV TABLES #############################################
# load sequencing data
nr1 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_NR1/"), 
                               paste0("atherton_NR1_", amplicon, 
                                      "_ASV_table_"), ".tsv")
nr2 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_NR2/"), 
                               paste0("atherton_NR2_", amplicon, 
                                      "_ASV_table_"), ".tsv")
r1 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run1/"), 
                              paste0("atherton_Run1_", amplicon, 
                                     "_ASV_table_"), ".tsv")
r2 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run2/"), 
                              paste0("atherton_Run2_", amplicon, 
                                     "_ASV_table_"), ".tsv")
r3 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run3/"), 
                              paste0("atherton_Run3_", amplicon, 
                                     "_ASV_table_"), ".tsv")
r4 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run4/"), 
                              paste0("atherton_Run4_", amplicon, 
                                     "_ASV_table_"), ".tsv")
r5 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run5/"), 
                              paste0("atherton_Run5_", amplicon, 
                                     "_ASV_table_"), ".tsv")
r6 <- read_in_dada2_asv_table(paste0("01_DADA2/", amplicon, "_Run6/"), 
                              paste0("atherton_Run6_", amplicon, 
                                     "_ASV_table_"), ".tsv")
if(amplicon == "16S"){
  rerun <- read_in_dada2_asv_table("01_DADA2/16S_ReRun/", 
                                   "atherton_ReRun_16S_ASV_table_", 
                                   ".tsv") 
}

# create dataframe of all sequence data together
if(amplicon == "16S"){
  data <- join_all(list(nr1, nr2, r1, r2, r3, r4, r5, r6, rerun), 
                   by = "ASV", type = "full")
  rm(list = c('nr1', 'nr2', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'rerun')) 
} else{
  data <- join_all(list(nr1 ,nr2, r1, r2, r3, r4, r5, r6), by = "ASV", 
                   type = "full")
  rm(list = c('nr1', 'nr2', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6'))
}

# change NA to 0 in sequence data
data[is.na(data)] <- 0

# make ASV names rownames, remove ASV column
rownames(data) <- data$ASV
data <- data[,-1]

### READ IN AND FORMAT TAXONOMY TABLES ########################################
# load taxonomic data
tax_nr1 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_NR1/"), 
                                  paste0("atherton_NR1_", amplicon, 
                                         "_taxonomy_"), ".tsv")
tax_nr2 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_NR2/"), 
                                  paste0("atherton_NR2_", amplicon, 
                                         "_taxonomy_"), ".tsv")
tax_r1 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run1/"), 
                                 paste0("atherton_Run1_", amplicon, 
                                        "_taxonomy_"), ".tsv")
tax_r2 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run2/"), 
                                 paste0("atherton_Run2_", amplicon, 
                                        "_taxonomy_"), ".tsv")
tax_r3 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run3/"), 
                                 paste0("atherton_Run3_", amplicon, 
                                        "_taxonomy_"), ".tsv")
tax_r4 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run4/"), 
                                 paste0("atherton_Run4_", amplicon, 
                                        "_taxonomy_"), ".tsv")
tax_r5 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run5/"), 
                                 paste0("atherton_Run5_", amplicon, 
                                        "_taxonomy_"), ".tsv")
tax_r6 <- read_in_dada2_taxonomy(paste0("01_DADA2/", amplicon, "_Run6/"), 
                                 paste0("atherton_Run6_", amplicon, 
                                        "_taxonomy_"), ".tsv")
if(amplicon == "16S"){
  tax_rerun <- read_in_dada2_taxonomy("01_DADA2/16S_ReRun/", 
                                      "atherton_ReRun_16S_taxonomy_", 
                                      ".tsv")
}


# create dataframe of all taxonomy information together
if(amplicon == "16S"){
  tax <- join_all(list(tax_nr1, tax_nr2, tax_r1, tax_r2, tax_r3, tax_r4, 
                       tax_r5, tax_r6, tax_rerun), 
                  by = c("Feature ID", "Taxon"), type = "full")
  rm(list = c('tax_nr1', 'tax_nr2', 'tax_r1', 'tax_r2', 'tax_r3', 'tax_r4',
              'tax_r5', 'tax_r6','tax_rerun')) 
} else{
  tax <- join_all(list(tax_nr1, tax_nr2, tax_r1, tax_r2, tax_r3, tax_r4, 
                       tax_r5, tax_r6), 
                  by = c("Feature ID", "Taxon"), type = "full")
  rm(list = c('tax_nr1', 'tax_nr2', 'tax_r1', 'tax_r2', 'tax_r3', 'tax_r4', 
              'tax_r5', 'tax_r6'))
}

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

# write files to csv
write.csv(data, paste0("02_DADA2_ASV_Tables/", amplicon, "/", yourname,
                       "_", amplicon, "_ASV_table_allsampletypes_raw",
                       date,".csv"))
write.csv(tax,paste0("02_DADA2_ASV_Tables/", amplicon, "/", yourname,
                     "_", amplicon, "_taxonomy_allsampletypes_raw",
                     date,".csv"))

### READ IN SAMPLE METADATA ###################################################
setwd(paste0(pwd,"01_Collect_Data/01_Sample_Metadata"))
# record raw dada2 read counts in metadata table
if(amplicon == "16S"){
  metadata <- read_in_file(getwd(), "atherton_sample_metadata_16S_", ".csv") 
} else{
  metadata <- read_in_file(getwd(), "atherton_sample_metadata_ITS_", ".csv") 
}

if(edit_metadata %in% c("Y", "y")){
  ### EDIT SAMPLE METADATA FILE ###############################################
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
setwd(paste0(pwd,"02_Clean_Data/02_DADA2_ASV_Tables/", amplicon))

# make phyloseq objects
if(amplicon == "ITS"){
  colnames(data) <- gsub("neg_control_1.1", "neg_control_1-1", 
                         colnames(data))
  colnames(data) <- gsub("neg_control_1.2", "neg_control_1-2", 
                         colnames(data))
}
tax <- tax_table(as.matrix(tax))
data <- otu_table(as.matrix(data), taxa_are_rows = TRUE)

ps <- phyloseq(data, tax)

meta <- as.data.frame(metadata)
meta <- meta[!is.na(meta$sample_name),]
rownames(meta) <- meta$sample_name
meta <- meta[,-1]

meta <- sample_data(meta)

ps_meta <- merge_phyloseq(ps, meta)

# make leaf dataset
if(amplicon == "16S"){
  leaf_ncs <- c("neg_control_1_16S", "neg_control_2_16S", "leaf_1", 
                "leaf_2", "leaf_3", "leaf_4", "leaf_5", "leaf_6") 
} else{
  leaf_ncs <- c("neg_control_1-1_ITS", "neg_control_1-2_ITS",
                "neg_control_2_ITS", "leaf_2", "leaf_3", "leaf_4", 
                "leaf_5", "leaf_6")
}
ps_leaf_w_nc <- subset_phyloseq(data, ps_meta, "Leaf", leaf_ncs)

# make root dataset
root_ncs <- c("1_root", "2_root", "root_1", "root_2", "root_3", "root_4", 
              "root_5", "root_6")
ps_root_w_nc <- subset_phyloseq(data, ps_meta, "Root", root_ncs)

# make M soil dataset
if(amplicon == "16S"){
  msoil_ncs <- c("neg_control_1_16S", "neg_control_2_16S", 
                 "neg_control_M_1", "neg_control_M_2", "neg_control_M_3", 
                 "neg_control_M_4", "neg_control_M_5", "neg_control_M_6")
} else{
  msoil_ncs <- c("neg_control_1-1_ITS", "neg_control_1-2_ITS", 
                 "neg_control_2_ITS", "neg_control_M_1", "neg_control_M_2", 
                 "neg_control_M_3", "neg_control_M_4", "neg_control_M_5", 
                 "neg_control_M_6")
}
ps_msoil_w_nc <- subset_phyloseq(data, ps_meta, "MSoil", msoil_ncs)

# make O soil dataset
if(amplicon == "16S"){
  osoil_ncs <- c("neg_control_1_soil", "neg_control_2_soil", 
                 "neg_control_O_1", "neg_control_O_2", "neg_control_O_3", 
                 "neg_control_O_4", "neg_control_O_5", "neg_control_O_6", 
                 "neg_control_o_rerun1_1", "neg_control_o_rerun1_2", 
                 "neg_control_o_rerun2_1", "neg_control_o_rerun2_2")
} else{
  osoil_ncs <- c("neg_control_1_soil", "neg_control_2_soil", 
                 "neg_control_O_1", "neg_control_O_2", "neg_control_O_3", 
                 "neg_control_O_4", "neg_control_O_5", "neg_control_O_6")
}
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
