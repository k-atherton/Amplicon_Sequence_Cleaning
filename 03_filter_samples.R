### LOAD IN PACKAGES ##########################################################
library(vroom)
library(phyloseq)
library(ggplot2)
library(robCompositions)
library(gridExtra)
library(vegan)

### SCRIPT SETUP ##############################################################
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"
amplicon <- "ITS" # options: 16S or ITS
yourname <- "atherton" # user's last name for file storage purposes
edit_metadata <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")

### READ IN FORMATTED ASV TABLES ###########################################
setwd(paste0("02_Clean_Data/02_DADA2_ASV_Tables/",amplicon))
ps_leaf <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_phyloseq_leaf_raw_withnegcontrols_"), 
                        ".RDS")
ps_root <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_phyloseq_root_raw_withnegcontrols_"), 
                        ".RDS")
ps_msoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                         "_phyloseq_msoil_raw_withnegcontrols_"), 
                         ".RDS")
ps_osoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                         "_phyloseq_osoil_raw_withnegcontrols_"), 
                         ".RDS")

# Save the metadata in an object
leaf_meta <- as.data.frame(as.matrix(ps_leaf@sam_data))
root_meta <- as.data.frame(as.matrix(ps_root@sam_data))
msoil_meta <- as.data.frame(as.matrix(ps_msoil@sam_data))
osoil_meta <- as.data.frame(as.matrix(ps_osoil@sam_data))

# Make sequencing depth numeric
leaf_meta$seq_count_dada2 <- as.numeric(leaf_meta$seq_count_dada2)
root_meta$seq_count_dada2 <- as.numeric(root_meta$seq_count_dada2)
msoil_meta$seq_count_dada2 <- as.numeric(msoil_meta$seq_count_dada2)
osoil_meta$seq_count_dada2 <- as.numeric(osoil_meta$seq_count_dada2)

# Bin sequencing depth
leaf_meta <- bin_seq_depth(leaf_meta)
root_meta <- bin_seq_depth(root_meta)
msoil_meta <- bin_seq_depth(msoil_meta)
osoil_meta <- bin_seq_depth(osoil_meta)

### VISUALIZE RAW DATA ########################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
             "/Figures"))

plot_prefilter_seq_depth(leaf_meta, "leaf", 10000, date)
plot_prefilter_seq_depth(root_meta, "root", 10000, date)
plot_prefilter_seq_depth(msoil_meta, "msoil", 8000, date)
plot_prefilter_seq_depth(osoil_meta, "osoil", 8000, date)

### DROP OUTLIERS AND SAMPLES WITH TOO FEW READS ##############################
### LEAF SAMPLES ##############################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
             "/Figures/leaf"))
leaf_raw <- otu_table(ps_leaf)

# see the data structure before removing outliers
id_outliers_evaluate_seq_depth(leaf_raw, leaf_meta, "leaf", yourname, 
                               amplicon, date, "no_outliers_removed") 

# run this step multiple times until you feel you have removed all outliers
# list outliers to remove
if(amplicon == "16S"){
  leaf_outliers <- c("HF069_D3_small_L2_16S_S181", "HF04_QR_Y_edge_L1_16S_S115",
                     "HF06_QR_M_int_L2_16S_S114", "HF06_QR_M_int_L1_16S_S104",
                     "BM03_QR_Y_int_PS_L2_16S_S102", "AB_QR_M_1_L_1_16S_S8", 
                     "HF06_QR_O_edge_L2_16S_S145", "SW08_QR_M_int_L1_16S_S141", 
                     "HF06_QR_M_int_L3_16S_S153", "BB_QR_O_1_L2_16S_S253", 
                     "HF04_QR_M_int_L2_16S_S297") 
} else{
  leaf_outliers <- c("SW08_QR_O_edge_L2_ITS_S66")
}
# filter outliers out of data and metadata
leaf_filter <- leaf_raw[,!colnames(leaf_raw) %in% leaf_outliers]
leaf_meta_filter <- leaf_meta[!rownames(leaf_meta) %in% leaf_outliers,]

# evaluate whether you have removed all outliers; re-run until you feel you have 
# removed all outliers in the above lines
id_outliers_evaluate_seq_depth(leaf_filter, leaf_meta_filter, "leaf", yourname, 
                               amplicon, date, "outliers_removed")

# test different drop thresholds
test_drop_threshold(leaf_filter, leaf_meta_filter, "leaf", yourname, amplicon, 
                    date, 5000)
test_drop_threshold(leaf_filter, leaf_meta_filter, "leaf", yourname, amplicon,
                    date, 10000)

# subset the filtered samples out of the phyloseq
ps_leaf <- subset_samples(ps_leaf, seq_count_dada2 > 10000) # same for ITS and 16S

# write to file
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
write.csv(otu_table(ps_leaf), paste0(yourname, "_", amplicon, 
                                     "_leaf_ASV_table_filteredsamples", date,
                                     ".csv"))
saveRDS(ps_leaf, paste0(yourname, "_", amplicon, 
                        "_leaf_phyloseq_filteredsamples", date, ".RDS"))

### ROOT SAMPLES ##############################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
             "/Figures/root"))
root_raw <- otu_table(ps_root)

# see the data structure before removing outliers
id_outliers_evaluate_seq_depth(root_raw, root_meta, "root", yourname, 
                               amplicon, date, "no_outliers_removed") 

# run this step multiple times until you feel you have removed all outliers
# list outliers to remove
if(amplicon == "16S"){
  root_outliers <- c("AA01_QR_Y_edge_R2_16S_S129", "HF06_QR_M_int_R2_16S_S124", 
                     "JP_QR_Y_pit_R2_16S_S195", "JP_QR_O_pit_R2_16S_S251", 
                     "HF04_QR_O_int_big_R1_16S_S183", 
                     "HF06_QR_O_edge_R2_16S_S138", "HF04_QR_Y_edge_R1_16S_S185", 
                     "HF06_QR_O_edge_R1_16S_S191", "HF06_QR_M_int_R1_16S_S120", 
                     "HF06_QR_Y_int_R1_16S_S196", 
                     "HF04_QR_O_int_big_R2_16S_S120", 
                     "HF069_E3_765_R1_16S_S242")
} else{
  root_outliers <- c("BB_QR_O_1_R2_ITS_S129", "AA01_QR_Y_edge_R1_ITS_S80", 
                     "AA01_QR_Y_edge_R2_ITS_S85")
}
# filter outliers out of data and metadata
root_filter <- root_raw[,!colnames(root_raw) %in% root_outliers]
root_meta_filter <- root_meta[!rownames(root_meta) %in% root_outliers,]

# evaluate whether you have removed all outliers; re-run until you feel you have 
# removed all outliers in the above lines
id_outliers_evaluate_seq_depth(root_filter, root_meta_filter, "root", yourname, 
                               amplicon, date, "outliers_removed")

# test different drop thresholds
test_drop_threshold(root_filter, root_meta_filter, "root", yourname, amplicon, 
                    date, 5000)
test_drop_threshold(root_filter, root_meta_filter, "root", yourname, amplicon,
                    date, 10000)

# subset the filtered samples out of the phyloseq
if(amplicon == "16S"){
  ps_root <- subset_samples(ps_root, seq_count_dada2 > 10000)
} else {
  ps_root <- subset_samples(ps_root, seq_count_dada2 > 5000)
}


# write to file
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
write.csv(otu_table(ps_root), paste0(yourname, "_", amplicon, 
                                     "_root_ASV_table_filteredsamples", date,
                                     ".csv"))
saveRDS(ps_root, paste0(yourname, "_", amplicon, 
                        "_root_phyloseq_filteredsamples", date, ".RDS"))

### MSOIL SAMPLES #############################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
             "/Figures/msoil"))
msoil_raw <- otu_table(ps_msoil)

# see the data structure before removing outliers
id_outliers_evaluate_seq_depth(msoil_raw, msoil_meta, "msoil", yourname, 
                               amplicon, date, "no_outliers_removed") 

# run this step multiple times until you feel you have removed all outliers
# list outliers to remove
if(amplicon == "16S"){
  msoil_outliers <- c() 
} else{
  msoil_outliers <- c("MC_QR_O_grass_M3_ITS_S97", "HF06_QR_Y_int_M2_ITS_S31",
                      "BM03_QR_O_int_M2_ITS_S50", "HF04_QR_Y_edge_M3_ITS_S41",
                      "BM03_QR_Y_int_PS_M2_ITS_S15")
}
# filter outliers out of data and metadata
msoil_filter <- msoil_raw[,!colnames(msoil_raw) %in% msoil_outliers]
msoil_meta_filter <- msoil_meta[!rownames(msoil_meta) %in% msoil_outliers,]

# evaluate whether you have removed all outliers; re-run until you feel you have 
# removed all outliers in the above lines
id_outliers_evaluate_seq_depth(msoil_filter, msoil_meta_filter, "root", 
                               yourname, amplicon, date, "outliers_removed")

# test different drop thresholds
test_drop_threshold(msoil_filter, msoil_meta_filter, "msoil", yourname, 
                    amplicon, date, 5000)
test_drop_threshold(msoil_filter, msoil_meta_filter, "msoil", yourname, 
                    amplicon, date, 8000)

# subset the filtered samples out of the phyloseq
ps_msoil <- subset_samples(ps_leaf, seq_count_dada2 > 8000) # same for 16S and ITS

# write to file
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
write.csv(otu_table(ps_msoil), paste0(yourname, "_", amplicon, 
                                      "_msoil_ASV_table_filteredsamples", date,
                                      ".csv"))
saveRDS(ps_msoil, paste0(yourname, "_", amplicon,
                         "_msoil_phyloseq_filteredsamples", date, ".RDS"))

### OSOIL SAMPLES #############################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon, 
             "/Figures/osoil"))
osoil_raw <- otu_table(ps_osoil)

# see the data structure before removing outliers
id_outliers_evaluate_seq_depth(osoil_raw, osoil_meta, "osoil", yourname, 
                               amplicon, date, "no_outliers_removed") 

# run this step multiple times until you feel you have removed all outliers
# list outliers to remove
if(amplicon == "16S"){
  osoil_outliers <- c("SW08_QR_O_edge_O1_16S_S128", "SW08_QR_M_int_O1_16S_S135",
                      "HW07_QR_Y_edge_O1_16S_S192", "SB_QR_O_pit_O2_16S_S39")
} else{
  osoil_outliers <- c()
}
# filter outliers out of data and metadata
osoil_filter <- osoil_raw[,!colnames(osoil_raw) %in% osoil_outliers]
osoil_meta_filter <- osoil_meta[!rownames(osoil_meta) %in% osoil_outliers,]

# evaluate whether you have removed all outliers; re-run until you feel you have 
# removed all outliers in the above lines
id_outliers_evaluate_seq_depth(osoil_filter, osoil_meta_filter, "root", 
                               yourname, amplicon, date, "outliers_removed")

# test different drop thresholds
test_drop_threshold(osoil_filter, osoil_meta_filter, "osoil", yourname, 
                    amplicon, date, 5000)
test_drop_threshold(osoil_filter, osoil_meta_filter, "osoil", yourname, 
                    amplicon, date, 8000)

# subset the filtered samples out of the phyloseq
if(amplicon == "16S"){
  ps_osoil <- subset_samples(ps_osoil, seq_count_dada2 > 8000)
} else{
  ps_osoil <- subset_samples(ps_osoil, seq_count_dada2 > 5000)
}

# write to file
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
write.csv(otu_table(ps_osoil), paste0(yourname, "_", amplicon, 
                                      "_osoil_ASV_table_filteredsamples", date,
                                      ".csv"))
saveRDS(ps_osoil, paste0(yourname, "_", amplicon, 
                         "_osoil_phyloseq_filteredsamples", date, ".RDS"))

### WRITE TO METADATA WHETHER SAMPLES WERE DROPPED OR NOT #####################
if(edit_metadata == "Y"){
  # read in metadata
  setwd(paste0(pwd,"01_Collect_Data/01_Sample_Metadata"))
  metadata <- read_in_file(getwd(), paste0("atherton_sample_metadata_", 
                                           amplicon, "_"), ".csv")
  # initialize column
  metadata$sequences_dropped <- "No"
  
  # record NA for samples without sequence data
  metadata$sequences_dropped[is.na(metadata$sample_name)] <- NA
  
  # record outlier for outlier samples
  metadata$sequences_dropped[metadata$sample_name %in% leaf_outliers] <- "Outlier"
  metadata$sequences_dropped[metadata$sample_name %in% root_outliers] <- "Outlier"
  metadata$sequences_dropped[metadata$sample_name %in% msoil_outliers] <- "Outlier"
  metadata$sequences_dropped[metadata$sample_name %in% osoil_outliers] <- "Outlier"
  
  # record filtered out for samples with low sequencing depth
  metadata$sequences_dropped[metadata$sample_type == "Leaf" & 
                               metadata$seq_count_dada2 < 10000] <- "Filtered out"
  if(amplicon == "16S"){
    metadata$sequences_dropped[metadata$sample_type == "Root" & 
                                 metadata$seq_count_dada2 < 10000] <- "Filtered out"
  } else{
    metadata$sequences_dropped[metadata$sample_type == "Root" & 
                                 metadata$seq_count_dada2 < 5000] <- "Filtered out"
  }
  metadata$sequences_dropped[metadata$sample_type == "Soil" & 
                               metadata$soil_horizon == "M" & 
                               metadata$seq_count_dada2 < 8000] <- "Filtered out"
  if(amplicon == "16S"){
    metadata$sequences_dropped[metadata$sample_type == "Soil" & 
                                 metadata$soil_horizon == "O" & 
                                 metadata$seq_count_dada2 < 8000] <- "Filtered out"
  } else{
    metadata$sequences_dropped[metadata$sample_type == "Soil" & 
                                 metadata$soil_horizon == "O" & 
                                 metadata$seq_count_dada2 < 5000] <- "Filtered out"
  }
  
  write.csv(metadata, paste0(getwd(), yourname, "_sample_metadata_", 
                             amplicon, date, ".csv"), row.names = FALSE)
}
