read_in_file <- function(path, prefix, extension){
  files <- list.files(path, pattern=prefix)
  files <- files[endsWith(files, extension)]
  if(length(files) > 1){
    dates <- gsub(prefix, "", files)
    dates <- gsub("\\..*", "", dates)
    dates <- as.numeric(dates)
    dates <- dates[!is.na(dates)]
    dates <- dates[which(dates == max(dates))]
    if(endsWith(path, "/")){
      file <- paste0(path,prefix,dates,extension) # keep the most recent version of the file to use
    } else{
      file <- paste0(path,"/",prefix,dates,extension) # keep the most recent version of the file to use
    }
  } else{
    if(endsWith(path, "/")){
      file <- paste0(path,files)
    } else{
      file <- paste0(path,"/",files)
    }
  }
  if(extension == ".RDS"){
    data <- readRDS(file)
  } else{
    data <- vroom(file)    
  }
  return(data)
}

read_in_dada2_asv_table <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  colnames(data)[1] <- "ASV" # rename the first column for merging purposes
  return(data)
}

read_in_dada2_taxonomy <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  data <- data[,c(1,2)] # keep only the ASV ID and taxonomy information
  return(data)
}

subset_phyloseq <- function(ps_otu, ps_meta, type, nc_list){
  metadata$seq_bin <- NA
  if(type == "Leaf"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Leaf") 
  } else if(type == "Root"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Root") 
  } else if(type == "MSoil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "M"))
  } else if(type == "OSoil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "O"))
  }
  # add negative controls
  nc <- colnames(ps_otu)[grep(paste(nc_list, collapse = "|"), colnames(ps_otu))]
  ps_nc <- prune_samples(nc, ps_meta)
  ps_subset_w_nc <- merge_phyloseq(ps_subset, ps_nc)
  return(ps_subset_w_nc)
}

bin_seq_depth <- function(metadata){
  for(i in 1:nrow(metadata)){
    depth <- metadata$seq_count_dada2[i]
    if(depth < 10000){
      metadata$seq_bin[i] <- "0-10k"
    } else if(depth < 20000){
      metadata$seq_bin[i] <- "10-20k"
    } else if(depth < 30000){
      metadata$seq_bin[i] <- "20-30k"
    } else if(depth < 40000){
      metadata$seq_bin[i] <- "30-40k"
    } else if(depth < 50000){
      metadata$seq_bin[i] <- "40-50k"
    } else if(depth < 60000){
      metadata$seq_bin[i] <- "50-60k"
    } else if(depth < 70000){
      metadata$seq_bin[i] <- "60-70k"
    } else if(depth < 80000){
      metadata$seq_bin[i] <- "70-80k"
    } else if(depth < 90000){
      metadata$seq_bin[i] <- "80-90k"
    } else if(depth < 100000){
      metadata$seq_bin[i] <- "90-100k"
    } else{
      metadata$seq_bin[i] <- "100k"
    }
  }
  return(metadata)
}

plot_prefilter_seq_depth <- function(metadata, sample_type, proposed_threshold,
                                     date){
  ggplot(metadata, 
         aes(as.numeric(seq_count_dada2), 
             fill = sequencing_batch)) + 
    geom_histogram(binwidth = 1000) + 
    theme_bw() + 
    ggtitle(paste0(sample_type, " sequencing depth histogram")) + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_dada2))), 
               color="blue", 
               linetype="dashed") + 
    geom_vline(aes(xintercept=proposed_threshold), 
               color="red", 
               linetype="dashed") +
    labs(x = "Read count after DADA2", y = "count", fill = "Sequencing batch")
  
  ggsave(paste0(sample_type, "/", yourname, "_", amplicon, "_", sample_type, 
                "_sequencing_depth_prefilter", date,".png"), width = 9, 
         height = 5, units = "in", dpi = 300)
}

id_outliers_evaluate_seq_depth <- function(data, metadata, sample_type, yourname, 
                                    amplicon, date, rep){
  # make a transposed verison of the leaf sequencing data
  data_t <- as.data.frame(t(data))

  # calculate Aitchison distance for leaf samples
  aitch_data <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch_data <- as.matrix(aitch_data)

  # test for sequencing batch and depth effect
  print(paste0("Testing for sequencing batch effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.factor(metadata$sequencing_batch)))
  print(paste0("Testing for sequencing depth effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.numeric(metadata$seq_count_dada2)))

  all.MDS <- metaMDS(aitch_data,k=2,zerodist="add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)

  # plot samples by different variables
  ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2)) +
    geom_text(label = rownames(coordinates),
              size = 2) +
    stat_ellipse() +
    theme_bw()
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_identify_outliers_", rep, date, ".png"), width = 7, 
         height = 5, units = "in", dpi = 300)

  by_batch <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(sequencing_batch))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_seqdepth <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(seq_bin))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_treeage <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_age))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_treepittype <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_pit_type))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_species <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_species))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_site <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(site_name))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()
  
  multipanel <- grid.arrange(by_batch, by_seqdepth, by_treeage, by_treepittype, 
                             by_species, by_site, nrow = 2, ncol = 3)
  
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_predrop_datastructure_", rep, date, ".png"), multipanel, 
         width = 21, height = 10, units = "in", dpi = 300)
}

test_drop_threshold <- function(data, metadata, sample_type, yourname, 
                                 amplicon, date, threshold){
  # drop samples < threshold
  metadata_drop <- metadata[metadata$seq_count_dada2 > threshold,]
  data_drop <- data[,colnames(data) %in% rownames(metadata_drop)]
  
  # make a transposed verison of the sequencing data
  data_t <- as.data.frame(t(data_drop))
  
  # calculate Aitchison distance for samples
  aitch <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print("Testing sequencing batch effect:")
  print(adonis2(aitch ~ as.factor(metadata_drop$sequencing_batch)))
  
  print("Testing sequencing depth effect:")
  adonis2(aitch ~ as.numeric(metadata_drop$seq_count_dada2))
  
  # calculate NMDS
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates<-data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata_drop)
  
  # plot samples by different variables
  by_batch <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                      col = as.factor(sequencing_batch))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_seqdepth <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                         col = as.factor(seq_bin))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_treeage <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                        col = as.factor(tree_age))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_treepittype <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                            col = as.factor(tree_pit_type))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_species <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                        col = as.factor(tree_species))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_site <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                     col = as.factor(site_name))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  multipanel <- grid.arrange(by_batch, by_seqdepth, by_treeage, by_treepittype, 
                             by_species, by_site, nrow = 2, ncol = 3)
  
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_drop", threshold,"_", date, ".png"), multipanel, 
         width = 21, height = 10, units = "in", dpi = 300)
  
}
