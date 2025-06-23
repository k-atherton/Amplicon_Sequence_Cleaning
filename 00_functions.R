ensure_directory_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    message("Directory created: ", folder_path)
  } else {
    message("Directory already exists: ", folder_path)
  }
}

process_file <- function(file_path) {
  print(paste0("Reading in file: ", file_path))
  if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    return(read.table(file_path, sep = "\t", header = TRUE, 
                      stringsAsFactors = FALSE, comment.char = ""))
  } else {
    warning(paste("Unsupported file type:", file_path))
    return(NULL)
  }
}

# Read the list of file paths from a .txt file
read_file_paths <- function(txt_file) {
  print("Reading in file names...")
  file_paths <- readLines(txt_file, warn = FALSE)
  file_paths <- trimws(file_paths)  # Remove any leading/trailing whitespace
  return(file_paths)
}

# Main function to load all files
load_files <- function(txt_file) {
  file_paths <- read_file_paths(txt_file)
  files_list <- lapply(file_paths, process_file)
  names(files_list) <- basename(file_paths)  # Name elements by file name
  return(files_list)
}

read_in_file <- function(path, prefix, extension){
  print(paste0("Looking for ", prefix," file ..."))
  files <- list.files(path, pattern=prefix)
  files <- files[endsWith(files, extension)]
  if(length(files) > 1){
    print("More than one file with that name, looking for most recent version.")
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
  print(paste0("Reading in ", file))
  if(extension == ".RDS"){
    data <- readRDS(file)
  } else{
    data <- vroom(file)    
  }
  return(data)
}

read_in_dada2_asv_table <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  print("Renaming the first column 'ASV' for merging purposes.")
  colnames(data)[1] <- "ASV" # rename the first column for merging purposes
  return(data)
}

read_in_dada2_taxonomy <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  print("Keeping only the ASV ID and taxonomy information columns.")
  data <- data[,c(1,2)] # keep only the ASV ID and taxonomy information
  return(data)
}

add_16s_guild <- function(tax, functional_guilds){
  print("Formatting the taxonomic rank text.")
  # format the taxonomic ranks
  tax$Kingdom <- gsub("d__", "", tax$Kingdom)
  tax$Phylum <- gsub(" p__", "", tax$Phylum)
  tax$Class <- gsub(" c__", "", tax$Class)
  tax$Order <- gsub(" o__", "", tax$Order)
  tax$Family <- gsub(" f__", "", tax$Family)
  tax$Genus <- gsub(" g__", "", tax$Genus)
  tax$Species <- gsub(" s__", "", tax$Species)
  tax$Species <- gsub("_", " ", tax$Species)
  
  print("Separating the functional groups by taxonomic level.")
  # separate the functional groups by the taxonomic level
  species_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Species"),]
  genus_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Genus"),]
  family_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Family"),]
  order_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Order"),]
  class_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Class"),]
  phylum_groups <- functional_guilds[which(functional_guilds$`Taxonomic level` == "Phylum"),]
  
  print("Keeping only the taxonomic level and guild information.")
  # keep only the taxonomic level and the guild information
  species_groups <- species_groups[,c(2,4)]
  genus_groups <- genus_groups[,c(2,4)]
  family_groups <- family_groups[,c(2,4)]
  order_groups <- order_groups[,c(2,4)]
  class_groups <- class_groups[,c(2,4)]
  phylum_groups <- phylum_groups[,c(2,4)]
  
  # keep only distinct information from the database and collapse any duplicates
  print("Ensuring distinct species-level guild information.")
  classification <- c()
  for(i in 1:length(unique(species_groups$Taxon))){
    taxa <- unique(species_groups$Taxon)[i]
    classes <- species_groups$Classification[which(species_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = "_"))
  }
  species_groups <- data.frame(unique(species_groups$Taxon), classification)
  colnames(species_groups) <- c("Taxon", "Classification")
  
  print("Ensuring distinct genus-level guild information.")
  genus_groups <- distinct(genus_groups)
  classification <- c()
  for(i in 1:length(unique(genus_groups$Taxon))){
    taxa <- unique(genus_groups$Taxon)[i]
    classes <- genus_groups$Classification[which(genus_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = ", "))
  }
  genus_groups <- data.frame(unique(genus_groups$Taxon), classification)
  colnames(genus_groups) <- c("Taxon", "Classification")
  
  print("Ensuring distinct family-level guild information.")
  family_groups <- distinct(family_groups)
  classification <- c()
  for(i in 1:length(unique(family_groups$Taxon))){
    taxa <- unique(family_groups$Taxon)[i]
    classes <- family_groups$Classification[which(family_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = ", "))
  }
  family_groups <- data.frame(unique(family_groups$Taxon), classification)
  colnames(family_groups) <- c("Taxon", "Classification")
  
  print("Ensuring distinct order-level guild information.")
  order_groups <- distinct(order_groups)
  classification <- c()
  for(i in 1:length(unique(order_groups$Taxon))){
    taxa <- unique(order_groups$Taxon)[i]
    classes <- order_groups$Classification[which(order_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = ", "))
  }
  order_groups <- data.frame(unique(order_groups$Taxon), classification)
  colnames(order_groups) <- c("Taxon", "Classification")
  
  print("Ensuring distinct class-level guild information.")
  class_groups <- distinct(class_groups)
  classification <- c()
  for(i in 1:length(unique(class_groups$Taxon))){
    taxa <- unique(class_groups$Taxon)[i]
    classes <- class_groups$Classification[which(class_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = ", "))
  }
  class_groups <- data.frame(unique(class_groups$Taxon), classification)
  colnames(class_groups) <- c("Taxon", "Classification")
  
  print("Ensuring distinct phylum-level guild information.")
  phylum_groups <- distinct(phylum_groups)
  classification <- c()
  for(i in 1:length(unique(phylum_groups$Taxon))){
    taxa <- unique(phylum_groups$Taxon)[i]
    classes <- phylum_groups$Classification[which(phylum_groups$Taxon == taxa)]
    classification <- c(classification, paste(classes, collapse = ", "))
  }
  phylum_groups <- data.frame(unique(phylum_groups$Taxon), classification)
  colnames(phylum_groups) <- c("Taxon", "Classification")
  
  tax$guild <- NA
  
  print("Matching guild information to the taxonomy table")
  for(i in 1:nrow(tax)){
    s <- tax$Species[i]
    g <- tax$Genus[i]
    f <- tax$Family[i]
    o <- tax$Order[i]
    c <- tax$Class[i]
    p <- tax$Phylum[i]
    k <- tax$Kingdom[i]
    if(s %in% species_groups$Taxon){
      tax$guild[i] <- species_groups$Classification[which(species_groups$Taxon == s)]
    } else if(g %in% genus_groups$Taxon){
      tax$guild[i] <- genus_groups$Classification[which(genus_groups$Taxon == g)]
    } else if(f %in% family_groups$Taxon){
      tax$guild[i] <- family_groups$Classification[which(family_groups$Taxon == f)]
    } else if(o %in% order_groups$Taxon){
      tax$guild[i] <- order_groups$Classification[which(order_groups$Taxon == o)]
    } else if(c %in% class_groups$Taxon){
      tax$guild[i] <- class_groups$Classification[which(class_groups$Taxon == c)]
    } else if(p %in% phylum_groups$Taxon){
      tax$guild[i] <- phylum_groups$Classification[which(phylum_groups$Taxon == p)]
    } else if(k == "Archaea"){
      tax$guild[i] <- "Archaea"
    }
  }
  return(tax)
}

add_its_guild <- function(tax, functional_guilds){
  # format the taxonomic ranks
  print("Formatting the taxonomic rank text.")
  tax$Kingdom <- gsub("k__", "", tax$Kingdom)
  tax$Phylum <- gsub("p__", "", tax$Phylum)
  tax$Class <- gsub("c__", "", tax$Class)
  tax$Order <- gsub("o__", "", tax$Order)
  tax$Family <- gsub("f__", "", tax$Family)
  tax$Genus <- gsub("g__", "", tax$Genus)
  tax$Species <- gsub("s__", "", tax$Species)
  tax$Species <- gsub("_", " ", tax$Species)
  
  # keep the necessary columns
  print("Keeping only the genus, primary lifestyle, plant pathogenic capacity, and animal biotrophic capacity columns.")
  functional_guilds <- functional_guilds[,which(colnames(functional_guilds) %in% 
                                                  c("GENUS", "primary_lifestyle",
                                                    "Plant_pathogenic_capacity_template",
                                                    "Animal_biotrophic_capacity_template"))]
  colnames(functional_guilds)[1] <- "Genus"
  
  # keep only distinct rows
  print("Keeping only distinct rows in the database.")
  functional_guilds <- distinct(functional_guilds)
  
  # keep the rownames (feature ids)
  print("Keeping the feature IDs (rownames) in a column for merging.")
  tax$Feature_ID <- row.names(tax)
  
  # merge functional guilds and tax information
  print("Merging functional guilds and taxonomic information.")
  tax <- merge(functional_guilds, tax, by = "Genus", all.y = TRUE)
  
  # keep feature ids in rownames
  print("Moving feature IDs back to rownames.")
  row.names(tax) <- tax$Feature_ID
  tax <- tax[,-which(colnames(tax) %in% c("Feature_ID"))]
  return(tax)
}

subset_phyloseq <- function(ps_otu, ps_meta, type, nc_list){
  print(paste0("Subsetting ", type, " samples from the whole dataset."))
  
  # Extract the metadata as a dataframe
  metadata <- ps_meta@sam_data
  metadata <- metadata[which(metadata$sample_type == type),]
  
  # Subset samples
  ps_subset <- prune_samples(rownames(metadata), ps_meta)
  
  # add negative controls
  print("Adding negative controls to the subsetted dataset.")
  nc <- colnames(ps_otu)[grep(paste(nc_list, collapse = "|"), colnames(ps_otu))]
  ps_nc <- prune_samples(nc, ps_meta)
  ps_subset_w_nc <- merge_phyloseq(ps_subset, ps_nc)
  return(ps_subset_w_nc)
}

bin_seq_depth <- function(metadata){
  metadata$seq_bin <- NA
  print("Binning the sequencing depths for visualization purposes.")
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
  print("Plotting the pre-filtered sequencing depth as a histogram.")
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
                "_sequencing_depth_prefilter_drop", proposed_threshold, date,
                ".png"), width = 9, height = 5, units = "in", dpi = 300)
}

id_outliers_evaluate_seq_depth <- function(data, metadata, sample_type, 
                                           yourname, amplicon, date, rep, 
                                           color_vars){
  # Transpose the sequencing data
  print("Transposing the sequencing data:")
  data_t <- as.data.frame(t(data))
  
  # Calculate Aitchison distance
  print("Calculating the Aitchison's distance for the samples")
  aitch_data <- aDist(data_t + 1)  # Add 1 to avoid zero issues
  aitch_data <- as.matrix(aitch_data)
  
  # Test for sequencing batch and depth effects
  print(paste0("Testing for sequencing batch effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data ~ as.factor(metadata$sequencing_batch)))
  print(paste0("Testing for sequencing depth effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data ~ as.numeric(metadata$seq_count_dada2)))
  
  # Perform NMDS
  print("Calculating NMDS for plotting.")
  all.MDS <- metaMDS(aitch_data, k = 2, zerodist = "add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)
  
  # Generate NMDS plots for each variable in color_vars
  print("Plotting samples by different variables")
  plot_list <- lapply(color_vars, function(var) {
    ggplot(coordinates, aes_string(x = "NMDS1", y = "NMDS2", col = var)) +
      geom_point() +
      stat_ellipse() +
      theme_bw() +
      ggtitle(paste("NMDS colored by", var))
  })
  
  # Arrange plots in a grid
  print("Arrange plots in a grid")
  multipanel <- do.call(grid.arrange, c(plot_list, 
                                        nrow = ceiling(length(plot_list) / 3), 
                                        ncol = min(3, length(plot_list))))
  
  # Save the multipanel figure
  print("Saving multipanel figures")
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_predrop_datastructure_", rep, date, ".png"), 
         multipanel, width = (7*ceiling(length(plot_list)/3)), 
         height = (5*min(3, length(plot_list))), units = "in", dpi = 300)
  
  coordinates$sample_name <- rownames(coordinates)
  
  # Make plot with sample name as point
  print("Make outlier ID figure")
  range_x <- range(coordinates$NMDS1, na.rm = TRUE)
  outlier_id <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                        label = rownames(coordinates))) +
    geom_text(size = 2) +
    stat_ellipse() +
    theme_bw() + xlim(c(range_x[1] - 25, range_x[2] + 25)) +
    ggtitle(paste0("Identify outliers in ", sample_type, " samples"))
  
  # Save the figure
  print("Save figure")
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_outlier_id_plot_", rep, date, ".png"), 
         outlier_id, width = (7*ceiling(length(plot_list)/3)), 
         height = (5*min(3, length(plot_list))), units = "in", dpi = 300)
}

test_drop_threshold <- function(data, metadata, sample_type, yourname, 
                                amplicon, date, threshold, color_vars){
  # drop samples < threshold
  print(paste0("Keeping samples with dada2 read count > ", threshold, "."))
  metadata_drop <- metadata[as.numeric(metadata$seq_count_dada2) > threshold,]
  data_drop <- data[,colnames(data) %in% rownames(metadata_drop)]
  print(nrow(metadata_drop))
  
  # make a transposed verison of the sequencing data
  print("Transposing the sequencing data.")
  data_t <- as.data.frame(t(data_drop))
  
  # calculate Aitchison distance for samples
  print("Calculating Aitchison's distance for samples.")
  aitch <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print("Testing sequencing batch effect:")
  print(unique(metadata_drop$sequencing_batch))
  print(adonis2(aitch ~ as.factor(metadata_drop$sequencing_batch)))
  
  print("Testing sequencing depth effect:")
  adonis2(aitch ~ as.numeric(metadata_drop$seq_count_dada2))
  
  # calculate NMDS
  print("Calculating NMDS for plotting purposes.")
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates<-data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata_drop)
  
  print("Plotting samples by different variables")
  plot_list <- lapply(color_vars, function(var) {
    ggplot(coordinates, aes_string(x = "NMDS1", y = "NMDS2", col = var)) +
      geom_point() +
      stat_ellipse() +
      theme_bw() +
      ggtitle(paste("NMDS colored by", var))
  })
  
  # Arrange plots in a grid
  multipanel <- do.call(grid.arrange, c(plot_list, 
                                        nrow = ceiling(length(plot_list) / 3), 
                                        ncol = min(3, length(plot_list))))
  
  # Save the multipanel figure
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, "_NMDS_drop_", 
                threshold, "_datastructure", date, ".png"), 
         multipanel, width = (7*ceiling(length(plot_list)/3)), 
         height = (5*min(3, length(plot_list))), units = "in", dpi = 300)
  
}

plot_filter_seq_depth <- function(metadata, sample_type, threshold, date){
  print("Plotting the pre-filtered sequencing depth as a histogram.")
  ggplot(metadata, 
         aes(as.numeric(seq_count_dada2), 
             fill = sequencing_batch)) + 
    geom_histogram(binwidth = 1000) + 
    theme_bw() + 
    ggtitle(paste0(sample_type, " sequencing depth histogram")) + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_dada2))), 
               color="blue", 
               linetype="dashed") + 
    geom_vline(aes(xintercept=threshold), 
               color="red", 
               linetype="dashed") +
    labs(x = "Read count after DADA2", y = "count", fill = "Sequencing batch")
  
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_sequencing_depth_filtered_drop", threshold, date,
                ".png"), width = 9, height = 5, units = "in", dpi = 300)
}

decontaminate_samples <- function(ps, sample_type, yourname, amplicon, date){
  print(sample_type)
  # get the batch IDs
  print("Getting the batch IDs.")
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  batches <- unique(metadata$sequencing_batch[which(metadata$sample_type != "Negative Control")])
  
  pre_decontam_pngs <- list()
  contaminants_pngs <- list()
  
  j <- 1
  for(i in 1:length(batches)){
    batch <- batches[i]
    print(batch)
    # subset the batch
    print("Subsetting the batch.")
    batch_data <- prune_samples(sample_data(ps)$sequencing_batch == batch, ps)
    meta_batch <- as.data.frame(as.matrix(batch_data@sam_data))
    
    # order for rank
    print("Ordering the data for plotting by rank.")
    meta_batch <- meta_batch[order(meta_batch$seq_count_dada2),]
    meta_batch$index <- seq(nrow(meta_batch))
    
    # plot the pre-decontam data
    print("Plotting the pre-decontam data.")
    pre_decontam_pngs[[i]] <- ggplot(meta_batch, aes(x = index, 
                                                     y = as.numeric(seq_count_dada2), 
                                                     color = is_control)) +
      geom_point() + theme_bw() + 
      ggtitle(paste0(sample_type, " ", batch))
    
    # run decontam
    print("Running decontam.")
    if(length(unique(meta_batch$is_control)) > 1){
      decontam_batch <- isContaminant(batch_data, neg = "is_control", 
                                      method = "prevalence")
      
      # visualize contaminants
      print("Visualizing contaminants.")
      pa <- transform_sample_counts(batch_data, function(abund) 1*(abund>0))
      neg <- prune_samples(sample_data(batch_data)$is_control == TRUE, 
                           batch_data)
      pos <- prune_samples(sample_data(batch_data)$is_control == FALSE, 
                           batch_data)
      df_batch <- data.frame(pos=taxa_sums(pos), neg=taxa_sums(neg), 
                             contaminant=decontam_batch$contaminant)
      contaminants_pngs[[j]] <- ggplot(df_batch, aes(x= neg, y = pos, 
                                                     color = contaminant)) +
        geom_point() + xlab("Prevalence (Negative Controls)") + 
        ylab("Prevalence (True Samples)") + theme_bw() +
        ggtitle(paste0(sample_type, " ", batch))
      
      j <- j + 1
    } else{
      decontam_batch <- isContaminant(batch_data, conc = "dna_conc", 
                                      method = "frequency")
      
      # visualize contaminants
      print("Visualizing contaminants.")
      freq_contaminants <- plot_frequency(batch_data, 
                                          taxa_names(batch_data)[sample(which(decontam_batch$contaminant),
                                                                        sum(decontam_batch$contaminant, 
                                                                            na.rm = TRUE))], 
                                          conc="dna_conc") + 
        xlab("DNA Concentration")
      
      png(filename=paste0("Figures/", sample_type, "/", 
                          yourname, "_", amplicon, "_", sample_type, "_",
                          batch, "_decontam_frequencymethod", date, 
                          ".png"))
      plot(freq_contaminants)
      dev.off()
    }
    # identify contaminants
    print("Identifying contaminant taxonomy.")
    taxonomy <- as.data.frame(as.matrix(batch_data@tax_table))
    contaminants <- taxonomy[which(rownames(taxonomy) %in% 
                                     rownames(decontam_batch)[which(decontam_batch$contaminant == TRUE)]),]
    
    # write information to table
    print("Writing contaminant information to table.")
    write.csv(contaminants, paste0("Contaminant_Taxonomy/", sample_type, "/", 
                                   yourname, "_", amplicon, "_", sample_type, 
                                   "_", batch, "_contaminants", date, 
                                   ".csv"))
    
    # remove contaminants
    print("Pruning contaminats.")
    decontam_data <- prune_taxa(!decontam_batch$contaminant, batch_data)
    
    # merge back into one dataset
    print("Merging back into one dataset.")
    if(i > 1){
      final_data <- merge_phyloseq(final_data, decontam_data)
    } else{
      final_data <- decontam_data
    }
  }
  # save multipanel figures
  print("Saving visualizations.")
  multipanel_pre_decontam <- do.call(grid.arrange, pre_decontam_pngs)
  multipanel_contaminants <- do.call(grid.arrange, contaminants_pngs)
  
  ggsave(paste0("Figures/", sample_type, "/",
                yourname, "_", amplicon, "_", sample_type, 
                "_predecontam", date, ".png"), multipanel_pre_decontam, 
         width = 21, height = 10, units = "in", dpi = 300)
  
  ggsave(paste0("Figures/", sample_type, "/",
                yourname, "_", amplicon, "_", sample_type, 
                "_decontam_prevalencemethod", date, ".png"), 
         multipanel_contaminants, width = 21, height = 10, units = "in", 
         dpi = 300)
  
  # remove negative controls
  print("Removing negative controls.")
  final_data <- subset_samples(final_data, is_control == "FALSE")
  
  if(amplicon == "16S"){
    # remove chloroplasts, mitochondria, and anything not a bacteria
    print("Removing non-bacterial/archeal reads, chloroplasts, and mitochondria.")
    final_data <- subset_taxa(final_data, Kingdom != "Unassigned")
    final_data <- subset_taxa(final_data, Kingdom != "Eukaryota")
    final_data <- subset_taxa(final_data, Family != "Mitochondria")
    final_data <- subset_taxa(final_data, Order != "Chloroplast")
  } else{
    # remove anything not a fungi
    print("Removing non-fungal reads.")
    final_data <- subset_taxa(final_data, Kingdom == "Fungi")
  }
  
  # write ASV table to CSV
  print("Writing decontaminated ASV table.")
  asv_table <- as.data.frame(as.matrix(final_data@otu_table))
  write.csv(asv_table, paste0(yourname, "_", amplicon, "_", sample_type, 
                              "_ASV_table_decontam", date, ".csv"))
  return(final_data)
}

plot_decontam_seq_depth <- function(metadata, sample_type, date){
  print(sample_type)
  print("Plotting post-decontam sequencing depth.")
  ggplot(metadata, 
         aes(as.numeric(seq_count_decontam), 
             fill = sequencing_batch)) + 
    geom_histogram(binwidth = 1000) + 
    theme_bw() + 
    ggtitle(paste0(sample_type, " sequencing depth histogram")) + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_dada2))), 
               color="blue", 
               linetype="dashed") + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_decontam), 
                                     na.rm = TRUE)), 
               color="red", 
               linetype="dashed") +
    labs(x = "Read count after decontam", y = "count", 
         fill = "Sequencing batch")
  
  ggsave(paste0(sample_type, "/", yourname, "_", amplicon, "_", sample_type, 
                "_sequencing_depth_postdecontam", date,".png"), width = 9, 
         height = 5, units = "in", dpi = 300)
}

check_batch_effect <- function(ps, amplicon, yourname, date, variables){
  # extract asv table
  print("Extract ASV table from phyloseq object.")
  sequence_data <- as.data.frame(as.matrix(ps@otu_table))
  
  # extract metadata
  print("Extract metadata from phyloseq object.")
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  
  # make a transposed version of asv table
  print("Make transposed version of ASV table.")
  sequence_data_t <- as.data.frame(t(sequence_data))
  
  # calculate Aitchison distance for samples
  print("Calculate Aitchison distance for samples.")
  aitch <- aDist(sequence_data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  lapply(variables, function(var) {
    l <- length(unique(metadata[,which(colnames(metadata) == var)]))
    if(l > 1){
      print(paste0("Calculating effect of ", var))
      print(adonis2(aitch ~ as.factor(metadata[,which(colnames(metadata) == var)])))
    }
  })
  
  print("Calculating NMDS for visualization purposes.")
  all.MDS <- metaMDS(aitch, k=2, zerodist="add", trymax = 100)
  
  coordinates <- as.data.frame(all.MDS$points)
  coordinates <- cbind(coordinates, metadata)
  
  # plot samples by different variables
  print("Plotting samples by different variables.")
  plot_list <- lapply(variables, function(var) {
    ggplot(coordinates, aes_string(x = "MDS1", y = "MDS2", col = var)) +
      geom_point() +
      stat_ellipse() +
      theme_bw() +
      ggtitle(paste("NMDS colored by", var))
  })
  
  multipanel <- do.call(grid.arrange, c(plot_list, 
                                        nrow = ceiling(length(plot_list) / 3), 
                                        ncol = min(3, length(plot_list))))
  
  ggsave(paste0(yourname, "_", amplicon, "_check_batch_effect", date, ".png"), 
         multipanel, width = 28, height = 10, units = "in", dpi = 300)
}

batch_correct <- function(ps, amplicon, yourname, date, batch, covars){
  # get asv table from phyloseq object
  print("Extracting ASV table from phyloseq object.")
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # remove any taxa with 0 counts and any samples under 1000 counts as 
  # ComBat_seq will not work if you have samples with 10k+ counts AND really 
  # low counts in the same dataset
  print("Removing taxa with 0 counts.")
  asv_table <- asv_table[which(rowSums(asv_table) > 0),]
  asv_table <- asv_table[,order(colnames(asv_table))]
  asv_table <- as.matrix(asv_table)
  
  # get metadata from phyloseq object
  print("Getting metadata from phyloseq object.")
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  
  # keep metadata samples which are still in the asv table
  print("Keeping metadata samples that are still in ASV table.")
  metadata <- metadata[which(rownames(metadata) %in% colnames(asv_table)),]
    print("Removing leaf samples from batch control -- too few reads.")
    metadata <- metadata[which(metadata$sample_type != "Leaf"),]
    asv_table <- asv_table[,which(colnames(asv_table) %in% rownames(metadata))]
    asv_table <- asv_table[,order(colnames(asv_table))]
  
  # order the metadata by sample name
  print("Ordering the metadata by sample name (same as the ASV table)")
  metadata <- metadata[order(row.names(metadata)),]
  
  # batch is defined as the sample type, as the samples were extracted and 
  # amplified with different methods, but the sequencing batch did not have
  # a batch effect
  print(paste0("Defining batch as ", batch))
  
  metadata$sample_type <- paste0(metadata$soil_horizon, metadata$sample_type)
  metadata$sample_type <- gsub(".*Root", "Root", metadata$sample_type)
  metadata$sample_type <- gsub("NALeaf", "Leaf", metadata$sample_type)
  
  batch_variable <- metadata[,which(colnames(metadata) == batch)]
  print(unique(batch_variable))
  
  # the covariates we want to maintain the signature of are tree age and tree 
  # pit type
  if(length(covars > 1)){
    # batch correct with ComBat_seq
    print(paste0("Defining covariates that should maintain data signature after batch control: ", covars))
    covar_mod <- as.data.frame(metadata[,which(colnames(metadata) %in% covars)])
    if(storage.mode(asv_table) == "double"){
      asv_table <- round(asv_table)
      storage.mode(asv_table) <- "integer"
    }
    batch_variable <- as.factor(batch_variable)
    print("Running batch control.")
    batch_corrected <- ComBat_seq(counts = asv_table, batch = batch_variable, 
                                  group=NULL, covar_mod = covar_mod)
  } else if(!is.na(covars)){
    print(paste0("Defining covariates that should maintain data signature after batch control: ", covars))
    covar_mod <- factor(metadata[,which(colnames(metadata) %in% covars)])
    print("Running batch control.")
    batch_corrected <- ComBat_seq(asv_table, batch_variable, group=covar_mod)
  } else{
    # batch correct with ComBat_seq
    print(paste0("No covariates defined."))
    print("Running batch control.")
    batch_corrected <- ComBat_seq(asv_table, batch_variable, group=NULL)
  }
  
  # visualize how ComBat_seq changed the sample sums
  print("Visualizing how batch control changed data.")
  pre_combatseq <- as.data.frame(colSums(asv_table))
  post_combatseq <- as.data.frame(colSums(batch_corrected))
  
  pre_figure <- ggplot(pre_combatseq, aes(x = `colSums(asv_table)`)) + 
    geom_histogram() + xlab("Pre-Batch Correction Sample Count") + 
    ylab("Count") + theme_bw() 
  post_figure <- ggplot(post_combatseq, aes(x = `colSums(batch_corrected)`)) +
    geom_histogram() + xlab("Post-Batch Correction Sample Count") + 
    ylab("Count") + theme_bw()
  
  multipanel <- grid.arrange(pre_figure, post_figure, ncol = 2, nrow = 1)
  
  ggsave(paste0(yourname, "_", amplicon, "_batchcorrection_change", date, 
                ".png"), multipanel, width = 14, height = 5, units = "in", 
         dpi = 300)
  
  return(batch_corrected)
}

rarefy_data <- function(ps, sample_type, amplicon, yourname, date, threshold){
  print(sample_type)
  # extract asv table from pre-rarefied data
  print("Extracting ASV table from phyloseq object.")
  pre_rare_asv <- as.data.frame(as.matrix(ps@otu_table))
  
  # Save the threshold as the minimum of the sample sums IF either a threshold 
  # isn't set or if the threshold set is less than the minimum of the sample
  # sums
  if(is.na(threshold)){
    threshold <- min(sample_sums(ps))
  } else if(threshold < min(sample_sums(ps))){
    threshold <- min(sample_sums(ps))
  }
  
  # rarefy data
  print("Rarefying data.")
  set.seed(1)
  ps_rare <- rarefy_even_depth(ps, rngseed=1, 
                                 sample.size=threshold, replace=F)
  
  # extract the asv table from the rarefied data
  print("Extracting ASV table from phyloseq object.")
  rare_asv <- as.data.frame(as.matrix(ps_rare@otu_table))
  
  # write to files
  print("Writing rarefied ASV table to file.")
  write.csv(rare_asv, paste0(yourname, "_", amplicon, "_", sample_type, 
                             "_rarefied_ASV_table", date, ".csv"))
  saveRDS(ps_rare, paste0(yourname, "_", amplicon, "_", sample_type, 
                          "_rarefied_phyloseq", date, ".RDS"))
  
  # format data for plotting
  print("Formatting data for plotting.")
  rare_asv$asv <- rownames(rare_asv)
  pre_rare_asv$asv <- rownames(pre_rare_asv)
  
  data_long <- pivot_longer(rare_asv, !asv, names_to = "sample", 
                            values_to = "rarefied_count")
  data_long_pre <- pivot_longer(pre_rare_asv, !asv, names_to = "sample", 
                                values_to = "pre_rarefied_count")
  
  data_compare <- merge(data_long, data_long_pre, by = c("asv", "sample"))
  
  tax <- as.data.frame(as.matrix(ps@tax_table))
  tax$asv <- rownames(tax)
  
  data_compare <- merge(data_compare, tax, by = "asv")
  
  # plot change
  print("Plotting changes due to rarefying in data.")
  by_phylum <- ggplot(data_compare, aes(x = pre_rarefied_count, 
                                        y = rarefied_count, col = Phylum)) + 
    geom_point() + theme_bw()
  
  by_sample <- ggplot(data_compare, aes(x = pre_rarefied_count, 
                                        y = rarefied_count, col = sample)) + 
    geom_point() + theme_bw()
  multipanel <- grid.arrange(by_phylum, by_sample, ncol = 1, nrow = 2)
  ggsave(paste0("Figures/", yourname, "_", amplicon, "_", sample_type, 
                "_pre_post_rarefaction_counts", date, ".png"), multipanel, 
         height = 10, width = 7, units = "in", dpi = 300)
  
  # plot rarefaction curves
  print("Plotting rarefaction curves.")
  if(sample_type == "leaf"){
    png(filename=paste0(yourname, "_", amplicon, "_", sample_type, 
                        "rarefaction_curve", date, ".png"))
    rarecurve(t(as.data.frame(as.matrix(ps@otu_table))), step = 20, 
              sample = 200, col = "blue", label = FALSE)
    dev.off()
  } else{
    png(filename=paste0("Figures/", yourname, "_", amplicon, "_", sample_type, 
                        "rarefaction_curve", date, ".png"))
    rarecurve(t(as.data.frame(as.matrix(ps@otu_table))), step = 20, 
              sample = min(colSums(as.data.frame(as.matrix(ps@otu_table)))), 
              col = "blue", label = FALSE)
    dev.off()
  }
}

clr_data <- function(ps, sample_type, amplicon, yourname, date){
  print(sample_type)
  # extract the ASV table
  print("Extracting ASV table.")
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # do the CLR transformation
  print("CLR-transforming data. ")
  clr_transformed <- as.data.frame(clr(asv_table+1))
  
  # write to file
  print("Writing transformation to file.")
  write.csv(clr_transformed, paste0(yourname, "_", amplicon, "_", sample_type, 
                                    "_CLRtransformed_ASV_table", date, ".csv"))
}

zscore_data <- function(ps, sample_type, amplicon, yourname, date){
  print(sample_type)
  # extract the ASV table
  print("Extracting ASV table.")
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # do the Z-Score transformation
  print("Z-score transforming data.")
  z_transformed <- as.data.frame(scale(asv_table, center = TRUE, scale = TRUE))
  
  # write to file
  print("Writing transformation to file.")
  write.csv(z_transformed, paste0(yourname, "_", amplicon, "_", sample_type, 
                                  "_ZScoretransformed_ASV_table", date, ".csv"))
}

aitchison_data <- function(ps, sample_type, amplicon, yourname, date){
  print(sample_type)
  
  if(class(ps)[[1]]=="phyloseq"){
    # extract the ASV table
    print("Extracting the ASV table.")
    asv_table <- as.data.frame(as.matrix(ps@otu_table))
  } else{
    asv_table <- ps
  }
  
  # calculate Aitchison distance matrix
  print("Calculating Aitchison's distance matrix for data.")
  aitchison_distance_matrix <- aDist(t(asv_table+1), y = NULL)
  aitchison_distance_matrix <- as.matrix(aitchison_distance_matrix)
  
  # write to file
  print("Writing to file.")
  write.csv(aitchison_distance_matrix, paste0(yourname, "_", amplicon, "_", 
                                              sample_type, 
                                              "_Aitchison_distance_matrix", 
                                              date, ".csv"))
}

aggregate_guilds_16S <- function(seq_data, perc_data, taxonomy){
  # merge relative abundance data with taxonomy data
  perc_data <- merge(perc_data, taxonomy, by.x = 0, 
                     by.y = "ASV_ID", all.x = TRUE)
  
  # keep the data we need for aggregating
  pd <- perc_data[,c(2:(ncol(seq_data)+1))]
  pd$guild <- perc_data$guild
  
  # aggregate by guild
  perc_data_guild <- aggregate(.~guild, perc_data, FUN = "sum", na.rm = F)
  return(perc_data_guild)
}

aggregate_guilds_ITS <- function(seq_data, perc_data, taxonomy){
  # merge relative abundance data with taxonomy data
  perc_data <- merge(perc_data, taxonomy, by.x = 0, 
                     by.y = "ASV_ID", all.x = TRUE)
  
  # keep the data we need for aggregating
  pd <- perc_data[,c(2:(ncol(seq_data)+1))]
  pd$primary_lifestyle <- perc_data$primary_lifestyle
  
  # aggregate by guild
  perc_data_guild <- aggregate(.~primary_lifestyle, 
                               pd, FUN = "sum", na.rm = F)
  return(perc_data_guild)
}

record_16s_guild_abund <- function(perc_data_guild){
  data <- perc_data_guild
  rownames(data) <- data$guild
  data <- data[,-1]
  # Archaea
  rows <- grep("Archaea", rownames(data)) 
  subset <- data[rows,]
  perc_archaea <- data.frame(colSums(subset))
  colnames(perc_archaea)[1] <- "perc_archaea"
  
  # C Fixation
  rows <- grep("C fixation", rownames(data)) 
  subset <- data[rows,]
  perc_c_fixation <- data.frame(colSums(subset))
  colnames(perc_c_fixation)[1] <- "perc_c_fixation"
  
  # Cellulolytic
  rows <- grep("Cellulolytic", rownames(data)) 
  subset <- data[rows,]
  perc_cellulolytic <- data.frame(colSums(subset))
  colnames(perc_cellulolytic)[1] <- "perc_cellulolytic"
  
  # Chitinolytic
  rows <- grep("Chitinolytic", rownames(data)) 
  subset <- data[rows,]
  perc_chitinolytic <- data.frame(colSums(subset))
  colnames(perc_chitinolytic)[1] <- "perc_chitinolytic"
  
  # Lignolytic
  rows <- grep("Lignolytic", rownames(data)) 
  subset <- data[rows,]
  perc_lignolytic <- data.frame(colSums(subset))
  colnames(perc_lignolytic)[1] <- "perc_lignolytic"
  
  # Carbon monoxide oxidation
  rows <- grep("Carbon monoxide oxidation", rownames(data)) 
  subset <- data[rows,]
  perc_carbon_monoxide_oxidation <- data.frame(colSums(subset))
  colnames(perc_carbon_monoxide_oxidation)[1] <- "perc_carbon_monoxide_oxidation"
  
  # N_fixation
  rows <- grep("N_fixation", rownames(data)) 
  subset <- data[rows,]
  perc_n_fixation <- data.frame(colSums(subset))
  colnames(perc_n_fixation)[1] <- "perc_n_fixation"
  
  # Copiotroph
  rows <- grep("Copiotroph", rownames(data)) 
  subset <- data[rows,]
  perc_copiotroph <- data.frame(colSums(subset))
  colnames(perc_copiotroph)[1] <- "perc_copiotroph"
  
  # Denitrification
  rows <- grep("Denitrification", rownames(data)) 
  subset <- data[rows,]
  perc_denitrification <- data.frame(colSums(subset))
  colnames(perc_denitrification)[1] <- "perc_denitrification"
  
  # Dissim_nitrate_reduction
  rows <- grep("Dissim_nitrate_reduction", rownames(data)) 
  subset <- data[rows,]
  perc_dissim_nitrate_reduction <- data.frame(colSums(subset))
  colnames(perc_dissim_nitrate_reduction)[1] <- "perc_dissim_nitrate_reduction"
  
  # Hydrocarbon degradation
  rows <- grep("Hydrocarbon degradation", rownames(data)) 
  subset <- data[rows,]
  perc_hydrocarbon_degradation <- data.frame(colSums(subset))
  colnames(perc_hydrocarbon_degradation)[1] <- "perc_hydrocarbon_degradation"
  
  # Sulfonate desulfurization
  rows <- grep("Sulfonate desulfurization", rownames(data)) 
  subset <- data[rows,]
  perc_sulfonate_desulfurization <- data.frame(colSums(subset))
  colnames(perc_sulfonate_desulfurization)[1] <- "perc_sulfonate_desulfurization"
  
  # Phosphate solubilization
  rows <- grep("foliar_endophyte", rownames(data))
  subset <- data[rows,]
  perc_phosphate_solubilization <- data.frame(colSums(subset))
  colnames(perc_phosphate_solubilization)[1] <- "perc_phosphate_solubilization"
  
  # Other N-cycling
  rows <- grep("lichen_parasite", rownames(data)) 
  subset <- data[rows,]
  perc_other_n_cycling <- data.frame(colSums(subset))
  colnames(perc_other_n_cycling)[1] <- "perc_other_n_cycling"
  
  # Other P-cycling
  rows <- grep("Other P-cycling", rownames(data)) 
  subset <- data[rows,]
  perc_other_p_cycling <- data.frame(colSums(subset))
  colnames(perc_other_p_cycling)[1] <- "perc_other_p_cycling"
  
  # Oxidize reduced sulfur
  rows <- grep("Oxidize reduced sulfur", rownames(data)) 
  subset <- data[rows,]
  perc_oxidize_reduced_sulfur <- data.frame(colSums(subset))
  colnames(perc_oxidize_reduced_sulfur)[1] <- "perc_oxidize_reduced_sulfur"
  
  # Methanotroph
  rows <- grep("Methanotroph", rownames(data)) 
  subset <- data[rows,]
  perc_methanotroph <- data.frame(colSums(subset))
  colnames(perc_methanotroph)[1] <- "perc_methanotroph"
  
  # Oligotroph
  rows <- grep("Oligotroph", rownames(data)) 
  subset <- data[rows,]
  perc_oligotroph <- data.frame(colSums(subset))
  colnames(perc_oligotroph)[1] <- "perc_oligotroph"
  
  # Partial_Nitrification
  rows <- grep("Partial_Nitrification", rownames(data))
  subset <- data[rows,]
  perc_partial_nitrification <- data.frame(colSums(subset))
  colnames(perc_partial_nitrification)[1] <- "perc_partial_nitrification"
  
  result_guilds <- cbind(perc_c_fixation, perc_cellulolytic, 
                           perc_chitinolytic, perc_lignolytic, 
                           perc_carbon_monoxide_oxidation, perc_n_fixation, 
                           perc_copiotroph, perc_denitrification, 
                           perc_dissim_nitrate_reduction, 
                           perc_hydrocarbon_degradation, 
                           perc_sulfonate_desulfurization, 
                           perc_phosphate_solubilization, 
                           perc_other_n_cycling, perc_other_p_cycling, 
                           perc_oxidize_reduced_sulfur, perc_methanotroph, 
                           perc_oligotroph, perc_partial_nitrification, 
                           perc_archaea)
  
  return(result_guilds)
}

record_its_guild_abund <- function(perc_data_guild, perc_data, taxonomy){
  data <- perc_data_guild
  rownames(data) <- data$primary_lifestyle
  data <- data[,-1]
  
  # Algal Ectosymbiont
  rows <- grep("algal_ectosymbiont", rownames(data)) 
  subset <- data[rows,]
  algal_ecto <- data.frame(colSums(subset))
  colnames(algal_ecto)[1] <- "algal_ecto"
  
  # Algal Parasite
  rows <- grep("algal_parasite", rownames(data)) 
  subset <- data[rows,]
  algal_para <- data.frame(colSums(subset))
  colnames(algal_para)[1] <- "algal_para"
  
  # Algivorous/Protistivorous
  rows <- grep("algivorous/protistivorous", rownames(data)) 
  subset <- data[rows,]
  algivorous <- data.frame(colSums(subset))
  colnames(algivorous)[1] <- "algivorous"
  
  # Animal Parasite
  animal_cap <- merge(perc_data, taxonomy, by.x = 0, by.y = "ASV_ID")
  rownames(animal_cap) <- animal_cap[,1]
  animal_cap <- animal_cap[,-1]
  animal_cap <- animal_cap[,c(1:ncol(perc_data), (ncol(perc_data) + 4))]
  animal_cap <- na.omit(animal_cap)
  animal_cap[animal_cap$Animal_biotrophic_capacity_template == "",] <-NA
  animal_cap <- na.omit(animal_cap)
  rows <- grep("parasite", animal_cap$Animal_biotrophic_capacity_template)
  subset <- animal_cap[rows,]
  subset <- subset[,-ncol(subset)]
  animal_parasite <- data.frame(colSums(subset))
  colnames(animal_parasite) [1] <- "Animal_parasite"
  
  rows <- grep(paste(c("animal_parasite", "animal-associated"), collapse="|"), 
                animal_cap$Animal_biotrophic_capacity_template)
  subset <- animal_cap[rows,]
  subset <- subset[,-ncol(subset)]
  animal_pathogen <- data.frame(colSums(subset))
  colnames(animal_pathogen) [1] <- "Animal_associated_pathogen"
  
  rows <- grep(c("animal_endosymbiont"), animal_cap$Animal_biotrophic_capacity_template)
  subset <- animal_cap[rows,]
  subset <- subset[,-ncol(subset)]
  animal_endosymbiont <- data.frame(colSums(subset))
  colnames(animal_endosymbiont) [1] <- "Animal_endosymbiont"
  
  rows <- grep("human", animal_cap$Animal_biotrophic_capacity_template)
  subset <- animal_cap[rows,]
  subset <- subset[,-ncol(subset)]
  human_pathogen <- data.frame(colSums(subset))
  colnames(human_pathogen) [1] <- "Human_pathogen"
  
  rows <- grep("opportunistic", animal_cap$Animal_biotrophic_capacity_template)
  subset <- animal_cap[rows,]
  subset <- subset[,-ncol(subset)]
  opportunistic_pathogen <- data.frame(colSums(subset))
  colnames(opportunistic_pathogen) [1] <- "Opportunistic_pathogen"
  
  # Arbuscular Mycorrhizal
  rows <- grep("arbuscular_mycorrhizal", rownames(data)) 
  subset <- data[rows,]
  amf <- data.frame(colSums(subset))
  colnames(amf)[1] <- "amf"
  
  # Arthropod Associated
  rows <- grep("arthropod-associated", rownames(data)) 
  subset <- data[rows,]
  arthropod <- data.frame(colSums(subset))
  colnames(arthropod)[1] <- "arthropod"
  
  # Bacterivorous
  rows <- grep("bacterivorous", rownames(data)) 
  subset <- data[rows,]
  bacterivorous <- data.frame(colSums(subset))
  colnames(bacterivorous)[1] <- "bacterivorous"
  
  # Dung Saprotroph
  rows <- grep("dung_saprotroph", rownames(data)) 
  subset <- data[rows,]
  dung <- data.frame(colSums(subset))
  colnames(dung)[1] <- "dung_sap"
  
  # ECM
  rows <- grep("ectomycorrhizal", rownames(data)) 
  subset <- data[rows,]
  ECM <- data.frame(colSums(subset))
  colnames(ECM)[1] <- "ECM"
  
  # Epiphyte
  rows <- grep("epiphyte", rownames(data)) 
  subset <- data[rows,]
  epiphyte <- data.frame(colSums(subset))
  colnames(epiphyte)[1] <- "epiphyte"
  
  # Foliar Endophyte
  rows <- grep("foliar_endophyte", rownames(data))
  subset <- data[rows,]
  foliar <- data.frame(colSums(subset))
  colnames(foliar)[1] <- "foliar_endophyte"
  
  # Lichen Parasite
  rows <- grep("lichen_parasite", rownames(data)) 
  subset <- data[rows,]
  lichen_parasite <- data.frame(colSums(subset))
  colnames(lichen_parasite)[1] <- "lichen_parasite"
  
  # Lichenized
  rows <- grep("lichenized", rownames(data)) 
  subset <- data[rows,]
  lichenized <- data.frame(colSums(subset))
  colnames(lichenized)[1] <- "lichenized"
  
  # Litter Saprotroph
  rows <- grep("litter_saprotroph", rownames(data)) 
  subset <- data[rows,]
  litter_sap <- data.frame(colSums(subset))
  colnames(litter_sap)[1] <- "litter_sap"
  
  # Moss Symbiont
  rows <- grep("moss_symbiont", rownames(data)) 
  subset <- data[rows,]
  moss_symb <- data.frame(colSums(subset))
  colnames(moss_symb)[1] <- "moss_symb"
  
  # Mycoparasite
  rows <- grep("mycoparasite", rownames(data)) 
  subset <- data[rows,]
  mycoparasite <- data.frame(colSums(subset))
  colnames(mycoparasite)[1] <- "mycoparasite"
  
  # Nectar/Tap Saprotroph
  rows <- grep("nectar/tap_saprotroph", rownames(data)) 
  subset <- data[rows,]
  nectar_sap <- data.frame(colSums(subset))
  colnames(nectar_sap)[1] <- "nectar_sap"
  
  # Plant Pathogen
  plant_cap <- merge(perc_data, taxonomy, by.x = 0, by.y = "ASV_ID")
  rownames(plant_cap) <- plant_cap[,1]
  plant_cap <- plant_cap[,-1]
  plant_cap <- plant_cap[,c(1:ncol(perc_data), (ncol(perc_data) + 3))]
  plant_path_cap <- na.omit(plant_cap)
  plant_path_cap[plant_path_cap$Plant_pathogenic_capacity_template == "",] <-NA
  plant_path_cap <- na.omit(plant_path_cap)
  subset <- plant_path_cap[,-ncol(plant_path_cap)]
  Plant_pathogenic_capacity <- data.frame(colSums (subset))
  colnames(Plant_pathogenic_capacity) [1] <- "Plant_pathogenic_capacity"
  
  rows <- grep("root", plant_path_cap$Plant_pathogenic_capacity_template)
  subset <- plant_path_cap[rows,]
  subset <- subset[,-ncol(subset)]
  root_pathogen <- data.frame(colSums(subset))
  colnames(root_pathogen) [1] <- "Root_associated_pathogen"
  
  rows <- grep("leaf", plant_path_cap$Plant_pathogenic_capacity_template)
  subset <- plant_path_cap[rows,]
  subset <- subset[,-ncol(subset)]
  leaf_pathogen <- data.frame(colSums(subset))
  colnames(leaf_pathogen) [1] <- "Leaf_fruit_seed_pathogen"
  
  rows <- grep("wood", plant_path_cap$Plant_pathogenic_capacity_template)
  subset <- plant_path_cap[rows,]
  subset <- subset[,-ncol(subset)]
  wood_pathogen <- data.frame(colSums(subset))
  colnames(wood_pathogen) [1] <- "Wood_pathogen"
  
  # Pollen Saprotroph
  rows <- grep("pollen_saprotroph", rownames(data)) 
  subset <- data[rows,]
  pollen_sap <- data.frame(colSums(subset))
  colnames(pollen_sap)[1] <- "pollen_sap"
  
  # Protistan Parasite
  rows <- grep("protistan_parasite", rownames(data)) 
  subset <- data[rows,]
  protistan <- data.frame(colSums(subset))
  colnames(protistan)[1] <- "protistan"
  
  # Root Endophyte
  rows <- grep("root_endophyte", rownames(data)) 
  subset <- data[rows,]
  root_endo <- data.frame(colSums(subset))
  colnames(root_endo)[1] <- "root_endo"
  
  # Soil Sapotroph
  rows <- grep("soil_saprotroph", rownames(data)) 
  subset <- data[rows,]
  soil_sap <- data.frame(colSums(subset))
  colnames(soil_sap)[1] <- "soil_sap"
  
  # Sooty Mold
  rows <- grep("sooty_mold", rownames(data)) 
  subset <- data[rows,]
  sooty_mold <- data.frame(colSums(subset))
  colnames(sooty_mold)[1] <- "sooty_mold"
  
  # Wood Saprotroph
  rows <- grep("wood_saprotroph", rownames(data)) 
  subset <- data[rows,]
  wood_sap <- data.frame(colSums(subset))
  colnames(wood_sap)[1] <- "wood_sap"
  
  # Saprotroph
  rows <- grep("saprotroph", rownames(data)) 
  subset <- data[rows,]
  saps_all <- data.frame(colSums(subset))
  colnames(saps_all)[1] <- "saps_all"
  
  # Pathogens
  rows <- grep("patho", rownames(data)) 
  subset <- data[rows,]
  pathos_all <- data.frame(colSums(subset))
  colnames(pathos_all)[1] <- "pathos_all"
  
  # Symbionts 
  rows <- grep(c("symbio"), rownames(data))
  subset <- data[rows,]
  symbionts_all <- data.frame(colSums(subset)) + animal_endosymbiont$Animal_endosymbiont
  colnames(symbionts_all)[1] <- "symbionts_all"
  
  result <- cbind(algal_ecto, algal_para, algivorous, animal_endosymbiont, 
                  animal_parasite, amf, arthropod, bacterivorous, dung, 
                  ECM, epiphyte, foliar, lichen_parasite, lichenized, 
                  litter_sap, moss_symb, mycoparasite, nectar_sap, 
                  Plant_pathogenic_capacity, pollen_sap, protistan, 
                  root_endo, soil_sap, sooty_mold, wood_sap, saps_all, 
                  pathos_all, symbionts_all, animal_pathogen, 
                  human_pathogen, opportunistic_pathogen, root_pathogen, 
                  leaf_pathogen, wood_pathogen)
  
  return(result)
}

process_16s_pathogens <- function(pathogens_file, pathogens_tax, taxonomy){
  print("Processing the 16S pathogens output files.")
  # find unique query ASV IDs
  asv_ids <- unique(pathogens$asv_id)
  
  # make a list of the rows to keep
  print("Making a list of the pathogen matches to keep")
  keep <- c()
  for(i in 1:length(asv_ids)){
    asv <- asv_ids[i]
    asv_matches <- pathogens[which(pathogens$asv_id == asv),]
    perc_id <- max(asv_matches$perc_id)
    keep <- c(keep, asv_matches$row[which(asv_matches$perc_id == perc_id)])
  }
  
  # keep the best matches for each unique ASV ID
  print("Keeping the best pathogen matches")
  pathogens <- pathogens[keep,]
  
  # remove the Kingdom information from the pathogen ID in the taxonomy file
  pathogens_tax$path_id <- gsub("\td__Bacteria", "", pathogens_tax$path_id)
  
  # rename the first column
  colnames(pathogens_tax)[1] <- "pathogen_id"
  
  # Keep the useful columns 
  pathogens_tax <- pathogens_tax[,c(1,7,8)]
  
  # Merge the pathogen taxonomy to the taxonomy information
  print("Merging pathogen type with ASV taxonomy")
  pathogens <- merge(pathogens_tax, pathogens, by = "pathogen_id", all.y = TRUE)
  
  # Keep only the useful columns
  pathogens <- pathogens[,c(2:4)]
  
  # Keep distinct rows
  pathogens <- distinct(pathogens)
  
  # find the unique ASV ids from the pathogens matched dataframe
  asvs <- unique(pathogens$asv_id)
  
  # create a new column for recording multiple pathogen types
  pathogens$pathogen_type <- NA
  
  # format the pathogen type column from the MBPD
  pathogens$type <- gsub("t__", "", pathogens$type)
  
  # collapse multiple pathogen types into one string
  print("Collapsing multiple pathogen types into one string")
  for(i in 1:length(asvs)){
    asv <- asvs[i]
    type <- unique(pathogens$type[which(pathogens$asv_id == asv)])
    pathogens$pathogen_type[which(pathogens$asv_id == asv)] <- paste(type, 
                                                                     collapse = ", ")
  }
  pathogens <- pathogens[,c(3,4)]
  pathogens <- distinct(pathogens)
  
  colnames(pathogens)[1] <- "ASV_ID"
  taxonomy <- merge(taxonomy, pathogens, all.x = TRUE)
  
  return(taxonomy)
}

aggregate_16s_pathogens <- function(seq_data, perc_data, taxonomy){
  # merge relative abundance data with taxonomy data
  perc_data <- merge(perc_data, taxonomy, by.x = 0, 
                     by.y = "ASV_ID", all.x = TRUE)
  
  # keep the data we need for aggregating
  perc_data <- perc_data[,c(2:(ncol(seq_data)+1), ncol(perc_data))]
  
  # aggregate by guild
  perc_data_guild <- aggregate(.~pathogen_type, perc_data, FUN = "sum", 
                               na.rm = F)
  return(perc_data_guild)
}

record_16s_pathogen_abund <- function(perc_data_guild){
  data <- perc_data_guild
  rownames(data) <- data$pathogen_type
  data <- data[,-1]
  
  # Plant Pathogen
  rows <- grep("Plant", rownames(data)) 
  subset <- data[rows,]
  perc_plant_pathogen <- data.frame(colSums(subset))
  colnames(perc_plant_pathogen)[1] <- "perc_plant_pathogen"
  
  # Animal Pathogen
  rows <- grep("Animal", rownames(data)) 
  subset <- data[rows,]
  perc_animal_pathogen <- data.frame(colSums(subset))
  colnames(perc_animal_pathogen)[1] <- "perc_animal_pathogen"
  
  # Zoonotic Pathogen
  rows <- grep("Zoonotic", rownames(data)) 
  subset <- data[rows,]
  perc_zoonotic_pathogen <- data.frame(colSums(subset))
  colnames(perc_zoonotic_pathogen)[1] <- "perc_zoonotic_pathogen"
  
  result <- cbind(perc_plant_pathogen, perc_animal_pathogen, 
                  perc_zoonotic_pathogen)
  
  return(result)
}
