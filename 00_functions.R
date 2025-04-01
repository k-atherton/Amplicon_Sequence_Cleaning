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

id_outliers_evaluate_seq_depth <- function(data, metadata, sample_type, yourname, 
                                           amplicon, date, rep){
  # make a transposed verison of the sequencing data
  print("Transposing the sequencing data:")
  data_t <- as.data.frame(t(data))
  
  # calculate Aitchison distance for samples
  print("Calculating the Aitchison's distance for the samples")
  aitch_data <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch_data <- as.matrix(aitch_data)
  
  # test for sequencing batch and depth effect
  print(paste0("Testing for sequencing batch effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.factor(metadata$sequencing_batch)))
  print(paste0("Testing for sequencing depth effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.numeric(metadata$seq_count_dada2)))
  
  print("Calculating NMDS for plotting.")
  all.MDS <- metaMDS(aitch_data,k=2,zerodist="add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)
  
  # plot samples by different variables
  print("Plotting samples by different variables")
  ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2)) +
    geom_text(label = rownames(coordinates),
              size = 2) +
    stat_ellipse() +
    theme_bw() + expand_limits(x = c(min(coordinates$NMDS1) - 10, 
                                     max(coordinates$NMDS1) + 10))
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
  print(paste0("Keeping samples with dada2 read count > ", threshold, "."))
  metadata_drop <- metadata[metadata$seq_count_dada2 > threshold,]
  data_drop <- data[,colnames(data) %in% rownames(metadata_drop)]
  
  # make a transposed verison of the sequencing data
  print("Transposing the sequencing data.")
  data_t <- as.data.frame(t(data_drop))
  
  # calculate Aitchison distance for samples
  print("Calculating Aitchison's distance for samples.")
  aitch <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print("Testing sequencing batch effect:")
  print(adonis2(aitch ~ as.factor(metadata_drop$sequencing_batch)))
  
  print("Testing sequencing depth effect:")
  adonis2(aitch ~ as.numeric(metadata_drop$seq_count_dada2))
  
  # calculate NMDS
  print("Calculating NMDS for plotting purposes.")
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates<-data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata_drop)
  
  # plot samples by different variables
  print("Plotting samples by different variables.")
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
                "_NMDS_drop", threshold, date, ".png"), multipanel, 
         width = 21, height = 10, units = "in", dpi = 300)
  
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
  if(sample_type == "leaf"){
    final_data <- subset_samples(final_data, sample_type == "Leaf")
  } else if(sample_type == "root"){
    final_data <- subset_samples(final_data, sample_type == "Root")
  } else {
    final_data <- subset_samples(final_data, sample_type == "Soil")
  }
  
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
                              "ASV_table_decontam", date, ".csv"))
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

check_batch_effect <- function(ps, amplicon, yourname, date){
  # extract asv table
  print("Extract ASV table from phyloseq object.")
  sequence_data <- as.data.frame(as.matrix(ps@otu_table))
  
  # extract metadata
  print("Extract metadata from phyloseq object.")
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  metadata$sample_type[which(metadata$sample_type == "Soil")] <- paste0(metadata$soil_horizon[which(metadata$sample_type == "Soil")], 
                                                                        " Soil")
  
  # make a transposed version of asv table
  print("Make transposed version of ASV table.")
  sequence_data_t <- as.data.frame(t(sequence_data))
  
  # calculate Aitchison distance for samples
  print("Calculate Aitchison distance for samples.")
  aitch <- aDist(sequence_data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print("Calculating batch effect:")
  print(adonis2(aitch ~ as.factor(metadata$sequencing_batch)))
  print("Calculating tree type effect:")
  print(adonis2(aitch ~ as.factor(metadata$tree_type)))
  print("Calculating tree age effect:")
  print(adonis2(aitch ~ as.factor(metadata$tree_age)))
  print("Calculating tree pit type effect:")
  print(adonis2(aitch ~ as.factor(metadata$tree_pit_type)))
  print("Calculating tree species effect:")
  print(adonis2(aitch ~ as.factor(metadata$tree_species)))
  print("Calculating site effect:")
  print(adonis2(aitch ~ as.factor(metadata$site_name)))
  print("Calculating sample type effect:")
  print(adonis2(aitch ~ as.factor(metadata$sample_type)))
  
  print("Calculating NMDS for visualization purposes.")
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)
  
  # plot samples by different variables
  print("Plotting samples by different variables.")
  seq_batch <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                       col = as.factor(sequencing_batch))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                       col = as.factor(tree_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_age <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                      col = as.factor(tree_age))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_pit_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                           col = as.factor(tree_pit_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_species <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                          col = as.factor(tree_species))) + 
    geom_point() + stat_ellipse() + theme_bw()
  site_name <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                       col = as.factor(site_name))) + 
    geom_point() + stat_ellipse() + theme_bw()
  sample_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                         col = as.factor(sample_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  
  multipanel <- grid.arrange(seq_batch, tree_type, tree_age, tree_pit_type,
                             tree_species, site_name, sample_type, nrow = 2, 
                             ncol = 4)
  
  ggsave(paste0(yourname, "_", amplicon, "_check_batch_effect", date, ".png"), 
         multipanel, width = 28, height = 10, units = "in", dpi = 300)
}

batch_correct <- function(ps, amplicon, yourname, date){
  # get asv table from phyloseq object
  print("Extracting ASV table from phyloseq object.")
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # remove any taxa with 0 counts and any samples under 1000 counts as 
  # ComBat_seq will not work if you have samples with 10k+ counts AND really 
  # low counts in the same dataset
  print("Removing taxa with 0 counts and samples with < 1000 counts.")
  asv_table <- asv_table[which(rowSums(asv_table) > 0), 
                         which(colSums(asv_table) > 1000)]
  asv_table <- asv_table[,order(colnames(asv_table))]
  asv_table <- as.matrix(asv_table)
  
  # get metadata from phyloseq object
  print("Getting metadata from phyloseq object.")
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  
  # keep metadata samples which are still in the asv table
  print("Keeping metadata samples that are still in ASV table.")
  metadata <- metadata[which(rownames(metadata) %in% colnames(asv_table)),]
  
  # order the metadata by sample name
  print("Ordering the metadata by sample name (same as the ASV table)")
  metadata <- metadata[order(row.names(metadata)),]
  
  # batch is defined as the sample type, as the samples were extracted and 
  # amplified with different methods, but the sequencing batch did not have
  # a batch effect
  print("Defining batch as sample type.")
  batch <- paste0(metadata$sample_type, " ", metadata$soil_horizon)
  # for leaves, remove the soil horizon NA
  batch <- gsub(" NA", "", batch) 
  # for roots, just keep roots as they were all amplified the same way
  batch <- gsub("Root.*", "Root", batch) 
  
  # the covariates we want to maintain the signature of are tree age and tree 
  # pit type
  print("Defining covariates that should maintain data signature after batch control: tree age and tree pit type.")
  age <- metadata$tree_age
  type <- metadata$tree_pit_type
  covar_mod <- data.frame(age, type)
  
  # batch correct with ComBat_seq
  print("Running batch control.")
  batch_corrected <- ComBat_seq(asv_table, batch, group=NULL, covar_mod)
  
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

rarefy_data <- function(ps, sample_type, amplicon, yourname, date){
  print(sample_type)
  # extract asv table from pre-rarefied data
  print("Extracting ASV table from phyloseq object.")
  pre_rare_asv <- as.data.frame(as.matrix(ps@otu_table))
  
  # rarefy data
  print("Rarefying data.")
  set.seed(1)
  if(sample_type == "leaf"){
    ps_rare <- rarefy_even_depth(ps, rngseed=1, sample.size=200, replace=F)
  }
  else{
    ps_rare <- rarefy_even_depth(ps, rngseed=1, 
                                 sample.size=min(sample_sums(ps)), replace=F)
  }
  
  # extract the asv table from the rarefied data
  print("Extracting ASV table from phyloseq object.")
  rare_asv <- as.data.frame(as.matrix(ps_rare@otu_table))
  
  # write to files
  print("Writing rarefied ASV table to file.")
  write.csv(rare_asv, paste0(yourname, "_", amplicon, "_", sample_type, 
                             "_rarefied_ASV_table", date, ".csv"))
  
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
  # extract the ASV table
  print("Extracting the ASV table.")
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
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
