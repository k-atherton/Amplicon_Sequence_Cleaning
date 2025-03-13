library(dplyr)

leaf_16s <- read.csv("Data_Cleaning/Rarefaction/16S/leaf_16s_asv_rare_wtax_guild_200.csv")
m_16s <- read.csv("Data_Cleaning/Rarefaction/16S/m_16s_asv_rare_wtax_guild.csv")
o_16s <- read.csv("Data_Cleaning/Rarefaction/16S/o_16s_asv_rare_wtax_guild.csv")
root_16s <- read.csv("Data_Cleaning/Rarefaction/16S/root_16s_asv_rare_wtax_guild.csv")

leaf_its <- read.csv("Data_Cleaning/Rarefaction/ITS/leaf_its_asv_rare_wtax_guild.csv")
m_its <- read.csv("Data_Cleaning/Rarefaction/ITS/m_its_asv_rare_wtax_guild.csv")
o_its <- read.csv("Data_Cleaning/Rarefaction/ITS/o_its_asv_rare_wtax_guild.csv")
root_its <- read.csv("Data_Cleaning/Rarefaction/ITS/root_its_asv_rare_wtax_guild.csv")

metadata_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_16s_200_new.csv")
metadata_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_its_new.csv")

metadata_16s <- metadata_16s[which(metadata_16s$site_name == "Harvard Forest"),]
metadata_its <- metadata_its[which(metadata_its$site_name == "Harvard Forest"),]

leaf_16s_hf <- leaf_16s[,colnames(leaf_16s) %in% metadata_16s$sample_name]
leaf_16s_hf <- cbind(leaf_16s[,c(3:9, ncol(leaf_16s))], leaf_16s_hf)

m_16s_hf <- m_16s[,colnames(m_16s) %in% metadata_16s$sample_name]
m_16s_hf <- cbind(m_16s[,c(3:9, ncol(m_16s))], m_16s_hf)

o_16s_hf <- o_16s[,colnames(o_16s) %in% metadata_16s$sample_name]
o_16s_hf <- cbind(o_16s[,c(3:9, ncol(o_16s))], o_16s_hf)

root_16s_hf <- root_16s[,colnames(root_16s) %in% metadata_16s$sample_name]
root_16s_hf <- cbind(root_16s[,c(3:9, ncol(root_16s))], root_16s_hf)


leaf_its_hf <- leaf_its[,colnames(leaf_its) %in% metadata_its$sample_name]
leaf_its_hf <- cbind(leaf_its[,c(6:10,2,11,3:5)], leaf_its_hf)
colnames(leaf_its_hf)[1:10] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Fungal_Functional_Guild", "Plant_pathogenic_capacity", "Animal_biotrophic_capacity")

m_its_hf <- m_its[,colnames(m_its) %in% metadata_its$sample_name]
m_its_hf <- cbind(m_its[,c(6:10,2,11,3:5)], m_its_hf)
colnames(m_its_hf)[1:10] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Fungal_Functional_Guild", "Plant_pathogenic_capacity", "Animal_biotrophic_capacity")

o_its_hf <- o_its[,colnames(o_its) %in% metadata_its$sample_name]
o_its_hf <- cbind(o_its[,c(6:10,2,11,3:5)], o_its_hf)
colnames(o_its_hf)[1:10] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Fungal_Functional_Guild", "Plant_pathogenic_capacity", "Animal_biotrophic_capacity")


root_its_hf <- root_its[,colnames(root_its) %in% metadata_its$sample_name]
root_its_hf <- cbind(root_its[,c(6:10,2,11,3:5)], root_its_hf)
colnames(root_its_hf)[1:10] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Fungal_Functional_Guild", "Plant_pathogenic_capacity", "Animal_biotrophic_capacity")

write.csv(leaf_16s_hf, "Data_Cleaning/Rarefaction/16S/leaf_16s_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(m_16s_hf, "Data_Cleaning/Rarefaction/16S/soil_15-30cm_16s_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(o_16s_hf, "Data_Cleaning/Rarefaction/16S/soil_0-15cm_16s_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(leaf_16s_hf, "Data_Cleaning/Rarefaction/16S/root_16s_asv_rare_wtax_guild_hf.csv", row.names = F)


write.csv(leaf_its_hf, "Data_Cleaning/Rarefaction/ITS/leaf_its_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(m_its_hf, "Data_Cleaning/Rarefaction/ITS/soil_15-30cm_its_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(o_its_hf, "Data_Cleaning/Rarefaction/ITS/soil_0-15cm_its_asv_rare_wtax_guild_hf.csv", row.names = F)
write.csv(leaf_its_hf, "Data_Cleaning/Rarefaction/ITS/root_its_asv_rare_wtax_guild_hf.csv", row.names = F)
