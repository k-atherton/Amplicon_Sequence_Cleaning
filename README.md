# Amplicon Sequence Cleaning Workflow
by Kathryn Atherton

The purpose of this project is to document the amplicon sequence cleaning process and decisions that I use in my research so that others can reproduce this process on my data or on their own datasets. 

My workflow starts after sequences have been clustered and annotated by DADA2 (I use the [BU16S pipeline](https://github.com/Boston-University-Microbiome-Initiative/BU16s) written by Dr. Michael Silverstein). It guides a user through testing for batch effect and sequence depth influence, identifying and dropping outliers, removing contaminants from samples based on negative control sequence data, and normalizing data appropriately for various downstream analyses, including rarefaction, CLR and z-score transformations, and calculating Aitchison distance. The workflow also allows for annotating fungal and bacterial functional guilds.

## Requirements
To run this pipeline, you will need the following R packages:

