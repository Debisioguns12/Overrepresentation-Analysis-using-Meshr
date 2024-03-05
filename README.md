# Overrepresentation-Analysis-using-Meshr
Meshr is used to confirm overrepresentation of the genes 
# MeSH (Medical Subject Headings) over-representation analysis
#Setting working directory
setwd("/Users/ogunbawoadebisi/Downloads/meta_analyse")

# Check if BiocManager package is installed, if not, install it
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install("biomaRt")
BiocManager::install("org.Bt.eg.db")
BiocManager::install("MeSH.db")
BiocManager::install("meshes")
BiocManager::install("MeSH.Bta.eg.db")

# Load required libraries
library(biomaRt)  # For accessing biological data through various databases
library(meshes)   # For MeSH enrichment analysis
library(readxl)   # For reading Excel files
library(AnnotationHub)  # For accessing and retrieving biological metadata

# Access AnnotationHub to query MeSHDb and B. Taurus data
ah <- AnnotationHub()
qr_hsa <- query(ah, c("MeSHDb", "Bos Taurus"))  # Querying MeSHDb and B. Taurus
filepath_hsa <- qr_hsa[[1]]  # Retrieve the path to the MeSH database file
db <- MeSHDbi::MeSHDb(filepath_hsa)  # Access the MeSH database

# Import dataset containing gene information
genes.modules <- read.table("Geneid_ra.csv",h = TRUE, sep = "\t", stringsAsFactors = FALSE)
#here in the dataset we need the genes in the ENSEMBL format: 
# ex:
# ENSBTAG00000011808
# ENSBTAG00000011873
# ENSBTAG00000013210
# ENSBTAG00000013578
# ENSBTAG00000015086
# ENSBTAG00000015204
# ENSBTAG00000018810
# Import dataset containing gene information
genes.modules <- read_xlsx("Geneid_ra.xlsx")

# Get Entrez IDs for the genes using biomaRt package
mart <- useMart("ENSEMBL_MART_ENSEMBL")  # Select the ENSEMBL database
mart <- useDataset("btaurus_gene_ensembl", mart)  # Choose the B. Taurus dataset
filter <- "ensembl_gene_id"  # Filter for ENSEMBL gene IDs
value <- unique(genes.modules$gene_id)  # Unique gene names from the dataset 
attributes <- c("ensembl_gene_id", "external_gene_name", "entrezgene_id")  # Attributes to retrieve
all.genes <- getBM(attributes = attributes, filters = filter, values = value, mart = mart)  # Get gene information

# Perform MeSH enrichment analysis for Anatomy (category A)
mesh.A <- enrichMeSH(all.genes[which(!is.na(all.genes$entrezgene_id)), "entrezgene_id"],
                     MeSHDb = db, database = 'gene2pubmed',
                     category = 'A',
                     pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
meshA.result <- mesh.A@result  # Store results
meshA.result$group <- "Anatomy"  # Add group label
head(mesh.A)
# Perform MeSH enrichment analysis for Diseases (category C)
mesh.C <- enrichMeSH(all.genes[which(!is.na(all.genes$entrezgene_id)), "entrezgene_id"],
                     MeSHDb = db, database = 'gene2pubmed',
                     category = 'C',
                     pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
meshC.result <- mesh.C@result  # Store results
meshC.result$group <- "Diseases"  # Add group label

# Perform MeSH enrichment analysis for Biological Sciences (category G)
mesh.G <- enrichMeSH(all.genes[which(!is.na(all.genes$entrezgene_id)), "entrezgene_id"],
                     MeSHDb = db, database = 'gene2pubmed',
                     category = 'G', pvalueCutoff = 1,
                     qvalueCutoff = 1, minGSSize = 1)
meshG.result <- mesh.G@result  # Store results
meshG.result$group <- "Biological Sciences"  # Add group label

# Combine results into a single dataframe
mesh.final <- rbind(meshA.result, meshC.result, meshG.result)

# Display the first few rows of the combined dataframe
head(mesh.final)

# Write the allgenes data to an Excel file
library(readxl)
library(openxlsx)
write.xlsx(all.genes, "all_genes_ra.xlsx")
# Write the meshfinal data to an Excel file
write.xlsx(mesh.final, "mesh_final_ra.xlsx")
