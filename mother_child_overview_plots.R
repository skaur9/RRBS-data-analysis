# Load the necessary library
library(gplots)


setwd("/home/simran/STENO_projects/DNBC_Multi_Omics/Methylation_data/Results/")


################### Child data overview ######################
child_data=read.table("120_child_pheno.txt", sep="\t", h=T)

# Create a subset of your data containing only the relevant columns
subset_data <- child_data[, c("Genotype", "Lipidomics", "Methylation", "Phenotype", "Trimester")]

# Create a matrix to represent the omics data
omics_matrix <- matrix(0, nrow = nrow(subset_data), ncol = ncol(subset_data) - 2)

# Fill in the matrix with 1 for non-missing data
for (i in 1:nrow(subset_data)) {
  for (j in 1:(ncol(subset_data) - 2)) {
    if (subset_data[i, j] == 1) {
      omics_matrix[i, j] <- j  # Assign column index as a unique value for each category
    }
  }
}

# Set row names to sample IDs and column names to omics categories
rownames(omics_matrix) <- rownames(subset_data)
colnames(omics_matrix) <- colnames(subset_data)[1:3]  # Exclude Phenotype and Trimester

# Create a color palette for the heatmap
colors <- c("grey", "steelblue", "springgreen4", "mediumpurple")

# Define colors for Genotype, Lipidomics, and Methylation
category_colors <- c("Genotype" = "steelblue", "Lipidomics" = "mediumpurple", "Methylation" = "springgreen4")

# Create a vector of colors for the columns based on the omics category
col_colors <- sapply(colnames(omics_matrix), function(category) category_colors[category])

# Define colors for Phenotype and Trimester
phenotype_colors <- c("MISSING" = "grey", "CONTROL" = "blue", "CASE" = "red")
trimester_colors <- c("1" = "yellow", "2" = "orange", "3" = "purple")

# Create vectors of colors for Phenotype and Trimester columns
col_phenotype <- phenotype_colors[as.character(subset_data$Phenotype)]
col_trimester <- trimester_colors[as.character(subset_data$Trimester)]

#pdf("DNBC_Methylation_Child.pdf", h=8, w=12)
tiff(filename="DNBC_Methylation_Child.tiff", width=14, height=8, units="in",res=300)
# Create the rotated heatmap with ColSideColors
heatmap.2(t(omics_matrix),  # Transpose the matrix to rotate the plot
          Colv = NA, # To suppress column clustering
          Rowv = NA, # To suppress row clustering
          dendrogram = "none", # No dendrogram
          col = colors,
          cexRow=0.7,
          key = TRUE, key.title = "Omics Category",
          key.xlab = "Grey: Missing, Blue: Genotype, Purple: Lipidomics, Green: Methylation",
          main =  "Multi-Omics DNBC Overview ",
          trace = "none", # No trace
          density.info = "none", # No density plot
          ColSideColors = col_phenotype)  # Assign colors to Phenotype, Trimester, and omics categories

dev.off()




################## Mother data overview########################3

mother_data=read.table("120_mother_pheno.txt", sep="\t", h=T)

# Create a subset of your data containing only the relevant columns
subset_data <- mother_data[, c("Genotype", "Lipidomics", "Methylation", "Phenotype", "Trimester")]

# Create a matrix to represent the omics data
omics_matrix <- matrix(0, nrow = nrow(subset_data), ncol = ncol(subset_data) - 2)

# Fill in the matrix with 1 for non-missing data
for (i in 1:nrow(subset_data)) {
  for (j in 1:(ncol(subset_data) - 2)) {
    if (subset_data[i, j] == 1) {
      omics_matrix[i, j] <- j  # Assign column index as a unique value for each category
    }
  }
}

# Set row names to sample IDs and column names to omics categories
rownames(omics_matrix) <- rownames(subset_data)
colnames(omics_matrix) <- colnames(subset_data)[1:3]  # Exclude Phenotype and Trimester

# Create a color palette for the heatmap
colors <- c("grey", "steelblue", "springgreen4", "mediumpurple")

# Define colors for Genotype, Lipidomics, and Methylation
category_colors <- c("Genotype" = "steelblue", "Lipidomics" = "mediumpurple", "Methylation" = "springgreen4")

# Create a vector of colors for the columns based on the omics category
col_colors <- sapply(colnames(omics_matrix), function(category) category_colors[category])

# Define colors for Phenotype and Trimester
phenotype_colors <- c("MISSING" = "grey", "CONTROL" = "blue", "CASE" = "red")
trimester_colors <- c("1" = "yellow", "2" = "orange", "3" = "purple")

# Create vectors of colors for Phenotype and Trimester columns
col_phenotype <- phenotype_colors[as.character(subset_data$Phenotype)]
col_trimester <- trimester_colors[as.character(subset_data$Trimester)]

#pdf("DNBC_Methylation_Mother.pdf", h=8, w=12)
tiff(filename="DNBC_Methylation_Mother.tiff", width=14, height=8, units="in",res=300)
# Create the rotated heatmap with ColSideColors
heatmap.2(t(omics_matrix),  # Transpose the matrix to rotate the plot
          Colv = NA, # To suppress column clustering
          Rowv = NA, # To suppress row clustering
          dendrogram = "none", # No dendrogram
          col = colors,
          cexRow=0.7,
          key = TRUE, key.title = "Omics Category",
          key.xlab = "Grey: Missing, Blue: Genotype, Purple: Lipidomics, Green: Methylation",
          main =  "Multi-Omics DNBC Overview ",
          trace = "none", # No trace
          density.info = "none", # No density plot
          ColSideColors = col_phenotype)  # Assign colors to Phenotype, Trimester, and omics categories

dev.off()


