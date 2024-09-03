library(dplyr)
library(mixOmics)
library(caret)



mirna_norm_counts <- read.table("deseq_mrna/diablo1.csv", header = TRUE, sep = ",")
normalized_counts <- read.table("deseq_mrna/normalised_counts_mRNA.csv", header = TRUE, sep = ",")
new_rownames <- c(
'control1', 'control2', 'X001T', 'X002T', 'X003T',
'X004T', 'X005T', 'X007T', 'X008T', 'X011T',
'X012T', 'X015T', 'X016T', 'X017T', 'X023T',
'X025T', 'X029T', 'X033T'
)

rownames(mirna_norm_counts) <- new_rownames
mirna_norm_counts$X <- NULL

rownames(normalized_counts) <- normalized_counts$X
normalized_counts$X <- NULL

#combined rownames for significant genes in all thre groups
significant_norm_counts <- normalized_counts[rownames(normalized_counts) %in% combined_rownames, ]
View(significant_norm_counts)
trans_sign_norm_counts <- t(significant_norm_counts)

trans_sign_norm_counts <- trans_sign_norm_counts[-1,]
trans_sign_norm_counts <- trans_sign_norm_counts[-5,]
trans_sign_norm_counts <- trans_sign_norm_counts[-3,]
trans_sign_norm_counts <- trans_sign_norm_counts[-2,]
trans_sign_norm_counts <- trans_sign_norm_counts[-7,]
trans_sign_norm_counts <- trans_sign_norm_counts[-1,]
trans_sign_norm_counts <- trans_sign_norm_counts[-3,]
trans_sign_norm_counts <- trans_sign_norm_counts[-3,]
trans_sign_norm_counts <- trans_sign_norm_counts[-3,]
rownames(trans_sign_norm_counts)[1] <- "control1"
rownames(trans_sign_norm_counts)[2] <- "control2"
colnames(mirna_norm_counts) <- gsub("\\.", "-", colnames(mirna_norm_counts))
View(mirna_norm_counts)


# Identify features with near-zero variance in the mRNA block
zero_var_features <- nearZeroVar(trans_sign_norm_counts)

# Remove the identified features
trans_sign_norm_counts_filtered <- trans_sign_norm_counts[, -zero_var_features]

# Update the data list
data_list <- list(mRNA = trans_sign_norm_counts_filtered, miRNA = mirna_norm_counts)


# Define the design matrix (for example, full integration)
design <- matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
diag(design) <- 0

factor_class_labels <- c('control1', 'control2', 'group3', 
                  'group3', 'group3', 'group3', 'group3', 'group3', 'group3', 
                  'group4', 'group2', 'group3', 'group3','group2', 
                  'group3', 'group4', 'group3', 'group3')
ncomp <- 2
diablo_result <- block.splsda(X = data_list, Y = factor_class_labels, design = design, ncomp = ncomp)
diablo_result <- wrapper.sgccda(X = data_list,
                                Y = factor_class_labels,
                                design = design,
                                keepX = list(mRNA = c(10, 10), miRNA = c(15, 15)),
                                ncomp = 2,
                                scheme = "horst")

# Circos Plot for Correlations
corMat <- circosPlot(diablo_result, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)


#Plot the sample plot
plotIndiv(diablo_result, legend = TRUE, title = 'DIABLO Sample Plot')

# Plot the correlation circle plot
plotVar(diablo_result, legend = TRUE, title = 'DIABLO Correlation Circle Plot')

# Reset par to default
par(mfrow = c(1, 1))
par(mar = c(5, 4, 4, 2) + 0.1)  # Default margin settings
# Adjust margins: bottom, left, top, right
par(mar = c(4, 4, 2, 1)) # This reduces the margin size

# Plot the network again
network(diablo_result, comp = list(mRNA = 1, miRNA = 1), cutoff = 0.7)
# Save the plot to a file
png("diablo_network_plot.png", width = 1500, height = 1200, res = 150)
par(mar = c(4, 4, 2, 1))
network(diablo_result, comp = list(mRNA = 1, miRNA = 1), cutoff = 0.7)
dev.off()

# Try plotting again
network(diablo_result, comp = list(mRNA = 1, miRNA = 1), cutoff = 1)
# Extract loadings for the first component
loadings_mrna <- diablo_result$loadings$mRNA[, 1]
loadings_mirna <- diablo_result$loadings$miRNA[, 1]

# Extract the latent components for each block (mRNA and miRNA)
latent_components_mrna <- diablo_result$variates$mRNA[, 1]
latent_components_mirna <- diablo_result$variates$miRNA[, 1]

correlation_matrix <- cor(latent_components_mrna, latent_components_mirna)

# Set a correlation threshold
threshold <- 0.7

# Filter the correlation matrix based on the threshold
high_correlations <- correlation_matrix[abs(correlation_matrix) > threshold]

# Save the filtered correlations to a CSV or text file
write.csv(high_correlations, "diablo_high_correlations.csv")
write.table(high_correlations, "diablo_high_correlations.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Combine the loadings into a data frame
loadings_df <- data.frame(mRNA = loadings_mrna, miRNA = loadings_mirna)


# Add a title to the plot
title(main = "DIABLO Network")
