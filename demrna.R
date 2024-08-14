library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)      
        
        
data <- read.table("Downloads/mRNA_for_miRNA.xlsx - sheet.csv", header = TRUE, sep = ",")

samples_to_remove <- c("X009T", "X010T", "X013T", "X014T", "X019T", "X021T")
count_data <- data[, !(colnames(data) %in% samples_to_remove)]

colnames(data)[1] <- "geneids"
rownames(count_data) <- count_data$X

count_data <- count_data[, -1]

count_data <- count_data[which(rowSums(count_data) > 5),]

stage_vector <- c('control', 'control','control', 'control', 'control', 'control', 
                  'control', 'control', 'control', 'control', 'control', 'group3', 
                  'group3', 'group3', 'group3', 'group3', 'group3', 'group3', 
                  'group4', 'group2', 'group3', 'group3','group2', 
                  'group3', 'group4', 'group3', 'group3')
row_names <- colnames(count_data)
sample_info <- data.frame(stage = stage_vector, row.names = row_names)
# Add pseudocount of 1 to count data
count_data_pseudocounts <- count_data + 1

dds <- DESeqDataSetFromMatrix(countData = count_data_pseudocounts, colData = sample_info, design = ~ stage)
sample_info$stage <- relevel(factor(sample_info$stage), ref = "control")
dds <- DESeq(dds)

# Extract results
res_group2_vs_control <- results(dds, contrast = c("stage", "group2", "control"))
res_group3_vs_control <- results(dds, contrast = c("stage", "group3", "control"))
res_group4_vs_control <- results(dds, contrast = c("stage", "group4", "control"))

# Filter results
res_group2_vs_control <- res_group2_vs_control[which(abs(res_group2_vs_control$log2FoldChange) > 1 & res_group2_vs_control$padj < 0.05), ]
res_group3_vs_control <- res_group3_vs_control[which(abs(res_group3_vs_control$log2FoldChange) > 1 & res_group3_vs_control$padj < 0.05), ]
res_group4_vs_control <- res_group4_vs_control[which(abs(res_group4_vs_control$log2FoldChange) > 1 & res_group4_vs_control$padj < 0.05), ]

# Convert to data frames
res_group2_vs_control_df <- as.data.frame(res_group2_vs_control)
res_group3_vs_control_df <- as.data.frame(res_group3_vs_control)
res_group4_vs_control_df <- as.data.frame(res_group4_vs_control)

# Add ensembl_gene_id for merging
res_group2_vs_control_df$ensembl_gene_id <- rownames(res_group2_vs_control_df)
res_group3_vs_control_df$ensembl_gene_id <- rownames(res_group3_vs_control_df)
res_group4_vs_control_df$ensembl_gene_id <- rownames(res_group4_vs_control_df)

# Retrieve gene mapping (Ensure you have the gene IDs for this)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = unique(c(rownames(res_group2_vs_control_df), rownames(res_group3_vs_control_df), rownames(res_group4_vs_control_df))), mart = ensembl)

# Merge with gene mapping
res_group2_vs_control_df <- merge(res_group2_vs_control_df, gene_mapping, by = "ensembl_gene_id", all.x = TRUE)
res_group3_vs_control_df <- merge(res_group3_vs_control_df, gene_mapping, by = "ensembl_gene_id", all.x = TRUE)
res_group4_vs_control_df <- merge(res_group4_vs_control_df, gene_mapping, by = "ensembl_gene_id", all.x = TRUE)

# Rename columns of res_group4_vs_control_df to add suffix "_group4"
colnames(res_group4_vs_control_df) <- paste0(colnames(res_group4_vs_control_df), "_group4")
colnames(res_group4_vs_control_df)[colnames(res_group4_vs_control_df) == "ensembl_gene_id_group4"] <- "ensembl_gene_id"

# Merge the data frames
combined_results <- merge(res_group2_vs_control_df, res_group3_vs_control_df, by = "ensembl_gene_id", suffixes = c("_group2", "_group3"))
combined_results <- merge(combined_results, res_group4_vs_control_df, by = "ensembl_gene_id")

# Check column names
colnames(combined_results)

# Calculate a combined score for ranking
combined_results$score <- with(combined_results, {

# Ensure there are no missing values in the necessary columns
padj_group2[is.na(padj_group2)] <- 1
padj_group3[is.na(padj_group3)] <- 1
padj_group4[is.na(padj_group4)] <- 1
  
  log_padj <- -log10(pmin(padj_group2, padj_group3, padj_group4))
  fold_change <- abs(pmin(log2FoldChange_group2, log2FoldChange_group3, log2FoldChange_group4))
  score <- log_padj * fold_change
  score
})

# Check if the score calculation worked
head(combined_results$score)

# Select top 100 DEGs based on the combined score
top100_combined <- combined_results %>%
  arrange(desc(score)) %>%
  head(100)

#compare with Table S6
top100_s6 <- read.table("Downloads/Table1_Transcriptome profiling and analysis of patients with esophageal squamous cell carcinoma from Kazakhstan.XLSX - Table S6.csv", header = TRUE,
                        sep = ",")


# Save top 100 DEGs to CSV file
write.csv(top100_combined, file = "top100_combined_DEGs.csv", row.names = FALSE)

# Save filtered results to CSV files
write.csv(res_group2_vs_control_df, file = "filtered_results_group2_vs_control.csv", row.names = FALSE)
write.csv(res_group3_vs_control_df, file = "filtered_results_group3_vs_control.csv", row.names = FALSE)
write.csv(res_group4_vs_control_df, file = "filtered_results_group4_vs_control.csv", row.names = FALSE)

# Plot MA plots
plotMA(res_group2_vs_control, ylim = c(-5, 5))
plotMA(res_group3_vs_control, ylim = c(-5, 5))
plotMA(res_group4_vs_control, ylim = c(-5, 5))

# Plot Volcano plots
EnhancedVolcano(res_group2_vs_control, lab = rownames(res_group2_vs_control), x = 'log2FoldChange', y = 'pvalue')
EnhancedVolcano(res_group3_vs_control, lab = rownames(res_group3_vs_control), x = 'log2FoldChange', y = 'pvalue')
EnhancedVolcano(res_group4_vs_control, lab = rownames(res_group4_vs_control), x = 'log2FoldChange', y = 'pvalue')