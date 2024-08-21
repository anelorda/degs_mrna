library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)      
library(VennDiagram)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)

        
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
# Assign rownames to ensembl_gene_id
res_group2_vs_control_df$ensembl_gene_id <- rownames(res_group2_vs_control_df)
res_group3_vs_control_df$ensembl_gene_id <- rownames(res_group3_vs_control_df)
res_group4_vs_control_df$ensembl_gene_id <- rownames(res_group4_vs_control_df)


# Retrieve gene mapping
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(c(rownames(res_group2_vs_control_df), rownames(res_group3_vs_control_df), rownames(res_group4_vs_control_df)))
gene_mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = gene_ids, mart = ensembl)

# Function to merge gene mapping and replace empty hgnc_symbol with NA
merge_and_clean <- function(df) {
  df <- merge(df, gene_mapping, by = "ensembl_gene_id", all.x = TRUE)
  df$hgnc_symbol[df$hgnc_symbol == ""] <- NA
  return(df)
}

# Apply function to each group
res_group2_vs_control_df <- merge_and_clean(res_group2_vs_control_df)
res_group3_vs_control_df <- merge_and_clean(res_group3_vs_control_df)
res_group4_vs_control_df <- merge_and_clean(res_group4_vs_control_df)


# Save data frames
write.csv(res_group2_vs_control_df, "deseq_mrna/res_group2_vs_control.csv", row.names = FALSE)
write.csv(res_group3_vs_control_df, "deseq_mrna/res_group3_vs_control.csv", row.names = FALSE)
write.csv(res_group4_vs_control_df, "deseq_mrna/res_group4_vs_control.csv", row.names = FALSE)

' Assign rownames to ensembl_gene_id
res_group2_vs_control_df$ensembl_gene_id <- rownames(res_group2_vs_control_df)
res_group3_vs_control_df$ensembl_gene_id <- rownames(res_group3_vs_control_df)
res_group4_vs_control_df$ensembl_gene_id <- rownames(res_group4_vs_control_df)
'

# Function to create and save DEGs
save_degs <- function(df, group_name) {
  upregulated <- df[df$log2FoldChange > 1,]
  downregulated <- df[df$log2FoldChange < -1,]
  
  # Save full upregulated and downregulated data frames
  write.csv(upregulated, paste0("deseq_mrna/upregulated_", group_name, ".csv"), row.names = FALSE)
  write.csv(downregulated, paste0("deseq_mrna/downregulated_", group_name, ".csv"), row.names = FALSE)
  
  # Save top 100 upregulated and downregulated data frames
  top100_upregulated <- head(upregulated[order(-upregulated$log2FoldChange), ], 100) %>%
    filter(!is.na(hgnc_symbol) & hgnc_symbol != "") 
  top100_downregulated <- head(downregulated[order(downregulated$log2FoldChange), ], 100) %>%
    filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
  write.csv(top100_upregulated, paste0("deseq_mrna/top100_upregulated_", group_name, ".csv"), row.names = FALSE)
  write.csv(top100_downregulated, paste0("deseq_mrna/top100_downregulated_", group_name, ".csv"), row.names = FALSE)
}

# Apply function to each group
save_degs(res_group2_vs_control_df, "group2")
save_degs(res_group3_vs_control_df, "group3")
save_degs(res_group4_vs_control_df, "group4")

# Filter and extract top 100 DEGs for each group
top100_group2 <- res_group2_vs_control_df %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(100)

top100_group3 <- res_group3_vs_control_df %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(100)

top100_group4 <- res_group4_vs_control_df %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(100)

# Extract HGNC symbols from top 100 DEGs for each group
hgnc_group2 <- top100_group2$hgnc_symbol
hgnc_group3 <- top100_group3$hgnc_symbol
hgnc_group4 <- top100_group4$hgnc_symbol

# Create a list of HGNC symbols
venn_list <- list(
  Group2 = hgnc_group2,
  Group3 = hgnc_group3,
  Group4 = hgnc_group4
)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Group2", "Group3", "Group4"),
  filename = NULL,
  output = TRUE
)

# Save to file
pdf("VennDiagram_HGNC_Symbols.pdf")
grid.draw(venn.plot)
dev.off()

# Display Venn diagram in R
grid.draw(venn.plot)


#Save top 100 for each group as csv
write.csv(top100_group2, "deseq_mrna/top100_group2_FC.csv", row.names = FALSE)
write.csv(top100_group3, "deseq_mrna/top100_group3_FC.csv", row.names = FALSE)
write.csv(top100_group4, "deseq_mrna/top100_group4_FC.csv", row.names = FALSE)



# Rename columns of res_group4_vs_control_df to add suffix "_group4"
colnames(res_group4_vs_control_df) <- paste0(colnames(res_group4_vs_control_df), "_group4")
colnames(res_group4_vs_control_df)[colnames(res_group4_vs_control_df) == "ensembl_gene_id_group4"] <- "ensembl_gene_id"

combined_results <- merge(res_group2_vs_control_df, res_group3_vs_control_df, by = "ensembl_gene_id", suffixes = c("_group2", "_group3"))
combined_results <- merge(combined_results, res_group4_vs_control_df, by = "ensembl_gene_id")

# Filter out rows where any of the hgnc_symbol columns contain NA
combined_results_filtered <- combined_results %>%
  filter(!is.na(hgnc_symbol))

# Select the maximum absolute log2FoldChange across the three groups for sorting
combined_results_filtered$max_abs_logFC <- with(combined_results_filtered, pmax(abs(log2FoldChange_group2), 
                                                                                abs(log2FoldChange_group3), 
                                                                                abs(log2FoldChange_group4)))

# Sort by the maximum absolute log2FoldChange
top100_combined <- combined_results_filtered %>%
  arrange(desc(max_abs_logFC)) %>%
  head(100)

'# Select top 100 DEGs based on the combined score
top100_combined <- combined_results %>%
  arrange(desc(score)) %>%
  head(100)
'

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


# Perform variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = stage)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# Save PCA plot
ggsave("PCA_plot.png")


#heatmap
top100_group2_counts <- count_data_pseudocounts[rownames(count_data_pseudocounts) %in% top100_group2$ensembl_gene_id, ]
top100_group3_counts <- count_data_pseudocounts[rownames(count_data_pseudocounts) %in% top100_group3$ensembl_gene_id, ]
top100_group4_counts <- count_data_pseudocounts[rownames(count_data_pseudocounts) %in% top100_group4$ensembl_gene_id, ]

top100_group2_counts$ensembl_gene_id <- rownames(top100_group2_counts)
top100_group2_counts <- merge_and_clean(top100_group2_counts)
rownames(top100_group2_counts) <- top100_group2_counts$hgnc_symbol
top100_group2_counts$ensembl_gene_id <- NULL
top100_group2_counts$hgnc_symbol <- NULL

top100_group3_counts$ensembl_gene_id <- rownames(top100_group3_counts)
top100_group3_counts <- merge_and_clean(top100_group3_counts)
rownames(top100_group3_counts) <- top100_group3_counts$hgnc_symbol
top100_group3_counts$ensembl_gene_id <- NULL
top100_group3_counts$hgnc_symbol <- NULL

top100_group4_counts$ensembl_gene_id <- rownames(top100_group4_counts)
top100_group4_counts <- merge_and_clean(top100_group4_counts)
rownames(top100_group4_counts) <- top100_group4_counts$hgnc_symbol
top100_group4_counts$ensembl_gene_id <- NULL
top100_group4_counts$hgnc_symbol <- NULL

# Assuming top100_group2, top100_group3, and top100_group4 are data frames containing 'hgnc_symbol' column
overlap_groups <- Reduce(intersect, list(top100_group2$hgnc_symbol, top100_group3$hgnc_symbol, top100_group4$hgnc_symbol))
overlap_groups


pheatmap(log2(top100_group2_counts), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         annotation_col = sample_info, 
         main = "Top 100 DEGs for Group 2 vs Control")

pheatmap(log2(top100_group3_counts), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         annotation_col = sample_info, 
         main = "Top 100 DEGs for Group 3 vs Control")

pheatmap(log2(top100_group4_counts), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         annotation_col = sample_info, 
         main = "Top 100 DEGs for Group 4 vs Control")



#GO
run_go_analysis <- function(results_df, title, ont) {
  diff_genes <- results_df$hgnc_symbol[results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1]
  
  ego <- enrichGO(gene = diff_genes,
                  OrgDb = org.Hs.eg.db,  # Human database
                  keyType = "SYMBOL",
                  ont = ont,  # Ontology type (BP, MF, CC)
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    print(dotplot(ego, showCategory = 20, title = paste("GO Enrichment Analysis -", title, "-", ont)))
  } else {
    print(paste("No significant GO terms found for", title, "-", ont))
  }
}

# Run GO analysis for Biological Process (BP)
run_go_analysis(group2_upregulated, "Non-stimulated", "BP")
run_go_analysis(top100_group2_upregulated, "Stimulated", "BP")
run_go_analysis(top100_group2_downregulated, "Stimulated Rest", "BP")
