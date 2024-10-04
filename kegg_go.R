library(VennDiagram)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)


stage2_mirna_targets_up <- intersect(mirtarbase_2_up$hgnc_symbol, tarbase_2_up$hgnc_symbol)
stage3_mirna_targets_up <- intersect(mirtarbase_3_up$hgnc_symbol, tarbase_3_up$hgnc_symbol)
stage4_mirna_targets_up <- intersect(mirtarbase_4_up$hgnc_symbol, tarbase_4_up$hgnc_symbol)
stage2_mirna_targets_down <- intersect(mirtarbase_2_down$hgnc_symbol, tarbase_2_down$hgnc_symbol)
stage3_mirna_targets_down <- intersect(mirtarbase_3_down$hgnc_symbol, tarbase_3_down$hgnc_symbol)
stage4_mirna_targets_down <- intersect(mirtarbase_4_down$hgnc_symbol, tarbase_4_down$hgnc_symbol)



# Create the Venn diagram function
create_three_group_venn <- function(stage2, stage3, stage4, title, output_filename) {
  
  # Define the list of groups
  venn_groups <- list(
    "Stage 2" = stage2,
    "Stage 3" = stage3,
    "Stage 4" = stage4
  )
  
  # Define colors for the circles
  fill_colors <- c("#999999", "#E69F00", "#56B4E9")
  
  # Create the Venn diagram
  venn.plot <- venn.diagram(
    x = venn_groups,
    category.names = names(venn_groups),
    main = title,
    
    # Customizations for Venn diagram
    main.cex = 2,                          # Increase title size
    main.fontface = "bold",
    lwd = 2,                               # Line width
    lty = 'blank',                         # Line type
    fill = fill_colors,                    # Fill colors for the circles
    alpha = 0.5,                           # Transparency level
    cex = 1.5,                             # Font size for numbers
    fontface = "italic",                   # Font style for numbers
    
    # Customizing label appearance
    cat.cex = 1.7,                         # Size of the category labels
    cat.fontface = "bold",                 # Font face for category labels
    cat.col = fill_colors,                 # Colors for the category labels
    cat.default.pos = "outer",             # Place labels outside the circles
    cat.dist = c(0.05, -0.55, 0.05),           # Adjust distance from circles
    cat.pos = c(-25, 25, 180),             # Adjust label positions
    
    height = 3000,                         # Increase diagram height
    width = 3000,                          # Increase diagram width
    resolution = 300,                      # Increase resolution
    filename = output_filename,            # Save to file
    output = TRUE,
    log = FALSE
  )
  
  # Plot the Venn diagram if the filename is NULL
  if (is.null(output_filename)) {
    grid.draw(venn.plot)
  }
}


# Example usage
create_three_group_venn(
  stage2 = stage2_mirna_targets_up,
  stage3 = stage3_mirna_targets_up,
  stage4 = stage4_mirna_targets_up,
  title = "Upregulated miRNA Targets Across Stages 2, 3, and 4",
  output_filename = "deseq_mrna//overlap_mirna_targets_up.png"
)

create_three_group_venn(
  stage2 = stage2_mirna_targets_down,
  stage3 = stage3_mirna_targets_down,
  stage4 = stage4_mirna_targets_down,
  title = "Downregulated miRNA Targets Across Stages 2, 3, and 4",
  output_filename = "deseq_mrna//overlap_mirna_targets_down.png"
)


targets_up <- intersect(stage4_mirna_targets_up, intersect(stage2_mirna_targets_up, stage3_mirna_targets_up))
targets_down <- intersect(stage2_mirna_targets_down, intersect(stage3_mirna_targets_down, stage4_mirna_targets_down))


targets <- c(targets_down, targets_up)
significant_genes <- read.table("deseq_mrna/significant_genes_across_samples.csv", header = TRUE, sep = ",")

filtered_df_targets <- significant_genes[significant_genes$hgnc_symbol %in% targets, ]

oncomir_genes <- read.table("merge_metatastatic/merge_oncomir_up.csv", header = TRUE, sep = ",")
oncomir_genes_list <- oncomir_genes$gene

run_go_analysis <- function(gene_symbols) {
  # Convert HGNC symbols to Entrez IDs
  entrez_ids <- bitr(gene_symbols, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  
  # Handle multiple mappings by keeping unique Entrez IDs
  unique_entrez_ids <- entrez_ids[!duplicated(entrez_ids$SYMBOL), ]
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene         = unique_entrez_ids$ENTREZID,
                            OrgDb        = org.Hs.eg.db,
                            keyType      = "ENTREZID",
                            ont          = "ALL",    # "BP", "CC", "MF", or "ALL"
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
  
  # Return the GO enrichment results
  return(go_enrichment)
}

# Assuming 'run_go_analysis' accepts a vector of gene symbols
go_results <- run_go_analysis(oncomir_genes_list)

# Dotplot for GO enrichment results
dotplot(go_results) + ggtitle("GO Enrichment Analysis - Dotplot")

# Barplot for GO enrichment results
barplot(go_results, showCategory=20) + ggtitle("GO Enrichment Analysis - Barplot")

# If you want to plot an enrichment map (useful for visualizing relationships)
emapplot(go_results) + ggtitle("GO Enrichment Analysis - Enrichment Map")

# If you want a network plot (requires significant terms with p.adjust < 0.05)
cnetplot(go_results, showCategory = 10) + ggtitle("GO Enrichment Analysis - Network Plot")
# View the GO analysis results
head(go_results)

# Barplot for GO enrichment results with increased text size
barplot(go_results, showCategory = 20) + 
  ggtitle("GO Enrichment Analysis - Barplot") + 
  theme(text = element_text(size = 14))  # Adjust the size as needed

edge_list <- data.frame(
  from = rep(oncomir_genes$miRNA, each = nrow(go_results)),
  to = go_results$Description
)

# Create a graph object
graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Plot the network
plot(graph, vertex.label = V(graph)$name, vertex.size = 5, edge.arrow.size = 0.5)

ggraph(graph, layout = 'fr') +
  geom_edge_link() +
  geom_node_point(aes(color = ifelse(V(graph)$name %in% oncomir_genes$miRNA, "miRNA", "GO Term"))) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

# Extend the function to include KEGG pathway analysis
run_kegg_analysis <- function(gene_symbols) {
  # Convert HGNC symbols to Entrez IDs
  entrez_ids <- bitr(gene_symbols, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  
  # Handle multiple mappings by keeping unique Entrez IDs
  unique_entrez_ids <- entrez_ids[!duplicated(entrez_ids$SYMBOL), ]
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene         = unique_entrez_ids$ENTREZID,
                                organism     = 'hsa',    # 'hsa' for human
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
  
  # Return the KEGG enrichment results
  return(kegg_enrichment)
}

# Run KEGG analysis
kegg_results <- run_kegg_analysis(oncomir_genes_list)

# Visualize the KEGG pathway enrichment results
# Dotplot for KEGG enrichment results
dotplot(kegg_results) + ggtitle("KEGG Pathway Enrichment - Dotplot")

# Barplot for KEGG enrichment results
barplot(kegg_results, showCategory=20) + ggtitle("KEGG Pathway Enrichment - Barplot")

kegg_results_with_similarity <- pairwise_termsim(kegg_results)

# Now you can plot the enrichment map
emapplot(kegg_results_with_similarity) + ggtitle("KEGG Pathway Enrichment - Enrichment Map")

# Barplot for KEGG enrichment results with increased text size
barplot(kegg_results, showCategory = 20) + 
  ggtitle("KEGG Pathway Enrichment - Barplot") + 
  theme(text = element_text(size = 14))  # Adjust the size as needed

# Network Plot for KEGG enrichment results
cnetplot(kegg_results, showCategory = 10) + ggtitle("KEGG Pathway Enrichment - Network Plot")



# Ensure oncomir_genes_list contains unique genes
oncomir_genes_list <- unique(oncomir_genes_list)

# Print some information about the data
cat("Number of unique genes:", length(oncomir_genes_list), "\n")
cat("First few genes:", head(oncomir_genes_list), "\n\n")

# Function to perform KEGG enrichment analysis
run_kegg_analysis <- function(gene_symbols, p_cutoff = 0.05, q_cutoff = 0.2) {
  # Convert HGNC symbols to Entrez IDs
  entrez_ids <- bitr(gene_symbols, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  
  cat("Number of genes mapped to Entrez IDs:", nrow(entrez_ids), "\n")
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                                organism     = 'hsa',
                                pvalueCutoff = p_cutoff,
                                qvalueCutoff = q_cutoff)
  
  return(kegg_enrichment)
}

# Run KEGG analysis
kegg_results <- run_kegg_analysis(filtered_df_targets$hgnc_symbol)

# Print summary of results
cat("\nKEGG Enrichment Results Summary:\n")
print(summary(kegg_results))

# If no results, try with more lenient cutoffs
if (nrow(as.data.frame(kegg_results)) == 0) {
  cat("\nNo enriched pathways found. Trying with more lenient cutoffs...\n")
  kegg_results <- run_kegg_analysis(oncomir_genes_list, p_cutoff = 0.1, q_cutoff = 0.3)
  cat("\nKEGG Enrichment Results Summary (lenient cutoffs):\n")
  print(summary(kegg_results))
}

# Visualize results if any are found
if (nrow(as.data.frame(kegg_results)) > 0) {
  # Dotplot
  png("kegg_dotplot.png", width = 800, height = 600)
  print(dotplot(kegg_results, showCategory = 20) + 
          ggtitle("KEGG Pathway Enrichment - Dotplot"))
  dev.off()
  
  # Barplot
  png("kegg_barplot.png", width = 800, height = 600)
  print(barplot(kegg_results, showCategory = 20) + 
          ggtitle("KEGG Pathway Enrichment - Barplot"))
  dev.off()
  
  cat("\nPlots saved as 'kegg_dotplot.png' and 'kegg_barplot.png'\n")
} else {
  cat("\nNo enriched pathways found even with lenient cutoffs.\n")
}




