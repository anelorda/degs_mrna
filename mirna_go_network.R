library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(ggraph)

gene_list <- unique(oncomir_genes$gene)
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Initialize a list to store GO terms for each miRNA
miRNA_GO_list <- list()

# Loop through each unique miRNA
for (mirna in unique(oncomir_genes$miRNA)) {
  # Get the target genes for this miRNA
  target_genes <- oncomir_genes$gene[oncomir_genes$miRNA == mirna]
  
  # Convert target genes to Entrez IDs
  target_entrez <- entrez_ids$ENTREZID[entrez_ids$SYMBOL %in% target_genes]
  
  # Perform GO enrichment analysis
  go_results <- enrichGO(gene = target_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
  
  # Store the GO terms in the list
  miRNA_GO_list[[mirna]] <- go_results
}


# Create an empty data frame for the edge list
edge_list <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

# Loop through each miRNA and extract the GO terms
for (mirna in names(miRNA_GO_list)) {
  go_terms <- miRNA_GO_list[[mirna]]@result$Description  # Extract GO terms
  
  # Add the miRNA-GO pairs to the edge list
  edge_list <- rbind(edge_list, data.frame(from = rep(mirna, length(go_terms)), to = go_terms))
}


# Create the graph object
graph <- graph_from_data_frame(edge_list, directed = FALSE)


# Plot the network using ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(color = ifelse(V(graph)$name %in% oncomir_genes$miRNA, "miRNA", "GO Term"))) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()



library(pheatmap)

# Function to get top N GO terms for each miRNA and their p-values
get_top_go_terms_with_pvalue <- function(mirna_go_list, n = 10) {
  top_terms <- lapply(names(mirna_go_list), function(mirna) {
    mirna_go_list[[mirna]]@result %>%
      arrange(p.adjust) %>%
      head(n) %>%
      select(Description, p.adjust) %>%
      mutate(miRNA = mirna)
  })
  return(bind_rows(top_terms))
}

# Get top GO terms with p-values
top_go_data <- get_top_go_terms_with_pvalue(miRNA_GO_list, n = 10)  # Adjust n as needed

# Create a matrix for the heatmap
heatmap_matrix <- top_go_data %>%
  pivot_wider(names_from = miRNA, values_from = p.adjust, values_fill = 1) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# Log-transform p-values for better color scaling
heatmap_matrix <- -log10(heatmap_matrix)

# Create the heatmap
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("white", "red"))(50),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         filename = "go_terms_heatmap.png",
         width = 15,
         height = 20)

library(ggraph)
library(ggplot2)
library(igraph)

plot <- ggraph(graph, layout = "kk") +
  geom_edge_link(alpha = 0.1, color = "grey80") +  # Lighter edge color
  geom_node_point(aes(color = ifelse(V(graph)$name %in% oncomir_genes$miRNA, "miRNA", "GO Term")), 
                  size = 1.5, alpha = 0.7) +  # Smaller, slightly transparent nodes
  geom_node_text(aes(label = name), repel = TRUE, size = 2, 
                 max.overlaps = 10, segment.size = 0.2) +  # Adjust text and connecting lines
  scale_color_manual(values = c("miRNA" = "red", "GO Term" = "blue")) +  # Custom colors
  theme_void() +  # Remove default theme elements
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # White background
    panel.background = element_rect(fill = "white", color = NA),  # White panel
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = NA),  # White legend background
    legend.key = element_rect(fill = "white", color = NA)  # White legend key background
  ) +
  labs(color = "Node Type") +  # Legend title
  guides(color = guide_legend(override.aes = list(size = 3)))  # Larger legend points

# Save the plot
ggsave("network_plot_improved.png", plot, width = 30, height = 30, units = "in", dpi = 300, bg = "white")
