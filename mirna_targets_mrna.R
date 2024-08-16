library(VennDiagram)
library(dplyr)


# Define the function to create the Venn diagram with title and save it
create_venn_diagram <- function(dataframes, labels, title, output_filename) {
  # Ensure the number of dataframes and labels match
  if (length(dataframes) != length(labels)) {
    stop("The number of data frames and labels must be the same.")
  }
  
  # Filter out NA values and extract hgnc_symbol columns
  symbol_lists <- lapply(dataframes, function(df) {
    df %>%
      filter(!is.na(hgnc_symbol)) %>%
      pull(hgnc_symbol)
  })
  
  # Create the Venn diagram
  venn.plot <- venn.diagram(
    x = symbol_lists,
    category.names = labels,
    main = title,
    filename = output_filename,  # Save to file
    output = TRUE
  )
  
  # Plot the Venn diagram if the filename is NULL
  if (is.null(output_filename)) {
    grid.draw(venn.plot)
  }
}

# Load your data
mirna_stage2_down <- read.table("Downloads/multimir_results_stage2_down.csv", header = TRUE, sep = ",")
mirna_stage3_down <- read.table("Downloads/multimir_results_stage3_down.csv", header = TRUE, sep = ",")
mirna_stage4_down <- read.table("Downloads/multimir_results_stage4_down.csv", header = TRUE, sep = ",")

# Filter for "mirtarbase" entries
mirtarbase_2_down <- mirna_stage2_down %>% filter(database == "mirtarbase")
mirtarbase_3_down <- mirna_stage3_down %>% filter(database == "mirtarbase")
mirtarbase_4_down <- mirna_stage4_down %>% filter(database == "mirtarbase")

mirtarbase_2_down$hgnc_symbol <- mirtarbase_2_down$target_symbol
mirtarbase_2_down$target_symbol <- NULL
mirtarbase_3_down$hgnc_symbol <- mirtarbase_3_down$target_symbol
mirtarbase_3_down$target_symbol <- NULL
mirtarbase_4_down$hgnc_symbol <- mirtarbase_4_down$target_symbol
mirtarbase_4_down$target_symbol <- NULL

# Create the Venn diagram for the three stages
create_venn_diagram(
  dataframes = list(mirtarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage2",
  output_filename = "deseq_mrna/2mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage3",
  output_filename = "deseq_mrna/3mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage4",
  output_filename = "deseq_mrna/4mirnatargets_vs_DEGs.png"
)
