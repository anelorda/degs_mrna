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
    output = TRUE,
    log = FALSE
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

# Filter for "mirtarbase" entries
tarbase_2_down <- mirna_stage2_down %>% filter(database == "tarbase")
tarbase_3_down <- mirna_stage3_down %>% filter(database == "tarbase")
tarbase_4_down <- mirna_stage4_down %>% filter(database == "tarbase")

mirtarbase_2_down$hgnc_symbol <- mirtarbase_2_down$target_symbol
mirtarbase_2_down$target_symbol <- NULL
mirtarbase_3_down$hgnc_symbol <- mirtarbase_3_down$target_symbol
mirtarbase_3_down$target_symbol <- NULL
mirtarbase_4_down$hgnc_symbol <- mirtarbase_4_down$target_symbol
mirtarbase_4_down$target_symbol <- NULL

# Rename target_symbol to hgnc_symbol and remove target_symbol column for targetscan
tarbase_2_down$hgnc_symbol <- tarbase_2_down$target_symbol
tarbase_2_down$target_symbol <- NULL
tarbase_3_down$hgnc_symbol <- tarbase_3_down$target_symbol
tarbase_3_down$target_symbol <- NULL
tarbase_4_down$hgnc_symbol <- tarbase_4_down$target_symbol
tarbase_4_down$target_symbol <- NULL

# Load downregulated gene lists
group2_downregulated <- read.table("deseq_mrna/downregulated_group2.csv", header = TRUE, sep = ",")
group3_downregulated <- read.table("deseq_mrna/downregulated_group3.csv", header = TRUE, sep = ",")
group4_downregulated <- read.table("deseq_mrna/downregulated_group4.csv", header = TRUE, sep = ",")

# Create the Venn diagram for the three stages with mirtarbase
create_venn_diagram(
  dataframes = list(mirtarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage 2",
  output_filename = "deseq_mrna/2mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage 3",
  output_filename = "deseq_mrna/3mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets found by mirtarbase and DEGs from mRNA-seq for stage 4",
  output_filename = "deseq_mrna/4mirnatargets_vs_DEGs.png"
)

# Create the Venn diagram for the three stages with tarbase
create_venn_diagram(
  dataframes = list(tarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (tarbase)","mRNA"),
  title = "miRNA targets found by tarbase and DEGs from mRNA-seq for stage 2",
  output_filename = "deseq_mrna/2tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets found by tarbase and DEGs from mRNA-seq for stage 3",
  output_filename = "deseq_mrna/3tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets found by tarbase and DEGs from mRNA-seq for stage 4",
  output_filename = "deseq_mrna/4tarbase_vs_DEGs.png"
)

 overlap2<- intersect(mirtarbase_2_down$hgnc_symbol, group2_downregulated$hgnc_symbol)
 overlap3<- intersect(mirtarbase_3_down$hgnc_symbol, group3_downregulated$hgnc_symbol)
 overlap4<- intersect(mirtarbase_4_down$hgnc_symbol, group4_downregulated$hgnc_symbol)
 


