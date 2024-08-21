library(VennDiagram)
library(dplyr)
library(grid)

# Define the function to create the Venn diagram with title and save it
create_venn_diagram <- function(dataframes, labels, title, output_filename) {
  if (length(dataframes) != length(labels)) {
    stop("The number of data frames and labels must be the same.")
  }
  
  symbol_lists <- lapply(dataframes, function(df) {
    df %>%
      filter(!is.na(hgnc_symbol)) %>%
      pull(hgnc_symbol)
  })
  
  venn.plot <- venn.diagram(
    x = symbol_lists,
    category.names = labels,
    main = title,
    main.cex = 2,
    main.pos = c(0.5, 1.1),
    cat.cex = 1.5,
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05),
    height = 3000,
    width = 3000,
    resolution = 300,
    cex = 1.5,
    fill = c("skyblue", "pink"),
    alpha = 0.5,
    filename = output_filename,
    output = TRUE,
    log = FALSE
  )
  
  if (is.null(output_filename)) {
    grid.draw(venn.plot)
  }
}

# Load your data
mirna_stage2_down <- read.table("Downloads/multimir_results_stage2_down.csv", header = TRUE, sep = ",")
mirna_stage3_down <- read.table("Downloads/multimir_results_stage3_down.csv", header = TRUE, sep = ",")
mirna_stage4_down <- read.table("Downloads/multimir_results_stage4_down.csv", header = TRUE, sep = ",")
mirna_stage2_up <- read.table("Downloads/multimir_results_stage2_up.csv", header = TRUE, sep = ",")
mirna_stage3_up <- read.table("Downloads/multimir_results_stage3_up.csv", header = TRUE, sep = ",")
mirna_stage4_up <- read.table("Downloads/multimir_results_stage4_up.csv", header = TRUE, sep = ",")

# Function to filter and rename for tarbase and mirtarbase
filter_and_rename <- function(data, database_name) {
  data_filtered <- data %>% filter(database == database_name)
  data_filtered$hgnc_symbol <- data_filtered$target_symbol
  data_filtered$target_symbol <- NULL
  return(data_filtered)
}

# Filter and rename columns for tarbase and mirtarbase
mirtarbase_2_down <- filter_and_rename(mirna_stage2_down, "mirtarbase")
mirtarbase_3_down <- filter_and_rename(mirna_stage3_down, "mirtarbase")
mirtarbase_4_down <- filter_and_rename(mirna_stage4_down, "mirtarbase")
mirtarbase_2_up <- filter_and_rename(mirna_stage2_up, "mirtarbase")
mirtarbase_3_up <- filter_and_rename(mirna_stage3_up, "mirtarbase")
mirtarbase_4_up <- filter_and_rename(mirna_stage4_up, "mirtarbase")

tarbase_2_down <- filter_and_rename(mirna_stage2_down, "tarbase")
tarbase_3_down <- filter_and_rename(mirna_stage3_down, "tarbase")
tarbase_4_down <- filter_and_rename(mirna_stage4_down, "tarbase")
tarbase_2_up <- filter_and_rename(mirna_stage2_up, "tarbase")
tarbase_3_up <- filter_and_rename(mirna_stage3_up, "tarbase")
tarbase_4_up <- filter_and_rename(mirna_stage4_up, "tarbase")

# Load gene lists
group2_downregulated <- read.table("deseq_mrna/downregulated_group2.csv", header = TRUE, sep = ",")
group3_downregulated <- read.table("deseq_mrna/downregulated_group3.csv", header = TRUE, sep = ",")
group4_downregulated <- read.table("deseq_mrna/downregulated_group4.csv", header = TRUE, sep = ",")
group2_upregulated <- read.table("deseq_mrna/upregulated_group2.csv", header = TRUE, sep = ",")
group3_upregulated <- read.table("deseq_mrna/upregulated_group3.csv", header = TRUE, sep = ",")
group4_upregulated <- read.table("deseq_mrna/upregulated_group4.csv", header = TRUE, sep = ",")

# Create Venn diagrams for mirtarbase DOWNREGULATED
create_venn_diagram(list(mirtarbase_2_down, group2_downregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 2, downregulated", "deseq_mrna/2_down_mirnatargets_vs_DEGs.png")
create_venn_diagram(list(mirtarbase_3_down, group3_downregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 3, downregulated", "deseq_mrna/3_down_mirnatargets_vs_DEGs.png")
create_venn_diagram(list(mirtarbase_4_down, group4_downregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 4, downregulated", "deseq_mrna/4_down_mirnatargets_vs_DEGs.png")

# Create Venn diagrams for tarbase DOWNREGULATED
create_venn_diagram(list(tarbase_2_down, group2_downregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 2, downregulated", "deseq_mrna/2_down_tarbase_vs_DEGs.png")
create_venn_diagram(list(tarbase_3_down, group3_downregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 3, downregulated", "deseq_mrna/3_down_tarbase_vs_DEGs.png")
create_venn_diagram(list(tarbase_4_down, group4_downregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 4, downregulated", "deseq_mrna/4_down_tarbase_vs_DEGs.png")

# Create Venn diagrams for mirtarbase UPREGULATED
create_venn_diagram(list(mirtarbase_2_up, group2_upregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 2, upregulated", "deseq_mrna/2_up_mirnatargets_vs_DEGs.png")
create_venn_diagram(list(mirtarbase_3_up, group3_upregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 3, upregulated", "deseq_mrna/3_up_mirnatargets_vs_DEGs.png")
create_venn_diagram(list(mirtarbase_4_up, group4_upregulated), c("miRNA targets (mirtarbase)", "mRNA"), "miRNA targets (mirtarbase) and DEGs for stage 4, upregulated", "deseq_mrna/4_up_mirnatargets_vs_DEGs.png")

# Create Venn diagrams for tarbase UPREGULATED
create_venn_diagram(list(tarbase_2_up, group2_upregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 2, upregulated", "deseq_mrna/2_up_tarbase_vs_DEGs.png")
create_venn_diagram(list(tarbase_3_up, group3_upregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 3, upregulated", "deseq_mrna/3_up_tarbase_vs_DEGs.png")
create_venn_diagram(list(tarbase_4_up, group4_upregulated), c("miRNA targets (tarbase)", "mRNA"), "miRNA targets (tarbase) and DEGs for stage 4, upregulated", "deseq_mrna/4_up_tarbase_vs_DEGs.png")

# Define the function to create the Venn diagram with four components
create_venn_diagram_four <- function(dataframes, labels, title, output_filename) {
  # Ensure the number of dataframes and labels match
  if (length(dataframes) != 4 || length(labels) != 4) {
    stop("This function requires exactly four dataframes and four labels.")
  }
  
  # Filter out NA values and extract hgnc_symbol columns
  symbol_lists <- lapply(dataframes, function(df) {
    df %>%
      filter(!is.na(hgnc_symbol)) %>%
      pull(hgnc_symbol)
  })
  
  # Create the Venn diagram with four components
  venn.plot <- venn.diagram(
    x = symbol_lists,
    category.names = labels,
    main = title,
    
    # Adjustments to ensure text isn't cut off
    main.cex = 2,                          # Increase title size
    main.pos = c(0.5, 1.1),                # Adjust title position
    cat.cex = 1.5,                         # Increase category label size
    cat.pos = c(-20, 20, 20, -20),         # Adjust category label positions
    cat.dist = c(0.05, 0.05, 0.05, 0.05),  # Adjust distance from circles
    height = 3000,                         # Increase diagram height
    width = 3000,                          # Increase diagram width
    resolution = 300,                      # Increase resolution
    cex = 1.5,                             # Increase font size for numbers
    fill = c("skyblue", "pink", "lightgreen", "orange"),  # Colors for the circles
    alpha = 0.5,                           # Set transparency level
    
    filename = output_filename,            # Save to file
    output = TRUE,
    log = FALSE
  )
  
  # Plot the Venn diagram if the filename is NULL
  if (is.null(output_filename)) {
    grid.draw(venn.plot)
  }
}

#STAGES vs TARBASE
create_venn_diagram_four(
  dataframes = list(group2_upregulated, tarbase_2_up, group2_downregulated, tarbase_2_down),
  labels = c("Group 2 Upregulated", "Tarbase 2 Up", "Group 2 Downregulated", "Tarbase 2 Down"),
  title = "Venn Diagram: Group 2 Upregulated, Tarbase 2 Up, Group 2 Downregulated, Tarbase 2 Down",
  output_filename = "deseq_mrna/group2_tarbase_up_down_venn.png"
)

create_venn_diagram_four(
  dataframes = list(group3_upregulated, tarbase_3_up, group3_downregulated, tarbase_3_down),
  labels = c("Group 3 Upregulated", "Tarbase 3 Up", "Group 3 Downregulated", "Tarbase 3 Down"),
  title = "Venn Diagram: Group 3 Upregulated, Tarbase 3 Up, Group 3 Downregulated, Tarbase 3 Down",
  output_filename = "deseq_mrna/group3_tarbase_up_down_venn.png"
)

create_venn_diagram_four(
  dataframes = list(group4_upregulated, tarbase_4_up, group4_downregulated, tarbase_4_down),
  labels = c("Group 4 Upregulated", "Tarbase 4 Up", "Group 4 Downregulated", "Tarbase 4 Down"),
  title = "Venn Diagram: Group 4 Upregulated, Tarbase 4 Up, Group 4 Downregulated, Tarbase 4 Down",
  output_filename = "deseq_mrna/group4_tarbase_up_down_venn.png"
)

#STAGES vs MIRTARBASE
create_venn_diagram_four(
  dataframes = list(group2_upregulated, mirtarbase_2_up, group2_downregulated, mirtarbase_2_down),
  labels = c("Group 2 Upregulated", "mirtarbase 2 Up", "Group 2 Downregulated", "mirtarbase 2 Down"),
  title = "Venn Diagram: Group 2 Upregulated, mirtarbase 2 Up, Group 2 Downregulated, mirTarbase 2 Down",
  output_filename = "deseq_mrna/group2_mirtarbase_up_down_venn.png"
)

create_venn_diagram_four(
  dataframes = list(group3_upregulated, mirtarbase_3_up, group3_downregulated, mirtarbase_3_down),
  labels = c("Group 3 Upregulated", "mirtarbase 3 Up", "Group 3 Downregulated", "mirtarbase 3 Down"),
  title = "Venn Diagram: Group 3 Upregulated, mirtarbase 3 Up, Group 3 Downregulated, mirTarbase 3 Down",
  output_filename = "deseq_mrna/group3_mirtarbase_up_down_venn.png"
)

create_venn_diagram_four(
  dataframes = list(group4_upregulated, mirtarbase_4_up, group4_downregulated, mirtarbase_4_down),
  labels = c("Group 4 Upregulated", "mirtarbase 4 Up", "Group 4 Downregulated", "mirtarbase 4 Down"),
  title = "Venn Diagram: Group 4 Upregulated, mirtarbase 4 Up, Group 4 Downregulated, mirTarbase 4 Down",
  output_filename = "deseq_mrna/group4_mirtarbase_up_down_venn.png"
)

# Find the overlap between group2_upregulated and tarbase_2_up``
overlap_uptar <- intersect(group2_upregulated$hgnc_symbol, tarbase_2_up$hgnc_symbol)
overlap_upmir <- intersect(group2_upregulated$hgnc_symbol, mirtarbase_2_up$hgnc_symbol)

# Find the overlap between group2_upregulated and tarbase_2_down
overlap_downtar <- intersect(group2_upregulated$hgnc_symbol, tarbase_2_down$hgnc_symbol)
overlap_downmir <- intersect(group2_upregulated$hgnc_symbol, mirtarbase_2_down$hgnc_symbol)

# Combine the two overlaps into one unique list
combined_overlap_tar <- unique(c(overlap_uptar, overlap_downtar))
combined_overlap_mir <- unique(c(overlap_upmir, overlap_downmir))

tar_mir_up_2 <- intersect(combined_overlap_tar, combined_overlap_mir)
# Convert to a data frame to save as a CSV
overlap_df <- data.frame(hgnc_symbol = combined_overlap)

# Save the overlapping genes to a CSV file
write.csv(overlap_df, "deseq_mrna/group2_overlap_tarbase_up_down.csv", row.names = FALSE)

