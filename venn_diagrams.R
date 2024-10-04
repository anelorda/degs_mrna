library(VennDiagram)
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)



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
overlap_downtar <- intersect(group2_downregulated$hgnc_symbol, tarbase_2_down$hgnc_symbol)
overlap_downmir <- intersect(group2_downregulated$hgnc_symbol, mirtarbase_2_down$hgnc_symbol)


create_venn_diagram_overlap <- function(dataframes, labels, title, output_filename) {
  # Ensure the number of dataframes/lists and labels match
  if (length(dataframes) != 4 || length(labels) != 4) {
    stop("This function requires exactly four dataframes/lists and four labels.")
  }
  
  # Process each item in the dataframes/lists
  symbol_lists <- lapply(dataframes, function(item) {
    if (is.data.frame(item)) {
      # If item is a data frame, filter out NA values and extract hgnc_symbol
      item %>%
        filter(!is.na(hgnc_symbol)) %>%
        pull(hgnc_symbol)
    } else if (is.vector(item) && is.character(item)) {
      # If item is a list of genes (character vector), use it directly
      item
    } else {
      stop("Each item in 'dataframes' must be either a data frame or a character vector of genes.")
    }
  })
  
  # Adjust labels to fit inside the plot area
  adjust_labels <- function(label) {
    if (nchar(label) > 15) {
      words <- strsplit(label, " ")[[1]]
      mid <- ceiling(length(words) / 2)
      label <- paste(words[1:mid], collapse = " ")
      label <- paste(label, "\n", paste(words[(mid+1):length(words)], collapse = " "))
    }
    return(label)
  }
  
  adjusted_labels <- sapply(labels, adjust_labels)
  
  # Define colors for the circles and labels
  fill_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  label_colors <- fill_colors
  
  # Create the Venn diagram with four components
  venn.plot <- venn.diagram(
    x = symbol_lists,
    category.names = adjusted_labels,  # Use adjusted labels here
    main = title,
    
    # Customizations for Venn diagram
    main.cex = 2.5,                          # Increase title size
    main.pos = c(0.2, 1),                # Adjust title position
    main.fontface = "bold",
    lwd = 2,                               # Line width
    lty = 'blank',                         # Line type
    fill = fill_colors,                    # Fill colors for the circles
    alpha = 0.5,                           # Transparency level
    cex = 1.5,                             # Font size for numbers
    fontface = "italic",                   # Font style for numbers
    
    # Customizing label appearance
    cat.cex = 1.5,                         # Size of the category labels
    cat.fontface = "bold",                 # Font face for category labels
    cat.col = label_colors,                # Colors for the category labels
    cat.default.pos = "outer",             # Place labels outside the circles
    cat.dist = c(0.21, 0.21, -0.1, -0.1),    # Adjust distance from circles
    cat.pos = c(-7, 7, 180, 180),        # Adjust label positions
    
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


# find overlap between mirtarbase and tarbase
overlap_db_2_down <- intersect(mirtarbase_2_down$hgnc_symbol, tarbase_2_down$hgnc_symbol)
overlap_db_2_up <- intersect(mirtarbase_2_up$hgnc_symbol, tarbase_2_up$hgnc_symbol)
overlap_db_3_down <- intersect(mirtarbase_3_down$hgnc_symbol, tarbase_3_down$hgnc_symbol)
overlap_db_3_up <- intersect(mirtarbase_3_up$hgnc_symbol, tarbase_3_up$hgnc_symbol)
overlap_db_4_down <- intersect(mirtarbase_4_down$hgnc_symbol, tarbase_4_down$hgnc_symbol)
overlap_db_4_up <- intersect(mirtarbase_4_up$hgnc_symbol, tarbase_4_up$hgnc_symbol)


# STAGE 2 OVERLAP
create_venn_diagram_overlap(
  dataframes = list(group2_upregulated, overlap_db_2_up, group2_downregulated, overlap_db_2_down),
  labels = c("DEGs upregulated", "upregulated miRNA targets", "DEGs downregulated", "upregulated miRNA targets"),
  title = "A. Stage 2",
  output_filename = "deseq_mrna/group2_overlap_up_down_venn.png"
)


create_venn_diagram_overlap(
  dataframes = list(group3_upregulated, overlap_db_3_up, group3_downregulated, overlap_db_3_down),
  labels = c("DEGs upregulated", "upregulated miRNA targets", "DEGs downregulated", "upregulated miRNA targets"),
  title = "B. Stage 3",
  output_filename = "deseq_mrna/group3_overlap_up_down_venn.png"
)

create_venn_diagram_overlap(
  dataframes = list(group4_upregulated, overlap_db_4_up, group4_downregulated, overlap_db_4_down),
  labels = c("DEGs upregulated", "upregulated miRNA targets", "DEGs downregulated", "upregulated miRNA targets"),
  title = "C. Stage 4",
  output_filename = "deseq_mrna/group4_overlap_up_down_venn.png"
)

# Combine the two overlaps into one unique list
combined_overlap_tar <- unique(c(overlap_uptar, overlap_downtar))
combined_overlap_mir <- unique(c(overlap_upmir, overlap_downmir))

tar_mir_up_2 <- intersect(combined_overlap_tar, combined_overlap_mir)

filtered_df <- res_group2_vs_control_df %>%
  filter(hgnc_symbol %in% tar_mir_up_2)
# Convert to a data frame to save as a CSV
overlap_df <- data.frame(hgnc_symbol = combined_overlap)

# Save the overlapping genes to a CSV file
write.csv(overlap_df, "deseq_mrna/group2_overlap_tarbase_up_down.csv", row.names = FALSE)


#TABLE------------------------------------------------------------------------------------------------

# Define the DEGs DataFrames
group_dataframes <- list(
  group2_upregulated = group2_upregulated$hgnc_symbol,
  group2_downregulated = group2_downregulated$hgnc_symbol,
  group3_upregulated = group3_upregulated$hgnc_symbol,
  group3_downregulated = group3_downregulated$hgnc_symbol,
  group4_upregulated = group4_upregulated$hgnc_symbol,
  group4_downregulated = group4_downregulated$hgnc_symbol
)

# Define the overlap DataFrames
overlap_dataframes <- list(
  overlap_db_2_up = overlap_db_2_up,
  overlap_db_2_down = overlap_db_2_down,
  overlap_db_3_up = overlap_db_3_up,
  overlap_db_3_down = overlap_db_3_down,
  overlap_db_4_up = overlap_db_4_up,
  overlap_db_4_down = overlap_db_4_down
)

# Initialize an empty result table
result_overlap <- data.frame(matrix(ncol = length(group_dataframes), nrow = length(overlap_dataframes)))

# Set column names based on the DEGs DataFrames
colnames(result_overlap) <- names(group_dataframes)

# Set row names based on the overlap DataFrames
rownames(result_overlap) <- names(overlap_dataframes)

# Fill in the table with counts of intersecting genes
for (i in seq_along(overlap_dataframes)) {
  for (j in seq_along(group_dataframes)) {
    intersection_count <- length(intersect(overlap_dataframes[[i]], group_dataframes[[j]]))
    result_overlap[i, j] <- ifelse(intersection_count > 0, intersection_count, "none")
  }
}

# View the final table
print(result_overlap)


#MIRTARBASE
# Define the overlap DataFrames
mirtarbase_dataframes <- list(
  mirtarbase_2_up = mirtarbase_2_up$hgnc_symbol,
  mirtarbase_2_down = mirtarbase_2_down$hgnc_symbol,
  mirtarbase_3_up = mirtarbase_3_up$hgnc_symbol,
  mirtarbase_3_down = mirtarbase_3_down$hgnc_symbol,
  mirtarbase_4_up = mirtarbase_4_up$hgnc_symbol,
  mirtarbase_4_down = mirtarbase_4_down$hgnc_symbol
)

# Initialize an empty result table
result_mirtarbase <- data.frame(matrix(ncol = length(group_dataframes), nrow = length(mirtarbase_dataframes)))

# Set column names based on the DEGs DataFrames
colnames(result_mirtarbase) <- names(group_dataframes)

# Set row names based on the overlap DataFrames
rownames(result_mirtarbase) <- names(mirtarbase_dataframes)

# Fill in the table with counts of intersecting genes
for (i in seq_along(mirtarbase_dataframes)) {
  for (j in seq_along(group_dataframes)) {
    intersection_count <- length(intersect(mirtarbase_dataframes[[i]], group_dataframes[[j]]))
    result_mirtarbase[i, j] <- ifelse(intersection_count > 0, intersection_count, "none")
  }
}

# View the final table
print(result_mirtarbase)

#TARBASE
# Define the overlap DataFrames
tarbase_dataframes <- list(
  tarbase_2_up = tarbase_2_up$hgnc_symbol,
  tarbase_2_down = tarbase_2_down$hgnc_symbol,
  tarbase_3_up = tarbase_3_up$hgnc_symbol,
  tarbase_3_down = tarbase_3_down$hgnc_symbol,
  tarbase_4_up = tarbase_4_up$hgnc_symbol,
  tarbase_4_down = tarbase_4_down$hgnc_symbol
)

# Initialize an empty result table
result_tarbase <- data.frame(matrix(ncol = length(group_dataframes), nrow = length(tarbase_dataframes)))

# Set column names based on the DEGs DataFrames
colnames(result_tarbase) <- names(group_dataframes)

# Set row names based on the overlap DataFrames
rownames(result_tarbase) <- names(tarbase_dataframes)

# Fill in the table with counts of intersecting genes
for (i in seq_along(tarbase_dataframes)) {
  for (j in seq_along(group_dataframes)) {
    intersection_count <- length(intersect(tarbase_dataframes[[i]], group_dataframes[[j]]))
    result_tarbase[i, j] <- ifelse(intersection_count > 0, intersection_count, "none")
  }
}

# View the final table
print(result_tarbase)

#TABLE------------------------------------------------------------------------------------------



#SAVE TABLE ***********************************************************************************

save_table_as_image <- function(result_table, output_filename, col_title, row_title) {

  
  # Create the table as a grid object
  table_grob <- tableGrob(result_table)
  
  # Add the column and row titles
  title_grob <- textGrob(col_title, gp = gpar(fontsize = 16))
  row_title_grob <- textGrob(row_title, gp = gpar(fontsize = 16), rot = 90)
  
  # Arrange the table with the titles
  final_plot <- grid.arrange(
    arrangeGrob(
      row_title_grob, 
      arrangeGrob(table_grob, nrow = 1),
      widths = c(0.2, 1)
    ),
    top = title_grob
  )
  
  # Save the plot as an image
  ggsave(output_filename, plot = final_plot, width = 10, height = 8)
}

# Example usage
save_table_as_image(result_overlap, "deseq_mrna/result_overlap_image.png", "DEGs from mRNA-seq", "miRNA targets from both databases miRTarbase and TarBase")
save_table_as_image(result_mirtarbase, "deseq_mrna/result_mirtarbase_image.png", "DEGs from mRNA-seq", "miRNA targets from miRTarbase")
save_table_as_image(result_tarbase, "deseq_mrna/result_tarbase_image.png", "DEGs from mRNA-seq", "miRNA targets from TarBase")

#SAVE TABLE ***********************************************************************************



#stage 2 down 
# Find overlapping genes between overlap_db_2_down and group4_upregulated
overlap_with_group2_up <- intersect(overlap_db_2_down, group2_upregulated$hgnc_symbol)
overlap_with_group2_down <- intersect(overlap_db_2_down, group2_downregulated$hgnc_symbol)

# Combine this with the existing overlap_2_db_down_downregulated list
com_overlap_2_down <- unique(c(overlap_with_group2_up, overlap_with_group2_down))

filtered_mirna_stage2_down <- mirna_stage2_down %>%
  filter(target_symbol %in% com_overlap_2_down)

unique_mature_mirna_ids2_down <- unique(filtered_mirna_stage2_down$mature_mirna_id)


uni#stage2 up
# Find overlapping genes between overlap_db_2_down and group4_upregulated
overlap_with_group2_1_up <- intersect(overlap_db_2_up, group2_upregulated$hgnc_symbol)
overlap_with_group2_1_down <- intersect(overlap_db_2_up, group2_downregulated$hgnc_symbol)

# Combine this with the existing overlap_2_db_down_downregulated list
comb_overlap_2_up <- unique(c(overlap_with_group2_1_up, overlap_with_group2_1_down))

filtered_mirna_stage2_up <- mirna_stage2_up %>%
  filter(target_symbol %in% comb_overlap_2_up)

unique_mature_mirna_ids2_up <- unique(filtered_mirna_stage2_up$mature_mirna_id)
