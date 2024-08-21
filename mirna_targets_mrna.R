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
  
  # Create the Venn diagram with adjustments
  venn.plot <- venn.diagram(
    x = symbol_lists,
    category.names = labels,
    main = title,
    
    # Adjustments to ensure text isn't cut off
    main.cex = 2,                          # Increase title size
    main.pos = c(0.5, 1.1),                # Adjust title position
    cat.cex = 1.5,                         # Increase category label size
    cat.pos = c(-20, 20),                  # Adjust category label positions
    cat.dist = c(0.05, 0.05),              # Adjust distance from circles
    height = 3000,                         # Increase diagram height
    width = 3000,                          # Increase diagram width
    resolution = 300,                      # Increase resolution
    cex = 1.5,                             # Increase font size for numbers
    fill = c("skyblue", "pink"),           # Colors for the circles
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
# Load your data
mirna_stage2_down <- read.table("Downloads/multimir_results_stage2_down.csv", header = TRUE, sep = ",")
mirna_stage3_down <- read.table("Downloads/multimir_results_stage3_down.csv", header = TRUE, sep = ",")
mirna_stage4_down <- read.table("Downloads/multimir_results_stage4_down.csv", header = TRUE, sep = ",")
mirna_stage2_up <- read.table("Downloads/multimir_results_stage2_up.csv", header = TRUE, sep = ",")
mirna_stage3_up <- read.table("Downloads/multimir_results_stage3_up.csv", header = TRUE, sep = ",")
mirna_stage3_up <- read.table("Downloads/multimir_results_stage4_up.csv", header = TRUE, sep = ",")

# Filter for "mirtarbase" entries
mirtarbase_2_down <- mirna_stage2_down %>% filter(database == "mirtarbase")
mirtarbase_3_down <- mirna_stage3_down %>% filter(database == "mirtarbase")
mirtarbase_4_down <- mirna_stage4_down %>% filter(database == "mirtarbase")
mirtarbase_2_up <- mirna_stage2_up %>% filter(database == "mirtarbase")
mirtarbase_3_up <- mirna_stage3_up %>% filter(database == "mirtarbase")
mirtarbase_4_up <- mirna_stage4_up %>% filter(database == "mirtarbase")

# Filter for "mirtarbase" entries
tarbase_2_down <- mirna_stage2_down %>% filter(database == "tarbase")
tarbase_3_down <- mirna_stage3_down %>% filter(database == "tarbase")
tarbase_4_down <- mirna_stage4_down %>% filter(database == "tarbase")
tarbase_2_up <- mirna_stage2_up %>% filter(database == "tarbase")
tarbase_3_up <- mirna_stage3_up %>% filter(database == "tarbase")
tarbase_4_up <- mirna_stage4_up %>% filter(database == "tarbase")

mirtarbase_2_down$hgnc_symbol <- mirtarbase_2_down$target_symbol
mirtarbase_2_down$target_symbol <- NULL
mirtarbase_3_down$hgnc_symbol <- mirtarbase_3_down$target_symbol
mirtarbase_3_down$target_symbol <- NULL
mirtarbase_4_down$hgnc_symbol <- mirtarbase_4_down$target_symbol
mirtarbase_4_down$target_symbol <- NULL
mirtarbase_2_up$hgnc_symbol <- mirtarbase_2_up$target_symbol
mirtarbase_2_up$target_symbol <- NULL
mirtarbase_3_up$hgnc_symbol <- mirtarbase_3_up$target_symbol
mirtarbase_3_up$target_symbol <- NULL
mirtarbase_4_up$hgnc_symbol <- mirtarbase_4_up$target_symbol
mirtarbase_4_up$target_symbol <- NULL

# Rename target_symbol to hgnc_symbol and remove target_symbol column for targetscan
tarbase_2_down$hgnc_symbol <- tarbase_2_down$target_symbol
tarbase_2_down$target_symbol <- NULL
tarbase_3_down$hgnc_symbol <- tarbase_3_down$target_symbol
tarbase_3_down$target_symbol <- NULL
tarbase_4_down$hgnc_symbol <- tarbase_4_down$target_symbol
tarbase_4_down$target_symbol <- NULL
tarbase_2_up$hgnc_symbol <- tarbase_2_up$target_symbol
tarbase_2_up$target_symbol <- NULL
tarbase_3_up$hgnc_symbol <- tarbase_3_up$target_symbol
tarbase_3_up$target_symbol <- NULL
tarbase_4_up$hgnc_symbol <- tarbase_4_up$target_symbol
tarbase_4_up$target_symbol <- NULL


# Load gene lists
group2_downregulated <- read.table("deseq_mrna/downregulated_group2.csv", header = TRUE, sep = ",")
group3_downregulated <- read.table("deseq_mrna/downregulated_group3.csv", header = TRUE, sep = ",")
group4_downregulated <- read.table("deseq_mrna/downregulated_group4.csv", header = TRUE, sep = ",")
group2_upregulated <- read.table("deseq_mrna/upregulated_group2.csv", header = TRUE, sep = ",")
group3_upregulated <- read.table("deseq_mrna/upregulated_group3.csv", header = TRUE, sep = ",")
group4_upregulated <- read.table("deseq_mrna/upregulated_group4.csv", header = TRUE, sep = ",")


# Create the Venn diagram for the three stages with mirtarbase DOWNREGULATED
create_venn_diagram(
  dataframes = list(mirtarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 2, downregulated",
  output_filename = "deseq_mrna/2_down_mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 3, downregulated",
  output_filename = "deseq_mrna/3_down_mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 4, downregulated",
  output_filename = "deseq_mrna/4_down_mirnatargets_vs_DEGs.png"
)

# Create the Venn diagram for the three stages with tarbase
create_venn_diagram(
  dataframes = list(tarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (tarbase)","mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 2, downregulated",
  output_filename = "deseq_mrna/2_down_tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 3, downregulated",
  output_filename = "deseq_mrna/3_down_tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 4, downregulated",
  output_filename = "deseq_mrna/4_down_tarbase_vs_DEGs.png"
)


#UPREGULATED
create_venn_diagram(
  dataframes = list(mirtarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 2, upregulated",
  output_filename = "deseq_mrna/2_down_mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 3, upregulated",
  output_filename = "deseq_mrna/3_down_mirnatargets_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(mirtarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (mirtarbase)", "mRNA"),
  title = "miRNA targets (mirtarbase) and DEGs for stage 4, upregulated",
  output_filename = "deseq_mrna/4_down_mirnatargets_vs_DEGs.png"
)

# Create the Venn diagram for the three stages with tarbase
create_venn_diagram(
  dataframes = list(tarbase_2_down, group2_downregulated),
  labels = c("miRNA targets (tarbase)","mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 2, upregulated",
  output_filename = "deseq_mrna/2_down_tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_3_down, group3_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 3, upregulated",
  output_filename = "deseq_mrna/3_down_tarbase_vs_DEGs.png"
)

create_venn_diagram(
  dataframes = list(tarbase_4_down, group4_downregulated),
  labels = c("miRNA targets (tarbase)", "mRNA"),
  title = "miRNA targets (tarbase) and DEGs for stage 4, upregulated",
  output_filename = "deseq_mrna/4_down_tarbase_vs_DEGs.png"
)
#mirtarbase vs tarbase
create_venn_diagram(
  dataframes = list(tarbase_2_up, mirtarbase_2_up),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 2 up",
  output_filename = "deseq_mrna/2_up_tarbase_vs_mirtarbase.png"
)

create_venn_diagram(
  dataframes = list(tarbase_3_up, mirtarbase_3_up),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 3 up",
  output_filename = "deseq_mrna/3_up_tarbase_vs_mirtarbase.png"
)

create_venn_diagram(
  dataframes = list(tarbase_3_down, mirtarbase_3_down),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 3 down",
  output_filename = "deseq_mrna/3_down_tarbase_vs_mirtarbase.png"
)

create_venn_diagram(
  dataframes = list(tarbase_2_down, mirtarbase_2_down),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 2 down",
  output_filename = "deseq_mrna/2_down_tarbase_vs_mirtarbase.png"
)

create_venn_diagram(
  dataframes = list(tarbase_4_up, mirtarbase_4_up),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 4 up",
  output_filename = "deseq_mrna/4_up_tarbase_vs_mirtarbase.png"
)

create_venn_diagram(
  dataframes = list(tarbase_4_down, mirtarbase_4_down),
  labels = c("tarbase", "mirtarbase"),
  title = "tarbase vs mirtarbase 4 down",
  output_filename = "deseq_mrna/4_down_tarbase_vs_mirtarbase.png"
)

top100_group2_upregulated <- group2_upregulated %>%
  arrange(desc(log2FoldChange)) %>%     
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

top100_group3_upregulated <- group3_upregulated %>%
  arrange(desc(log2FoldChange)) %>%           
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

top100_group4_upregulated <- group4_upregulated %>%
  arrange(desc(log2FoldChange)) %>%           
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

top100_group2_downregulated <- group2_downregulated %>%
  arrange(desc(log2FoldChange)) %>%           
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

top100_group3_downregulated <- group3_downregulated %>%
  arrange(desc(log2FoldChange)) %>%           
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

top100_group4_downregulated <- group4_downregulated %>%
  arrange(desc(log2FoldChange)) %>%           
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  head(100) 

 overlap2<- intersect(mirtarbase_2_down$hgnc_symbol, group2_downregulated$hgnc_symbol)
 overlap3<- intersect(mirtarbase_3_down$hgnc_symbol, group3_downregulated$hgnc_symbol)
 overlap4<- intersect(mirtarbase_4_down$hgnc_symbol, group4_downregulated$hgnc_symbol)
 overlap2_3_up <- intersect(top100_group2_upregulated$hgnc_symbol, top100_group3_upregulated$hgnc_symbol)
 overlap2_3_up 
 overlap3_4_up <- intersect(top100_group3_upregulated$hgnc_symbol, top100_group4_upregulated$hgnc_symbol)
 overlap3_4_up 
 overlap2_4_up <- intersect(top100_group2_upregulated$hgnc_symbol, top100_group4_upregulated$hgnc_symbol)
 overlap2_4_up 
 overlap_groups_up <- Reduce(intersect, list(top100_group2_upregulated$hgnc_symbol, top100_group3_upregulated$hgnc_symbol, top100_group4_upregulated$hgnc_symbol))
 overlap_groups_up
 
 # Create a list of HGNC symbols
 venn_list <- list(
   top100_group2_upregulated$hgnc_symbol,
   top100_group3_upregulated$hgnc_symbol,
   top100_group4_upregulated$hgnc_symbol
 )
 
 # Create Venn diagram
 venn.plot <- venn.diagram(
   x = venn_list,
   category.names = c("Group2", "Group3", "Group4"),
   filename = NULL,
   output = TRUE
 )

 grid.draw(venn.plot)

 
 