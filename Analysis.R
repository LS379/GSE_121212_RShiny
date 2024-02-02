# Load count matrix from a text file using fread for fast reading, and convert it to a data frame.
count_mat <- as.data.frame(data.table::fread("~/path/GSE121212_readcount.txt"))
# Gene names to be excluded from the analysis bc of wrong format name
gene_exclusion_list <- c(
  "1-Dec",
  "1-Mar",
  "1-Sep",
  "10-Mar",
  "10-Sep",
  "11-Sep",
  "12-Sep",
  "15-Sep",
  "2-Mar",
  "2-Sep",
  "3-Mar",
  "3-Sep",
  "4-Mar",
  "4-Sep",
  "5-Mar",
  "5-Sep",
  "6-Mar",
  "6-Sep",
  "7-Mar",
  "7-Sep",
  "8-Mar",
  "8-Sep",
  "9-Mar",
  "9-Sep"
)
# Exclude rows from count_mat that match any gene name in the exclusion list.
count_mat <- count_mat[-which(count_mat$V1 %in% gene_exclusion_list),]
# Set the rownames of count_mat to the first column, which contains gene names, then remove this column.
rownames(count_mat) <- count_mat$V1
count_mat$V1 <- NULL
# Create a metadata frame with sample IDs extracted from the column names of count_mat.
meta <- data.frame(sampID = colnames(count_mat))
# Extract disease state from sampID and assign it to a new column 'disease'.
meta$disease <- unlist(lapply(strsplit(meta$sampID,split = "_", fixed = TRUE), "[[", 1))
meta <- meta[-1,]
# Extract tissue type from sampID, with some string manipulation to ensure consistent formatting
meta$tissue <- unlist(lapply(strsplit(gsub("chronic_lesion", "chronic-lesion", meta$sampID, fixed = TRUE),split = "_", fixed = TRUE), "[[", 3))
meta$subcat <- paste0(meta$disease, "_", meta$tissue)
# Save the processed count matrix and metadata as an RDS file
saveRDS(list(count_mat, meta), "~/path/GSE121212_preformatted.rds")
# Function to perform differential expression analysis
get_DEGs <- function(c1,c2) {
  # c1 <- "CTRL_healthy"
  # c2 <- "AD_lesional"
  data <- readRDS("~/path/GSE121212_preformatted.rds")
  count_mat <- data[[1]]
  meta <- data[[2]]
  # Filter metadata and count data for the two conditions of interest
  smeta <- meta[meta$subcat %in% c(c1,c2),] 
  scount <- count_mat[,colnames(count_mat) %in% smeta$sampID]
  smeta$subcat <- factor(smeta$subcat, levels=c(c1,c2))
  
  dim(scount)
  dim(smeta)
  
  library(DESeq2)
  # Initialize a DESeq2 dataset object
  dds <- DESeqDataSetFromMatrix(countData = scount,
                                colData = smeta,
                                design= ~ subcat)
  #Run DESeq2
  dds <- DESeq(dds)
  
  # Get differential expression results
  results <- as.data.frame(results(dds))
  print(results)
  results$gene <- rownames(results)
  
  return(results)
}


data <- readRDS("~/path/GSE121212_preformatted.rds")
count_mat <- data[[1]]
meta <- data[[2]]
possible_options <- unique(meta$subcat)
# Perform DEG analysis between the first two options.
results <- get_DEGs(possible_options[1],possible_options[2]) 

library(ggplot2)
# Visualize the results using ggplot2 and convert it to an interactive plot
p <- ggplot(results, aes(x=log2FoldChange, y=-log10(padj),
                         color = abs(log2FoldChange)>0.58 & padj <0.05,
                         label=gene )) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = c(-0.58, 0.58))
plotly::ggplotly(
  p = p
)



