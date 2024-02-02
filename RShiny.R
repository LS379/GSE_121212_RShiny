library(shiny)
library(plotly)
library(DESeq2)
library(edgeR)
library(tidyr)


# Define UI for application
ui <- fluidPage(
  titlePanel("Differential Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("option1", "Select Option 1", choices = NULL), # Choices will be updated server-side
      selectInput("option2", "Select Option 2", choices = NULL), # Choices will be updated server-side
      selectizeInput("selected_genes", "Select Genes", choices = c("IL13", "IL22", "IL33", "CCL17", "CCL18", "IL10", "IL19", "IL20", "IL17A", "IFNG", "IL36G", "CXCL1", "CXCL6", "TNIP3", "FOXE1", "CLEC3A"), multiple = TRUE),
      sliderInput("p_value_cutoff", "P-Value Cutoff", min = 0, max = 1, value = 0.05, step = 0.01),
      sliderInput("expression_cutoff", "Expression Cutoff", min = 0, max = 2, value = 0.58, step = 0.01)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Gene Expression Data", plotlyOutput("gene_plot")),
        tabPanel("Contrastive Analysis", plotlyOutput("contrastive_plot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  data <- readRDS("~/path/GSE121212_preformatted.rds")
  count_mat <- data[[1]]
  meta <- data[[2]]
  
  # Update choices for selectInput based on 'meta'
  updateSelectInput(session, "option1", choices = unique(meta$subcat))
  updateSelectInput(session, "option2", choices = unique(meta$subcat))
  
  get_gene_expression_data <- function(selected_genes, option1, option2) {
    # Assuming 'meta' and 'count_mat' are available here
    smeta <- meta[meta$subcat %in% c(option1, option2),]
    scount <- count_mat[, colnames(count_mat) %in% smeta$sampID]
    smeta$subcat <- factor(smeta$subcat, levels = c(option1, option2))
    
    # Normalize gene expression
    scount_norm <- edgeR::cpm(scount)
    
    # Subset on selected genes
    scount_norm <- as.data.frame(scount_norm[rownames(scount_norm) %in% selected_genes, ])
    scount_norm$Gene <- rownames(scount_norm)
    
    # Prepare data for plotting
    gel <- tidyr::pivot_longer(scount_norm, cols = -Gene, names_to = "Sample", values_to = 'Expression')
    gel$Contrast <- smeta$subcat[match(gel$Sample, smeta$sampID)]
    return(gel)
  }
  
  # Gene expression plot
  output$gene_plot <- renderPlotly({
    gel <- get_gene_expression_data(input$selected_genes, input$option1, input$option2)
    plot_ly(gel, x = ~Gene, y = ~Expression, type = 'bar', color = ~Contrast) %>% layout(barmode = 'group', showlegend = TRUE)
  })
  get_contrastive_analysis_data <- function(option1, option2) {
    data <- readRDS("~/path/GSE121212_preformatted.rds")
    count_mat <- data[[1]]
    meta <- data[[2]]
    smeta <- meta[meta$subcat %in% c(option1,option2),] 
    scount <- count_mat[,colnames(count_mat) %in% smeta$sampID]
    smeta$subcat <- factor(smeta$subcat, levels = c(option1,option2))
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = scount,
                                  colData = smeta,
                                  design= ~ subcat)
    #Run DESeq2
    dds <- DESeq(dds)
    
    # Get differential expression results
    results <- as.data.frame(results(dds))
    #print(results)
    results$gene <- rownames(results)
    # Highlight genes based on p-value and expression cutoffs
    #data$Color <- ifelse(data$PValue < p_value_cutoff & abs(data$Expression) > expression_cutoff, "red", "grey")
    
    return(results)
  }
  
  # Adjusted contrastive analysis plot section
  output$contrastive_plot <- renderPlotly({
    # Call 'get_contrastive_analysis_data' to get the results based on user input
    results <- get_contrastive_analysis_data(input$option1, input$option2)
    
    # Ensure 'results' dataframe has the expected structure
    if ("padj" %in% names(results) && "log2FoldChange" %in% names(results) && "gene" %in% names(results)) {
      # Generate the ggplot
      p <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj),
                               color = abs(log2FoldChange) > input$expression_cutoff & padj < input$p_value_cutoff,
                               label = gene)) +
        geom_point() +
        geom_hline(yintercept = -log10(input$p_value_cutoff), color = "red") +
        geom_vline(xintercept = c(-input$expression_cutoff, input$expression_cutoff), color = "red") +
        guides(color = "none") +
        labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
        theme_minimal()
      
      # Convert the ggplot object to a Plotly object for interactivity
      plotly::ggplotly(p)
    } else {
      # If 'results' doesn't have the expected structure, return an informative message
      plotly::plot_ly() %>%
        add_annotations(text = "Data not available or unexpected format", x = 0.5, y = 0.5, showarrow = F, font = list(size = 20))
    }
  })
}

shinyApp(ui, server)
