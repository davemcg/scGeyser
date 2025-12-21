# TABLE EXPLORER MODULE | R/module_table_explorer.R ----

#' UI for Table Explorer Controls
#' @description Creates the UI for the table explorer controls in the sidebar.
#' @param id A namespace ID.
#' @importFrom shiny NS tagList h4 selectInput numericInput hr h5 actionButton icon
#'
tableExplorerControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("2. Table Options"),
    
    # Grouping variable selection
    selectInput(
      ns("groupingVars"),
      "1. Group cells by (order matters):",
      choices = NULL,
      multiple = TRUE
    ),
    
    # Threshold for filtering results
    numericInput(
      ns("ratioThreshold"),
      "2. Minimum Ratio Threshold:",
      value = 0.01, min = 0, max = 1, step = 0.01
    ),
    
    hr(),
    h5("3. Select a Gene (optional)"),
    # Gene selection table (re-using the single-select gene table module)
    geneTableUI(ns("gene_selector")),
    
    hr(),
    # Action button to trigger table generation
    actionButton(
      ns("generateTables"),
      "Generate Tables",
      icon = icon("table-cells"),
      class = "btn-primary w-100"
    )
  )
}

#' UI for Table Explorer Output
#' @description Creates the UI for the table explorer output in the main panel.
#' @param id A namespace ID.
#' @importFrom shiny NS fluidRow column h4 p em
#' @importFrom DT DTOutput
#'
tableExplorerOutputUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    # Column for the metadata summary table
    column(
      width = 6,
      h4("Metadata Summary"),
      p(em("Shows cell counts and ratios for the selected groupings across all cells.")),
      DT::DTOutput(ns("metadataSummaryTable"))
    ),
    # Column for the gene expression summary table
    column(
      width = 6,
      h4("Gene Expression Summary"),
      p(em("Shows the same summary, but only for cells expressing the selected gene.")),
      DT::DTOutput(ns("geneSummaryTable"))
    )
  )
}

#' Server for the Table Explorer Module
#' @description Handles the logic for summarizing metadata and gene expression in tables.
#' @param id A namespace ID.
#' @param loaded_data A reactive list from the dataLoaderServer module.
#' @importFrom shiny moduleServer NS observeEvent req updateSelectInput eventReactive showNotification validate need
#' @importFrom DT renderDataTable datatable formatRound
#' @importFrom rlang syms
#' @importFrom dplyr group_by summarise n arrange left_join mutate
#' @importFrom stats setNames
#'
tableExplorerServer <- function(id, loaded_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- Module Wiring ---
    # Use the geneTable module for single gene selection
    selected_gene_r <- geneTableServer("gene_selector", reactive(loaded_data$gene_table))
    
    # --- UI Updates ---
    # Update grouping choices when new data is loaded
    observeEvent(loaded_data$obs_data, {
      req(loaded_data$obs_data)
      categorical_cols <- names(loaded_data$obs_data)[sapply(loaded_data$obs_data, function(x) is.character(x) || is.factor(x))]
      # Pre-select the first two categorical columns as a default
      updateSelectInput(session, "groupingVars", choices = categorical_cols, selected = '')
    })
    
    # --- Core Reactive Logic ---
    # A reactive to hold the results for the metadata summary table
    metadata_summary_r <- eventReactive(input$generateTables, {
      req(loaded_data$obs_data, input$groupingVars)
      
      showNotification("Generating metadata summary table...", type = "message", duration = 3)
      # Use the provided function to summarize the full dataset
      summary_df <- summarize_grouped_data(
        data = loaded_data$obs_data,
        !!!syms(input$groupingVars),
        threshold = input$ratioThreshold
      )
      
      return(summary_df)
    }, ignoreNULL = FALSE)
    
    # A reactive to hold the results for the gene summary table
    gene_summary_r <- eventReactive(input$generateTables, {
      req(
        loaded_data$obs_data,
        loaded_data$config,
        input$groupingVars,
        selected_gene_r()
      )
      
      gene_name <- selected_gene_r()
      validate(need(gene_name, "Please select a gene from the list."))
      showNotification(paste("Generating summary for", gene_name, "..."), type = "message", duration = 3)
      
      # Load the data for the selected gene
      gene_data <- loaded_data$get_gene_data(gene_name)
      validate(need(!is.null(gene_data), paste("Could not load data for gene:", gene_name)))
      
      total_counts_col <- loaded_data$config$obs$columns$total_counts
      
      # Join metadata with gene data, keeping all cells
      merged_data <- loaded_data$obs_data %>%
        left_join(gene_data, by = setNames(loaded_data$config$quant$barcode, loaded_data$config$obs$columns$barcode))
      
      # Conditionally calculate scaled expression
      slot_type <- loaded_data$config$quant$slot
      if (!is.null(slot_type) && slot_type == "counts") {
        validate(need(total_counts_col %in% names(loaded_data$obs_data), "Total counts column not found for normalization."))
        merged_data <- merged_data %>%
          mutate(
            expression = ifelse(is.na(expression), 0, expression),
            scaled_expression = (expression / .data[[total_counts_col]]) * 1e4
          )
      } else {
        merged_data <- merged_data %>%
          mutate(
            expression = ifelse(is.na(expression), 0, expression),
            scaled_expression = expression # Use value directly
          )
      }
      
      # Group and summarize
      summary_df <- merged_data %>%
        group_by(!!!syms(input$groupingVars)) %>%
        summarise(
          `Cells # Detected` = sum(expression > 0), 
          `Total Cells` = n(),
          mean_expression = mean(scaled_expression, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(!!!syms(input$groupingVars))
      
      return(summary_df)
      
    }, ignoreNULL = FALSE)
    
    # --- Render Outputs ---
    output$metadataSummaryTable <- DT::renderDataTable({
      DT::datatable(
        metadata_summary_r(),
        options = list(pageLength = 25, scrollX = TRUE, info = TRUE),
        rownames = FALSE,
        caption = "Summary of all cells"
      )
    })
    
    output$geneSummaryTable <- DT::renderDataTable({
      DT::datatable(
        gene_summary_r(),
        options = list(pageLength = 25, scrollX = TRUE, info = TRUE),
        rownames = FALSE,
        caption = paste("Summary for cells expressing", selected_gene_r())
      ) %>% DT::formatRound('mean_expression', 3) 
    })
    
  })
}