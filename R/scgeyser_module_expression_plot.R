# EXPRESSION PLOT MODULE | R/module_expression_plot.R ----

#' UI for Expression Plot Controls
#' @description Creates the UI for the expression plot controls in the sidebar.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 selectInput checkboxInput hr h5 p em actionButton icon br downloadButton
#' @importFrom DT dataTableOutput
#'
expressionPlotControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("2. Plotting Options"),
    
    # Main plot settings
    selectInput(ns("xaxis"), "1. Group by (X-axis):", choices = NULL),
    selectInput(ns("color"), "2. Color by:", choices = NULL),
    selectInput(ns("facet"), "3. Facet by:", choices = NULL, multiple = TRUE),
    selectInput(ns("plotType"), "4. Plot type:",
                choices = c("Violin", "Boxplot"),
                selected = "Boxplot" 
    ),
    
    checkboxInput(ns("flipAxes"), "Flip Axes", value = FALSE),
    
    hr(),
    h5("Gene Selection"),
    p(em("Click a row to add a gene to the selection below.")),
    DT::dataTableOutput(ns("geneTable")),
    
    hr(),
    h5("Genes to Plot"),
    DT::dataTableOutput(ns("selectionTable")),
    actionButton(ns("removeGenes"), "Remove from Selection", icon = icon("trash-alt"), class = "btn-light w-100 mt-2"),
    
    hr(),
    h5("Generate & Download"),
    actionButton(
      ns("generatePlot"),
      "Generate Plot",
      class = "btn-primary w-100",
      icon = icon("chart-bar")
    ),
    br(),
    downloadButton(ns("downloadPlot"), "Download Plot", class = "w-100 mt-2")
  )
}

#' UI for Expression Plot Output Area
#' @description Creates the UI for the expression plot output in the main panel.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList plotOutput
#'
expressionPlotOutputUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("expressionPlot"), height = "800px", width = "100%")
  )
}

#' Server for the Expression Plot Module
#' @description Handles the logic for generating and displaying violin or box plots.
#' @param id Namespace ID.
#' @param loaded_data A reactive list from the dataLoaderServer module.
#' @importFrom shiny moduleServer NS reactiveVal observeEvent req updateSelectInput validate need showNotification eventReactive renderPlot downloadHandler
#' @importFrom DT renderDataTable datatable
#' @importFrom purrr map_dfr
#' @importFrom dplyr left_join mutate select
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot facet_wrap coord_flip scale_fill_manual labs theme element_line element_text element_rect
#' @importFrom pals trubetskoy glasbey
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
expressionPlotServer <- function(id, loaded_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    selected_genes_rv <- reactiveVal(data.frame())
    
    # --- UI Updates ---
    observeEvent(loaded_data$obs_data, {
      req(loaded_data$obs_data, loaded_data$config)
      
      # Get categorical columns for UI selectors
      categorical_cols <- names(loaded_data$obs_data)[sapply(loaded_data$obs_data, function(x) is.character(x) || is.factor(x))]
      
      selected_xaxis <- loaded_data$config$obs$columns$celltype
      if (is.null(selected_xaxis) || !selected_xaxis %in% categorical_cols) {
        selected_xaxis <- categorical_cols[1]
      }
      
      updateSelectInput(session, "xaxis", choices = categorical_cols, selected = selected_xaxis)
      updateSelectInput(session, "color", choices = c("None", categorical_cols), selected = "None")
      updateSelectInput(session, "facet", choices = c("None", categorical_cols), selected = "None")
      
      selected_genes_rv(data.frame())
    })
    
    # --- Gene Selection Logic ---
    output$geneTable <- DT::renderDataTable({
      req(loaded_data$gene_table)
      DT::datatable(loaded_data$gene_table, selection = 'single', rownames = FALSE,
                    options = list(pageLength = 25, scrollY = "150px", info = FALSE))
    })
    
    output$selectionTable <- DT::renderDataTable({
      validate(need(nrow(selected_genes_rv()) > 0, "No genes selected."))
      DT::datatable(selected_genes_rv(), selection = 'multiple', rownames = FALSE,
                    options = list(pageLength = 25, scrollY = "150px", info = FALSE, searching = FALSE))
    })
    
    observeEvent(input$geneTable_rows_selected, {
      req(loaded_data$gene_table, input$geneTable_rows_selected)
      newly_selected <- loaded_data$gene_table[input$geneTable_rows_selected, ]
      current_selection <- selected_genes_rv()
      combined <- rbind(current_selection, newly_selected)
      unique_selection <- combined[!duplicated(combined[[1]]), ]
      selected_genes_rv(unique_selection)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    observeEvent(input$removeGenes, {
      req(nrow(selected_genes_rv()) > 0, input$selectionTable_rows_selected)
      rows_to_remove <- input$selectionTable_rows_selected
      current_selection <- selected_genes_rv()
      updated_selection <- current_selection[-rows_to_remove, ]
      selected_genes_rv(updated_selection)
    })
    
    # --- Core Data Processing ---
    plot_data_r <- eventReactive(input$generatePlot, {
      req(loaded_data$obs_data, input$xaxis, loaded_data$config)
      
      # Validation for UI inputs
      # ... (omitted for brevity)
      
      selection_df <- selected_genes_rv()
      validate(
        need(nrow(selection_df) > 0, "Please select at least one gene."),
        need("var_names" %in% names(selection_df), "Gene table requires 'var_names' column.")
      )
      
      genes_to_load <- selection_df$var_names
      showNotification("Loading and processing gene data...", type = "message", duration = 4)
      
      grouping_vars <- unique(c(input$xaxis, input$color[input$color != "None"], input$facet[input$facet != "None"]))
      
      all_gene_data <- purrr::map_dfr(genes_to_load, function(gene) {
        gene_df <- loaded_data$get_gene_data(gene)
        if (is.null(gene_df)) return(NULL)
        
        joined_data <- loaded_data$obs_data %>%
          left_join(gene_df, by = setNames(loaded_data$config$quant$barcode, loaded_data$config$obs$columns$barcode)) %>%
          mutate(expression = ifelse(is.na(expression), 0, expression))
        
        # Conditionally normalize data based on the selected slot
        slot_type <- loaded_data$config$quant$slot
        if (!is.null(slot_type) && slot_type == "counts") {
          final_data <- joined_data %>%
            mutate(scount = log2((expression / .data[[loaded_data$config$obs$columns$total_counts]]) * 1e4 + 1))
        } else {
          # Assume data from "data" or "scale.data" slots is already processed
          final_data <- joined_data %>% mutate(scount = expression)
        }
        
        final_data %>%
          select(all_of(grouping_vars), scount) %>%
          mutate(gene = gene)
      })
      
      validate(need(nrow(all_gene_data) > 0, "Could not load data for selected genes."))
      all_gene_data$gene <- factor(all_gene_data$gene, levels = genes_to_load)
      
      return(all_gene_data)
    })
    
    # --- Plot Rendering ---
    reactive_exp_plot <- eventReactive(input$generatePlot, {
      plot_df <- plot_data_r()
      
      base_aes <- aes(x = .data[[input$xaxis]], y = scount)
      if (input$color != "None") {
        base_aes$fill <- as.name(input$color)
      }
      
      p <- ggplot(plot_df, base_aes)
      
      # Add geom layer
      geom_args <- list(lwd = 0.1)
      if (input$plotType == "Violin") {
        geom_args$scale <- "width"; geom_args$trim <- TRUE
        if(input$color == "None") geom_args$fill <- "grey80" 
        p <- p + do.call(geom_violin, geom_args)
      } else {
        geom_args$outlier.shape <- NA; geom_args$fatten <- 2
        if(input$color == "None") geom_args$fill <- "grey80"
        p <- p + do.call(geom_boxplot, geom_args)
      }
      
      # Add facets
      facet_vars <- input$facet[input$facet != "None"]
      p <- p + facet_wrap(c("gene", facet_vars), scales = "free_y")
      
      if (isTRUE(input$flipAxes)) p <- p + coord_flip()
      if (input$color != "None") p <- p + scale_fill_manual(values = unname(rep(c(pals::trubetskoy(), pals::glasbey()), 10)))
      
      # Dynamic axis label
      slot_type <- loaded_data$config$quant$slot
      y_label <- if (!is.null(slot_type) && slot_type == "counts") {
        "log2(CP10K + 1)"
      } else {
        paste("Expression (", slot_type, ")", sep="")
      }
      
      caption_text <- paste("scGeyser", loaded_data$data_source_name)
      plot_labs <- if (isTRUE(input$flipAxes)) {
        labs(y = input$xaxis, x = y_label, caption = caption_text)
      } else {
        labs(x = input$xaxis, y = y_label, caption = caption_text)
      }
      if (input$color != "None") plot_labs$fill <- input$color
      
      p <- p + plot_labs + get_core_plot_theme() +
        theme(
          axis.line = element_line(linewidth = 0.2),
          axis.ticks = element_line(linewidth = 0.2),
          axis.text = element_text(size = 3),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "right"
        )
      
      return(p)
    })
    
    output$expressionPlot <- renderPlot({
      reactive_exp_plot()
    }, res = 300)
    
    output$downloadPlot <- downloadHandler(
      filename = function() { paste0("expression_plot_", Sys.Date(), ".png") },
      content = function(file) {
        req(reactive_exp_plot())
        ggsave(file, plot = reactive_exp_plot(), width = 12, height = 8, dpi = 300)
      }
    )
  })
}