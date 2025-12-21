# DOT PLOT MODULE | R/module_dot_plot.R ----

#' UI for Dot Plot Controls
#' @description Creates the UI for the dot plot controls in the sidebar.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 selectInput checkboxInput h5 hr p em actionButton icon br downloadButton
#' @importFrom DT dataTableOutput
#' 
dotPlotControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("2. Plotting Options"),
    
    # Step 1: Select grouping variable and axis orientation
    selectInput(
      ns("groupingVar"),
      "1. Group cells by:",
      choices = NULL # Choices will be updated from the server
    ),
    checkboxInput(ns("swapAxes"), "Swap Axes (Genes on X-axis)", value = FALSE),
    
    # Clustering Options
    h5("Clustering Options"),
    checkboxInput(ns("clusterRows"), "Cluster Genes", value = FALSE),
    checkboxInput(ns("clusterCols"), "Cluster Groups", value = FALSE),
    
    hr(),
    
    # Step 2: Select genes from the master list
    h5("2. Select Genes from Master List"),
    p(em("Click a row to add a gene to the selection below.")),
    DT::dataTableOutput(ns("geneTable")),
    
    hr(),
    
    # Step 3: Manage the list of genes to plot
    h5("3. Genes Selected for Plot"),
    p(em("This table contains the genes that will be plotted.")),
    DT::dataTableOutput(ns("selectionTable")),
    actionButton(ns("removeGenes"), "Remove from Selection", icon = icon("trash-alt"), class = "btn-light w-100 mt-2"),
    
    hr(),
    
    # Step 4: Generate and download the plot
    h5("4. Generate Plot"),
    actionButton(
      ns("generateDotPlot"),
      "Generate Dot Plot",
      class = "btn-primary w-100",
      icon = icon("chart-bar")
    ),
    br(),
    downloadButton(ns("downloadDotPlot"), "Download Plot", class = "w-100 mt-2")
  )
}

#' UI for Dot Plot Output Area
#' @description Creates the UI for the dot plot output in the main panel.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList plotOutput
#'
dotPlotOutputUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("dotPlot"), height = "800px", width = "100%")
  )
}


#' Server for the Dot Plot Module.
#' @description Handles the logic for generating and displaying a dot plot based on user selections.
#' @param id A namespace ID.
#' @param loaded_data A reactive list from the dataLoaderServer module.
#' @importFrom shiny moduleServer NS reactiveVal observeEvent req updateSelectInput validate need showNotification reactive renderPlot downloadHandler
#' @importFrom DT renderDataTable datatable
#' @importFrom purrr map_dfr
#' @importFrom dplyr left_join mutate group_by summarise n
#' @importFrom data.table as.data.table dcast
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit grid.circle
#' @importFrom ComplexHeatmap Heatmap Legend draw
#' @importFrom grDevices png dev.off
#' @importFrom stats setNames
#' 
dotPlotServer <- function(id, loaded_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    selected_genes_rv <- reactiveVal(data.frame())
    
    observeEvent(loaded_data$obs_data, {
      req(loaded_data$obs_data, loaded_data$config)
      
      categorical_cols <- names(loaded_data$obs_data)[sapply(loaded_data$obs_data, function(x) is.character(x) || is.factor(x))]
      selected_col <- if (!is.null(loaded_data$config$obs$columns$celltype) && loaded_data$config$obs$columns$celltype %in% categorical_cols) {
        loaded_data$config$obs$columns$celltype
      } else {
        categorical_cols[1]
      }
      updateSelectInput(session, "groupingVar", choices = categorical_cols, selected = selected_col)
      
      selected_genes_rv(data.frame())
    })
    
    output$geneTable <- DT::renderDataTable({
      req(loaded_data$gene_table)
      DT::datatable(
        loaded_data$gene_table,
        selection = 'single',
        rownames = FALSE,
        options = list(pageLength = 50, scrollY = "200px", info = FALSE)
      )
    })
    
    output$selectionTable <- DT::renderDataTable({
      validate(need(nrow(selected_genes_rv()) > 0, "No genes selected. Please add genes from the master list above."))
      DT::datatable(
        selected_genes_rv(),
        selection = 'multiple',
        rownames = FALSE,
        options = list(pageLength = 50, scrollY = "150px", searching = FALSE, info = FALSE)
      )
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
    
    plot_data_r <- eventReactive(input$generateDotPlot, {
      req(loaded_data$obs_data, input$groupingVar, loaded_data$config)
      
      selection_df <- selected_genes_rv()
      
      validate(
        need(nrow(selection_df) > 1, "Validation Error: Please select at least two genes to generate the plot."),
        need("var_names" %in% names(selection_df), "Error: The gene table must contain a 'var_names' column to load data.")
      )
      
      genes_to_load <- selection_df$var_names 
      
      showNotification("Processing genes...", type = "message", duration = 3)
      
      dot_plot_data <- purrr::map_dfr(genes_to_load, function(gene) {
        gene_data <- loaded_data$get_gene_data(gene)
        
        if (is.null(gene_data)) {
          showNotification(paste("Warning: Could not load data for gene", gene), type = "warning")
          return(NULL)
        }
        
        joined_data <- loaded_data$obs_data %>%
          left_join(gene_data, by = setNames(loaded_data$config$quant$barcode, loaded_data$config$obs$columns$barcode))
        
        # Handle different expression data types
        slot_type <- loaded_data$config$quant$slot
        if (!is.null(slot_type) && slot_type == "counts") {
          processed_data <- joined_data %>%
            mutate(scount = log2((expression / .data[[loaded_data$config$obs$columns$total_counts]]) * 1e4) + 1) %>%
            mutate(gene_count = ifelse(is.na(scount), 0, scount))
        } else {
          processed_data <- joined_data %>%
            mutate(gene_count = ifelse(is.na(expression), 0, expression))
        }
        
        processed_data %>%
          group_by(group = .data[[input$groupingVar]]) %>%
          dplyr::summarise(
            pct_expressing = (sum(gene_count > 0) / n()) * 100,
            avg_exp = mean(gene_count[gene_count > 0]),
            .groups = 'drop'
          ) %>%
          mutate(avg_exp = ifelse(is.nan(avg_exp), 0, avg_exp)) %>%
          mutate(gene = gene)
      })
      
      validate(need(nrow(dot_plot_data) > 0, "Could not generate plot. Verify that 'var_names' match the gene data filenames."))
      
      return(dot_plot_data)
    })
    
    reactive_dot_plot <- reactive({
      req(plot_data_r())
      df <- plot_data_r()
      
      # 1. Reshape data into matrices
      dt <- data.table::as.data.table(df)
      exp_mat_dt <- data.table::dcast(dt, gene ~ group, value.var = "avg_exp", fill = 0)
      pct_mat_dt <- data.table::dcast(dt, gene ~ group, value.var = "pct_expressing", fill = 0)
      
      gene_order <- exp_mat_dt$gene
      exp_mat <- as.matrix(exp_mat_dt[, -1]); rownames(exp_mat) <- gene_order
      pct_mat <- as.matrix(pct_mat_dt[, -1]); rownames(pct_mat) <- gene_order
      pct_mat <- pct_mat[rownames(exp_mat), colnames(exp_mat)] # Ensure perfect order
      
      # 2. Handle axis swapping if requested
      if (isTRUE(input$swapAxes)) {
        exp_mat <- t(exp_mat)
        pct_mat <- t(pct_mat)
      }
      
      # 3. Handle clustering logic based on final matrix dimensions
      n_rows_in_mat <- nrow(exp_mat)
      n_cols_in_mat <- ncol(exp_mat)
      
      if (isTRUE(input$swapAxes)) {
        cluster_the_rows <- isTRUE(input$clusterCols) && (n_rows_in_mat > 1)
        cluster_the_cols <- isTRUE(input$clusterRows) && (n_cols_in_mat > 1)
      } else {
        cluster_the_rows <- isTRUE(input$clusterRows) && (n_rows_in_mat > 1)
        cluster_the_cols <- isTRUE(input$clusterCols) && (n_cols_in_mat > 1)
      }
      
      # 4. Define aesthetics
      slot_type <- loaded_data$config$quant$slot
      legend_title <- "Avg. Expression"
      
      if (!is.null(slot_type) && slot_type == "scale.data") {
        legend_title <- "Avg. Scaled Exp"
        max_abs_val <- max(abs(exp_mat), na.rm = TRUE)
        col_fun <- circlize::colorRamp2(c(-max_abs_val, 0, max_abs_val), c("blue", "grey90", "red"))
      } else {
        if (!is.null(slot_type) && slot_type == "data") legend_title <- "Avg. Log-Norm Exp"
        max_exp <- max(exp_mat, na.rm = TRUE)
        if (max_exp == 0) max_exp <- 1
        col_fun <- circlize::colorRamp2(c(0, max_exp), c("grey90", "#440154FF")) # Viridis
      }
      
      # 5. Custom cell drawing function
      cell_fun <- function(j, i, x, y, width, height, fill) {
        if (exp_mat[i, j] > 0 || pct_mat[i, j] > 0) {
          grid::grid.circle(
            x = x, y = y,
            r = unit(sqrt(pct_mat[i, j] / 100) * 1.5, "mm"),
            gp = grid::gpar(fill = col_fun(exp_mat[i, j]), col = NA)
          )
        }
      }
      
      # 6. Create Heatmap object
      ht <- ComplexHeatmap::Heatmap(
        exp_mat,
        cluster_rows = cluster_the_rows,
        cluster_columns = cluster_the_cols,
        show_row_dend = cluster_the_rows,
        show_column_dend = cluster_the_cols,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        column_names_rot = 45,
        show_heatmap_legend = FALSE
      )
      
      # 7. Create custom legends
      lgd_exp <- ComplexHeatmap::Legend(
        title = legend_title, 
        col_fun = col_fun,
        grid_width = unit(1, "mm"),
        grid_height = unit(2, "mm"),
        title_gp = gpar(fontsize = 3),
        labels_gp = gpar(fontsize = 3)
      )
      lgd_pct <- ComplexHeatmap::Legend(
        title = "% Expressing",
        at = c(0, 25, 50, 75, 100),
        type = "points",
        pch = 16,
        size = unit(sqrt(c(0, 25, 50, 75, 100) / 100) * 3, "mm"),
        legend_gp = gpar(col = "black"),
        title_gp = gpar(fontsize = 3),
        labels_gp = gpar(fontsize = 3)
      )
      
      # 8. Draw the heatmap and legends together
      return(ComplexHeatmap::draw(ht,
                                  heatmap_legend_list = list(lgd_exp, lgd_pct),
                                  merge_legend = TRUE,
                                  heatmap_legend_side = "right"
      ))
    })
    
    output$dotPlot <- renderPlot({
      reactive_dot_plot()
    }, res = 300)
    
    output$downloadDotPlot <- downloadHandler(
      filename = function() {
        paste0("dot_plot_", input$groupingVar, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(reactive_dot_plot())
        p <- reactive_dot_plot()
        
        final_mat_dims <- dim(p@ht_list[[1]]@matrix)
        n_rows <- final_mat_dims[1]
        n_cols <- final_mat_dims[2]
        
        height <- max(6, 0.35 * n_rows)
        width <- max(8, 0.4 * n_cols)
        
        png(file, width = width, height = height, units = "in", res = 300)
        grid::grid.draw(p)
        dev.off()
      }
    )
  })
}