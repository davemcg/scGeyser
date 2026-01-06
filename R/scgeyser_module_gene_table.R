# GENE TABLE MODULE | R/module_gene_table.R ----

#' UI for the Gene Table Module
#' @description A reusable UI component that displays a searchable, selectable gene table.
#' @param id Namespace ID.
#' @importFrom shiny NS
#' @importFrom DT DTOutput
#' 
geneTableUI <- function(id) {
  ns <- NS(id)
  DT::DTOutput(ns("geneTable"))
}


#' Server for the Gene Table Module
#' @description A reusable server component for the gene table that returns a list of the selected gene (both the file name and the displayed name)
#' @param id Namespace ID.
#' @param gene_data A reactive data.table of gene information.
#' @importFrom shiny moduleServer req reactive
#' @importFrom DT renderDataTable datatable
#' @return A reactive containing the name of the gene selected from the table.
#' 
geneTableServer <- function(id, gene_data, config) {
  moduleServer(id, function(input, output, session) {
    
    output$geneTable <- DT::renderDataTable({
      req(gene_data())
      DT::datatable(
        gene_data(),
        options = list(pageLength = 15, scrollX = TRUE, info = FALSE),
        rownames = FALSE,
        filter = 'top',
        selection = 'single'
      )
    })
    
    # returns a list instead of a single string
    selected_gene <- reactive({
      req(gene_data(), input$geneTable_rows_selected)
      
      # Always use 'var_names' for file lookups
      file_id <- gene_data()[["var_names"]][input$geneTable_rows_selected]
      
      # Use active_gene_column for display, fallback to var_names
      display_col <- if (!is.null(config()$var$active_gene_column)) {
        config()$var$active_gene_column
      } else {
        "var_names"
      }
      
      display_name <- gene_data()[[display_col]][input$geneTable_rows_selected]
      
      return(list(id = file_id, display = display_name))
    })
    
    return(selected_gene)
  })
}