#' @title scGeyser
#'
#' @description Run the shiny app for single-cell data visualization.
#'
#' @export
#'
#' @importFrom shiny NS actionButton br checkboxInput column conditionalPanel
#' @importFrom shiny downloadButton downloadHandler em eventReactive fluidRow
#' @importFrom shiny h4 h5 hr HTML icon isolate moduleServer need
#' @importFrom shiny numericInput fileInput  uiOutput
#' @importFrom shiny observeEvent p plotOutput 
#' @importFrom shiny reactive reactiveVal renderPlot renderText req
#' @importFrom shiny selectInput tagList tags updateSelectInput
#' @importFrom shiny validate verbatimTextOutput showNotification removeNotification
#' @importFrom shiny nearPoints reactiveValues callModule shinyApp
#' @importFrom shinyjs useShinyjs
#' @import tidyverse
#' @importFrom magrittr "%>%"
#' @import yaml
#' @importFrom data.table fread data.table
#' @import scattermore
#' @import ggrepel
#' @import ggh4x
#' @import pals
#' @import cowplot
#' @importFrom bslib nav_menu nav_panel nav_spacer navset_card_tab page_fluid 
#' @importFrom bslib page_navbar sidebar bs_theme
#' @importFrom DT renderDataTable datatable dataTableOutput
#' @import ComplexHeatmap
#' @import rlang
#' @import shinyFiles
#'
scGeyser <- function(mgif_dir = NULL) {
  # ----------------- UI -----------------
  ui <- bslib::page_navbar(
    shinyjs::useShinyjs(),
    id = "main_nav",
    title = "scGeyser",
    tags$style(HTML('table.dataTable tr.selected td, table.dataTable tr.selected {box-shadow: inset 0 0 0 9999px #400A91 !important;}')),
    theme = bs_theme(primary = "#400a91",
                     secondary = "#1F0940",
                     fg = "#400a91",bg = "#fff",
                     base_font = "Arial",
                     code_font = "Courier",
                     font_scale = 0.8,
                     spacer = "0.5rem", preset = "cosmo"),
    sidebar = sidebar(
      width = 350,
      conditionalPanel(
        condition = "input.main_nav == 'Data Load'",
        dataLoaderUI("data_loader"),
        hr()
      ),
      conditionalPanel(
        condition = "input.main_nav == '2D Explorer'",
        umapPlotControlsUI("umap_plotter")
      ),
      conditionalPanel(
        condition = "input.main_nav == 'Dot Plot'",
        dotPlotControlsUI("dot_plotter")
      ),
      conditionalPanel(
        condition = "input.main_nav == 'Expression Plot'",
        expressionPlotControlsUI("expression_plotter")
      ),
      conditionalPanel(
        condition = "input.main_nav == 'Table Explorer'",
        tableExplorerControlsUI("table_explorer")
      )
    ),
    nav_panel(
      title = "Data Load",
      value = "Data Load",
      dataLoaderOutputUI("data_loader")
    ),
    nav_panel(
      title = "2D Explorer",
      value = "2D Explorer",
      umapPlotOutputUI("umap_plotter")
    ),
    nav_panel(
      title = "Dot Plot",
      value = "Dot Plot",
      dotPlotOutputUI("dot_plotter")
    ),
    nav_panel(
      title = "Expression Plot",
      value = "Expression Plot",
      expressionPlotOutputUI("expression_plotter")
    ),
    nav_panel(
      title = "Table Explorer",
      value = "Table Explorer",
      tableExplorerOutputUI("table_explorer")
    )
  )
  
  # ----------------- Server -----------------
  server <- function(input, output, session) {
    observeEvent(input$inputType, {
      if (input$inputType == "rds") {
        if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
          showModal(modalDialog(
            title = "Missing Packages",
            "To load Seurat or SingleCellExperiment objects, please install the following packages: 'Seurat', 'SingleCellExperiment'",
            easyClose = TRUE
          ))
        }
      } else if (input$inputType == "hdf5") {
        if (!requireNamespace("HDF5Array", quietly = TRUE)) {
          showModal(modalDialog(
            title = "Missing Package",
            "To load HDF5-backed SCE objects, please install the 'HDF5Array' package.",
            easyClose = TRUE
          ))
        }
      }
    })
    
    loaded_data <- dataLoaderServer("data_loader", mgif_dir = mgif_dir)
    umapPlotServer("umap_plotter", loaded_data)
    dotPlotServer("dot_plotter", loaded_data)
    expressionPlotServer("expression_plotter", loaded_data)
    tableExplorerServer("table_explorer", loaded_data)
  }
  
  
  # ----------------- Run App -----------------
  shinyApp(ui = ui, server = server)
}