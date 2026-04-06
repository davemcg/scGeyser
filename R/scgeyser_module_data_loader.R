# DATA LOADER MODULE | R/module_data_loader.R ----

#' UI for the Data Loader Module (Sidebar Controls)
#' @description Creates the UI for the data loader module, which appears in the sidebar.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 selectInput conditionalPanel fileInput br tags verbatimTextOutput div hr h5 p em actionButton
#' @importFrom shinyjs hidden
#' @importFrom shinyFiles shinyDirButton
dataLoaderUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("User Data"),
    selectInput(ns("inputType"), "Select Input Type",
                choices = c("MGIF Project" = "mgif",
                            "Seurat or SCE Object (.rds)" = "rds",
                            "HDF5-backed SCE (Folder)" = "hdf5")),
    
    conditionalPanel(
      condition = "input.inputType == 'mgif'", ns = ns,
      fileInput(ns("configFile"), "Upload MGIF YAML File", accept = c(".yaml", ".yml"))
    ),
    
    conditionalPanel(
      condition = "input.inputType == 'rds'", ns = ns,
      fileInput(ns("rdsFile"), "Upload .rds File", accept = c(".rds"))
    ),
    
    conditionalPanel(
      condition = "input.inputType == 'hdf5'", ns = ns,
      shinyFiles::shinyDirButton(ns("hdf5Dir"), "Select Folder", "Select HDF5-backed SCE Folder"),
      br(), br(),
      tags$strong("Selected folder:"),
      verbatimTextOutput(ns("hdf5Path"), placeholder = TRUE)
    ),
    
    uiOutput(ns("local_files_ui")),
    
    shinyjs::hidden(
      div(id = ns("rds_mapping_div"),
          hr(),
          h5("Map Object Data"),
          p(em("Map the data from your object to the required fields.")),
          selectInput(ns("rdsEmbedding1"), "Embedding 1", choices = NULL),
          selectInput(ns("rdsEmbedding2"), "Embedding 2", choices = NULL),
          selectInput(ns("rdsBarcode"), "Cell Barcode Column", choices = NULL),
          selectInput(ns("rdsTotalCounts"), "Total Counts Column", choices = NULL),
          selectInput(ns("rdsDefaultCluster"), "Default Annotation Column", choices = NULL),
          selectInput(ns("rdsAssaySlot"), "Assay/Slot for Expression", choices = NULL),
          actionButton(ns("processRds"), "Load Object Data", class = "btn-primary w-100")
      )
    )
  )
}

#' UI for the Data Loader Output Area (Main Panel)
#' @description Creates the UI for the data loader's output, which appears in the main panel.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 hr
#' @importFrom DT DTOutput
dataLoaderOutputUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("mgif_header")), 
    h4("Cell Metadata Preview (First 10 Rows)"),
    div(style = "margin-bottom: 20px;",
        DT::DTOutput(ns("obsDataPreview"))),
    hr(),
    h4("Gene/Variable Metadata Preview (First 10 Rows)"),
    div(style = "margin-bottom: 20px;",
        DT::DTOutput(ns("varDataPreview"))),
    hr(),
    h4("YAML File"),
    div(style = "margin-bottom: 20px;",
        DT::DTOutput(ns("yaml_contents")))
  )
}


#' Server for the Data Loader Module
#' @description Handles the logic for loading data from different sources like MGIF, RDS, and HDF5.
#' @param id Namespace ID.
#' @param mgif_dir Optional local directory for MGIF files.
#' @return A reactive list containing `$config`, `$obs_data`, `$gene_table`, and a data access function `$get_gene_data`.
dataLoaderServer <- function(id, mgif_dir = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to store loaded data and intermediate objects
    data_out <- reactiveValues(
      config = NULL,
      obs_data = NULL,
      gene_table = NULL,
      get_gene_data = NULL,
      data_source_name = NULL
    )
    raw_rds_object <- reactiveVal(NULL) 
    temp_data_source_name <- reactiveVal(NULL)
    
    # Render local files UI if directory provided
    output$local_files_ui <- renderUI({
      req(mgif_dir)
      if (!dir.exists(mgif_dir)) return(NULL)
      
      available_files <- list.files(mgif_dir, 
                                    pattern = "\\.(yaml|yml|mgif)$", 
                                    full.names = TRUE, 
                                    ignore.case = TRUE)
      
      if (length(available_files) == 0) return(NULL)
      
      names(available_files) <- basename(available_files)
      
      tagList(
        h4("Computer Datasets"),
        selectInput(ns("local_file_select"), "Select YAML:", choices = available_files),
        actionButton(ns("load_local_btn"), "Load Local MGIF", 
                     class = "btn-primary w-100", icon = icon("hdd"))
      )
    })
    
    # Logic to handle local file loading
    observeEvent(input$load_local_btn, {
      req(input$local_file_select)
      path <- input$local_file_select
      
      tryCatch({
        config <- yaml::read_yaml(path)
        data_out$config <- config
        data_out$obs_data <- data.table::fread(config$obs$file)
        data_out$gene_table <- data.table::fread(config$var$file)
        
        mgif_path <- config$mgifPath
        data_out$get_gene_data <- function(gene_name) {
          data.table::fread(paste0(mgif_path, '/', gene_name, '.csv.gz'))
        }
        
        data_out$data_source_name <- basename(path)
        showNotification(paste("Loaded:", basename(path)), type = "message")
      }, error = function(e) {
        showModal(modalDialog(title = "Error Loading Local File", tags$pre(e$message)))
      })
    })
    
    # Helper function to extract metadata and reductions
    extract_full_metadata <- function(obj) {
      if (inherits(obj, "Seurat")) {
        meta_data <- obj@meta.data
        reduction_names <- names(obj@reductions)
        
        all_embeddings <- lapply(reduction_names, function(reduc_name) {
          embeddings <- as.data.frame(obj@reductions[[reduc_name]]@cell.embeddings)
          n_dims <- min(ncol(embeddings), 20)
          if (n_dims == 0) return(NULL)
          embeddings_to_add <- embeddings[, 1:n_dims, drop = FALSE]
          colnames(embeddings_to_add) <- paste0(toupper(reduc_name), "_", 1:n_dims)
          return(embeddings_to_add)
        })
        
      } else if (inherits(obj, "SingleCellExperiment")) {
        meta_data <- as.data.frame(SummarizedExperiment::colData(obj))
        reduction_names <- SingleCellExperiment::reducedDimNames(obj)
        
        all_embeddings <- lapply(reduction_names, function(reduc_name) {
          embeddings <- as.data.frame(SingleCellExperiment::reducedDim(obj, reduc_name))
          n_dims <- min(ncol(embeddings), 20)
          if (n_dims == 0) return(NULL)
          embeddings_to_add <- embeddings[, 1:n_dims, drop = FALSE]
          colnames(embeddings_to_add) <- paste0(toupper(reduc_name), "_", 1:n_dims)
          return(embeddings_to_add)
        })
      } else {
        stop("Object is not a recognized Seurat or SingleCellExperiment.")
      }
      
      all_embeddings <- all_embeddings[!sapply(all_embeddings, is.null)]
      if (length(all_embeddings) > 0) {
        full_meta <- do.call(cbind, c(list(meta_data), all_embeddings))
        colnames(full_meta) <- make.unique(colnames(full_meta))
        return(full_meta)
      } else {
        return(meta_data)
      }
    }
    
    # MGIF Loading Logic
    observeEvent(input$configFile, {
      req(input$configFile)
      tryCatch({
        config <- read_yaml(input$configFile$datapath)
        obs_data <- data.table::fread(config$obs$file)
        
        for (col_name in names(obs_data)) {
          if (is.numeric(obs_data[[col_name]])) {
            num_unique <- data.table::uniqueN(obs_data[[col_name]])
            if (num_unique > 1 && num_unique < 1000) {
              obs_data[, (col_name) := as.factor(get(col_name))]
            }
          }
        }
        
        gene_table <- data.table::fread(config$var$file)
        mgif_path <- config$mgifPath
        
        get_gene_data_mgif <- function(gene_name) {
          file_path <- paste0(mgif_path, '/', gene_name, '.csv.gz')
          tryCatch({ data.table::fread(file_path) }, error = function(e) return(NULL))
        }
        
        data_out$config <- config
        data_out$obs_data <- obs_data
        data_out$gene_table <- gene_table
        data_out$get_gene_data <- get_gene_data_mgif
        data_out$data_source_name <- input$configFile$name
        
        shinyjs::hide("rds_mapping_div")
        raw_rds_object(NULL)
        showNotification("MGIF configuration loaded.", type = "message")
      }, error = function(e) {
        showModal(modalDialog(title = "Error Loading MGIF", tags$pre(e$message), easyClose = TRUE))
      })
    })
    
    # Helper for processing Seurat/SCE objects
    process_loaded_object <- function(obj) {
      meta_data <- extract_full_metadata(obj)
      all_cols <- colnames(meta_data)
      numeric_cols <- names(meta_data)[sapply(meta_data, is.numeric)]
      factor_cols <- names(meta_data)[sapply(meta_data, function(x) is.character(x) || is.factor(x))]
      
      find_default <- function(choices, patterns, fallback) {
        match <- unlist(lapply(patterns, function(p) grep(p, choices, value = TRUE, ignore.case = TRUE)))
        if (length(match) > 0) return(match[1])
        if (fallback %in% choices) return(fallback)
        return(choices[1])
      }
      
      updateSelectInput(session, "rdsEmbedding1", choices = numeric_cols, 
                        selected = find_default(numeric_cols, c("UMAP_1", "UMAP1"), "UMAP_1"))
      updateSelectInput(session, "rdsEmbedding2", choices = numeric_cols, 
                        selected = find_default(numeric_cols, c("UMAP_2", "UMAP2"), "UMAP_2"))
      
      updateSelectInput(session, "rdsTotalCounts", choices = numeric_cols)
      updateSelectInput(session, "rdsDefaultCluster", choices = factor_cols)
      updateSelectInput(session, "rdsBarcode", choices = c("rownames", all_cols), selected = "rownames")
      
      if (inherits(obj, "Seurat")) {
        updateSelectInput(session, "rdsAssaySlot", choices = c("counts", "data", "scale.data"), selected = "data")
      } else {
        updateSelectInput(session, "rdsAssaySlot", choices = SummarizedExperiment::assayNames(obj))
      }
      
      shinyjs::show("rds_mapping_div")
    }
    
    # RDS File Observer
    observeEvent(input$rdsFile, {
      req(input$rdsFile)
      shinyjs::hide("rds_mapping_div")
      tryCatch({
        obj <- readRDS(input$rdsFile$datapath)
        temp_data_source_name(input$rdsFile$name)
        raw_rds_object(obj)
        process_loaded_object(obj)
      }, error = function(e) {
        showModal(modalDialog(title = "RDS Error", tags$pre(e$message)))
      })
    })
    
    # HDF5 Logic
    volumes <- c(Home = fs::path_home(), getVolumes()())
    shinyFiles::shinyDirChoose(input, "hdf5Dir", roots = volumes, session = session)
    
    hdf5_path_r <- reactive({
      req(input$hdf5Dir)
      shinyFiles::parseDirPath(volumes, input$hdf5Dir)
    })
    
    output$hdf5Path <- renderText({ hdf5_path_r() })
    
    observeEvent(hdf5_path_r(), {
      path <- hdf5_path_r()
      req(length(path) > 0)
      tryCatch({
        obj <- HDF5Array::loadHDF5SummarizedExperiment(dir = path)
        temp_data_source_name(path)
        raw_rds_object(obj)
        process_loaded_object(obj)
      }, error = function(e) {
        showModal(modalDialog(title = "HDF5 Error", tags$pre(e$message)))
      })
    })
    
    # Final Object Processing
    observeEvent(input$processRds, {
      req(raw_rds_object())
      obj <- raw_rds_object()
      tryCatch({
        obs_df <- extract_full_metadata(obj)
        barcode_col <- if (input$rdsBarcode == "rownames") {
          obs_df[[".__barcode"]] <- rownames(obs_df)
          ".__barcode"
        } else { input$rdsBarcode }
        
        data_out$obs_data <- data.table::as.data.table(obs_df)
        data_out$gene_table <- if (inherits(obj, "Seurat")) {
          data.table::as.data.table(obj[[Seurat::DefaultAssay(obj)]]@meta.features, keep.rownames = "var_names")
        } else {
          data.table::as.data.table(as.data.frame(SummarizedExperiment::rowData(obj)), keep.rownames = "var_names")
        }
        
        data_out$get_gene_data <- function(gene_name) {
          # Retrieval logic (Simplified for space)
          expr <- if(inherits(obj, "Seurat")) { Seurat::GetAssayData(obj, slot=input$rdsAssaySlot)[gene_name,] }
          else { SummarizedExperiment::assay(obj, input$rdsAssaySlot)[gene_name,] }
          data.table::data.table(barcode = colnames(obj), expression = as.numeric(expr))
        }
        data_out$config <- list(
          obs = list(
            columns = list(
              celltype = input$rdsDefaultCluster,
              embedding1 = input$rdsEmbedding1,
              embedding2 = input$rdsEmbedding2,
              barcode = barcode_col,
              total_counts = input$rdsTotalCounts
            )
          ),
          quant = list(
            barcode = "barcode", # Hardcoded to match get_gene_data's return column name
            slot = input$rdsAssaySlot
          )
        )
        data_out$data_source_name <- temp_data_source_name()
      }, error = function(e) {
        showModal(modalDialog(title = "Process Error", tags$pre(e$message)))
      })
    })
    
    # Previews
    output$obsDataPreview <- DT::renderDataTable({
      req(data_out$obs_data)
      DT::datatable(head(data_out$obs_data, 10), options = list(scrollX = TRUE), rownames = FALSE)
    })
    
    output$varDataPreview <- DT::renderDataTable({
      req(data_out$gene_table)
      DT::datatable(head(data_out$gene_table, 10), options = list(scrollX = TRUE), rownames = FALSE)
    })
    
    output$yaml_contents <- DT::renderDataTable({
      req(data_out$config)
      yaml_table <- as.data.frame(unlist(data_out$config))
      DT::datatable(yaml_table, options = list(pageLength = 50, scrollX = TRUE))
    })
    
    return(data_out)
  })
}