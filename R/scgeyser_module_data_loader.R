# DATA LOADER MODULE | R/module_data_loader.R ----

#' UI for the Data Loader Module (Sidebar Controls)
#' @description Creates the UI for the data loader module, which appears in the sidebar.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 selectInput conditionalPanel fileInput br tags verbatimTextOutput div hr h5 p em actionButton
#' @importFrom shinyjs hidden
#' @importFrom shinyFiles shinyDirButton
#' 
dataLoaderUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("1. Data Source"),
    selectInput(ns("inputType"), "Select Input Type",
                choices = c("MGIF Project" = "mgif",
                            "Seurat or SCE Object (.rds)" = "rds",
                            "HDF5-backed SCE (Folder)" = "hdf5")),
    
    # Conditional UI for MGIF
    conditionalPanel(
      condition = "input.inputType == 'mgif'", ns = ns,
      fileInput(ns("configFile"), "Upload MGIF YAML File",
                accept = c(".yaml", ".yml", ".mgif"))
    ),
    
    # Conditional UI for RDS
    conditionalPanel(
      condition = "input.inputType == 'rds'", ns = ns,
      fileInput(ns("rdsFile"), "Upload .rds File", accept = c(".rds"))
    ),
    
    # Conditional UI for HDF5-backed SCE
    conditionalPanel(
      condition = "input.inputType == 'hdf5'", ns = ns,
      shinyFiles::shinyDirButton(ns("hdf5Dir"), "Select Folder", "Select HDF5-backed SCE Folder"),
      br(),br(),
      tags$strong("Selected folder:"),
      verbatimTextOutput(ns("hdf5Path"), placeholder = TRUE)
    ),
    
    # This UI section will be shown after an RDS or HDF5 file is loaded
    shinyjs::hidden(
      div(id = ns("rds_mapping_div"),
          hr(),
          h5("Map Object Data"),
          p(em("Map the data from your object to the required fields.")),
          selectInput(ns("rdsEmbedding1"), "Embedding 1 (e.g., UMAP_1)", choices = NULL),
          selectInput(ns("rdsEmbedding2"), "Embedding 2 (e.g., UMAP_2)", choices = NULL),
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
#' 
dataLoaderOutputUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Cell Metadata Preview (First 10 Rows)"),
    DT::DTOutput(ns("obsDataPreview")),
    hr(),
    h4("Gene/Variable Metadata Preview (First 10 Rows)"),
    DT::DTOutput(ns("varDataPreview"))
  )
}


#' Server for the Data Loader Module
#' @description Handles the logic for loading data from different sources like MGIF, RDS, and HDF5.
#' @param id Namespace ID.
#' @return A reactive list containing `$config`, `$obs_data`, `$gene_table`, and a data access function `$get_gene_data`.
#' @importFrom shiny moduleServer reactiveValues reactiveVal observeEvent req updateSelectInput showNotification showModal modalDialog tags renderText removeNotification reactive
#' @importFrom DT renderDataTable datatable
#' @importFrom yaml read_yaml
#' @importFrom data.table fread as.data.table data.table uniqueN
#' @importFrom shinyFiles shinyDirChoose parseDirPath getVolumes
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData assayNames assay
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom utils head
#' @importFrom fs path_home
#' 
dataLoaderServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store loaded data and intermediate objects
    data_out <- reactiveValues(
      config = NULL,
      obs_data = NULL,
      gene_table = NULL,
      get_gene_data = NULL, # Unified data access function
      data_source_name = NULL
    )
    raw_rds_object <- reactiveVal(NULL) # To hold the unprocessed RDS/HDF5 object
    temp_data_source_name <- reactiveVal(NULL)
    
    # ---- Helper function to extract metadata and all reductions ----
    extract_full_metadata <- function(obj) {
      if (inherits(obj, "Seurat")) {
        meta_data <- obj@meta.data
        reductions_list <- obj@reductions
        reduction_names <- names(reductions_list)
        
        all_embeddings <- lapply(reduction_names, function(reduc_name) {
          embeddings <- as.data.frame(reductions_list[[reduc_name]]@cell.embeddings)
          n_dims <- min(ncol(embeddings), 20)
          if (n_dims == 0) return(NULL)
          embeddings_to_add <- embeddings[, 1:n_dims, drop = FALSE]
          colnames(embeddings_to_add) <- paste0(toupper(reduc_name), "_", 1:n_dims)
          return(embeddings_to_add)
        })
        
      } else if (inherits(obj, "SingleCellExperiment")) {
        meta_data <- as.data.frame(colData(obj))
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
      
      # Filter out any NULLs (from reductions with 0 columns)
      all_embeddings <- all_embeddings[!sapply(all_embeddings, is.null)]
      
      # If there's anything to add, cbind it
      if (length(all_embeddings) > 0) {
        # Combine the original metadata with all extracted embeddings
        full_meta <- do.call(cbind, c(list(meta_data), all_embeddings))
        # Ensure final column names are unique
        colnames(full_meta) <- make.unique(colnames(full_meta))
        return(full_meta)
      } else {
        return(meta_data)
      }
    }
    
    # ---- MGIF Loading Logic ----
    observeEvent(input$configFile, {
      req(input$configFile)
      tryCatch({
        config <- read_yaml(input$configFile$datapath)
        
        # Load observation data (cell metadata)
        obs_data <- data.table::fread(config$obs$file)
        
        # Convert low-cardinality numeric metadata to factors
        for (col_name in names(obs_data)) {
          if (is.numeric(obs_data[[col_name]])) {
            num_unique <- data.table::uniqueN(obs_data[[col_name]])
            if (num_unique > 1 && num_unique < 1000) {
              obs_data[, (col_name) := as.factor(get(col_name))]
            }
          }
        }
        
        # Load variable data (gene metadata)
        gene_table <- data.table::fread(config$var$file)
        
        # Define the data retrieval function for the MGIF format
        mgif_path <- config$mgifPath
        get_gene_data_mgif <- function(gene_name) {
          file_path <- paste0(mgif_path, '/', gene_name, '.csv.gz')
          tryCatch({
            data.table::fread(file_path)
          }, error = function(e) {
            warning(paste("Could not load data for:", gene_name))
            return(NULL)
          })
        }
        
        # Store results in our reactiveValues object
        data_out$config <- config
        data_out$obs_data <- obs_data
        data_out$gene_table <- gene_table
        data_out$get_gene_data <- get_gene_data_mgif
        data_out$data_source_name <- input$configFile$name
        
        # Reset the RDS UI state if a user switches back to MGIF
        shinyjs::hide("rds_mapping_div")
        raw_rds_object(NULL)
        
        showNotification("MGIF configuration and data loaded.", type = "message", duration = 5)
        
      }, error = function(e) {
        showModal(modalDialog(title = "Error Loading MGIF Data", tags$pre(e$message), easyClose = TRUE))
      })
    })
    
    # ---- Helper for processing Seurat/SCE objects ----
    process_loaded_object <- function(obj) {
      # This function takes a loaded object (Seurat or SCE) and prepares the UI for column mapping.
      meta_data <- tryCatch({
        extract_full_metadata(obj)
      }, error = function(e) {
        stop("Failed to extract metadata and reductions: ", e$message)
      })
      
      # Get available column names for mapping
      all_cols <- colnames(meta_data)
      numeric_cols <- names(meta_data)[sapply(meta_data, is.numeric)]
      factor_cols <- names(meta_data)[sapply(meta_data, function(x) is.character(x) || is.factor(x))]
      
      # Update UI selectors with smart defaults
      updateSelectInput(session, "rdsEmbedding1", choices = numeric_cols, selected = grep("UMAP_1", numeric_cols, value = TRUE, ignore.case = TRUE)[1])
      updateSelectInput(session, "rdsEmbedding2", choices = numeric_cols, selected = grep("UMAP_2", numeric_cols, value = TRUE, ignore.case = TRUE)[1])
      common_counts <- c("nCount_RNA", "total_counts", "sum")
      updateSelectInput(session, "rdsTotalCounts", choices = numeric_cols, selected = intersect(common_counts, numeric_cols)[1])
      common_clusters <- c("seurat_clusters", "celltype", "cluster")
      updateSelectInput(session, "rdsDefaultCluster", choices = factor_cols, selected = intersect(common_clusters, factor_cols)[1])
      updateSelectInput(session, "rdsBarcode", choices = c("rownames", all_cols), selected = "rownames")
      
      # Update Assay/Slot selector based on non-empty matrix dimensions
      if (inherits(obj, "Seurat")) {
        assay_obj <- obj[[DefaultAssay(obj)]]
        available_slots <- c()
        
        # Check 'counts' slot
        counts_matrix <- try(assay_obj@counts, silent = TRUE)
        if (!inherits(counts_matrix, "try-error") && all(dim(counts_matrix) > 0)) {
          available_slots <- c(available_slots, "counts")
        }
        # Check 'data' slot
        data_matrix <- try(assay_obj@data, silent = TRUE)
        if (!inherits(data_matrix, "try-error") && all(dim(data_matrix) > 0)) {
          available_slots <- c(available_slots, "data")
        }
        # Check 'scale.data' slot
        scale_data_matrix <- try(assay_obj@scale.data, silent = TRUE)
        if (!inherits(scale_data_matrix, "try-error") && all(dim(scale_data_matrix) > 0)) {
          available_slots <- c(available_slots, "scale.data")
        }
        updateSelectInput(session, "rdsAssaySlot", choices = available_slots, selected = "data")
        
      } else { # Assumes SCE
        all_assays <- assayNames(obj)
        available_assays <- c()
        if (length(all_assays) > 0) {
          for (aname in all_assays) {
            assay_matrix <- try(assay(obj, aname), silent = TRUE)
            if (!inherits(assay_matrix, "try-error") && all(dim(assay_matrix) > 0)) {
              available_assays <- c(available_assays, aname)
            }
          }
        }
        updateSelectInput(session, "rdsAssaySlot", choices = available_assays, selected = "counts")
      }
      
      shinyjs::show("rds_mapping_div")
      showNotification("Object read successfully. Please map columns.", type = "message")
    }
    
    # ---- RDS Pre-processing Logic ----
    observeEvent(input$rdsFile, {
      req(input$rdsFile)
      shinyjs::hide("rds_mapping_div")
      data_out$config <- NULL # Clear previous data
      temp_data_source_name(input$rdsFile$name)
      tryCatch({
        obj <- readRDS(input$rdsFile$datapath)
        raw_rds_object(obj)
        process_loaded_object(obj)
      }, error = function(e) {
        showModal(modalDialog(title = "Error Reading RDS File", tags$pre(e$message), easyClose = TRUE))
      })
    })
    
    # ---- HDF5-backed SCE Loading Logic ----
    # Set up the directory chooser server logic
    volumes <- c(Home = fs::path_home(), "R_LIBS" = .libPaths()[1], getVolumes()())
    shinyDirChoose(input, "hdf5Dir", roots = volumes, session = session)
    
    # A reactive to hold the parsed path
    hdf5_path_r <- reactive({
      req(input$hdf5Dir)
      parseDirPath(volumes, input$hdf5Dir)
    })
    
    # Render the chosen path in the UI
    output$hdf5Path <- renderText({
      hdf5_path_r()
    })
    
    # Observe the path and load data when it's available
    observeEvent(hdf5_path_r(), {
      path <- hdf5_path_r()
      req(length(path) > 0) # Ensure a path is selected
      
      shinyjs::hide("rds_mapping_div")
      data_out$config <- NULL # Clear previous data
      temp_data_source_name(path)
      
      tryCatch({
        showNotification("Loading HDF5-backed SummarizedExperiment...", type = "message", duration = NULL, id = "hdf5_loading")
        obj <- HDF5Array::loadHDF5SummarizedExperiment(dir = path)
        removeNotification(id = "hdf5_loading")
        raw_rds_object(obj)
        process_loaded_object(obj)
      }, error = function(e) {
        removeNotification(id = "hdf5_loading")
        showModal(modalDialog(title = "Error Loading HDF5 Data", tags$pre(e$message), easyClose = TRUE))
      })
    }, ignoreInit = TRUE)
    
    # ---- RDS & HDF5 Final Processing Logic ----
    observeEvent(input$processRds, {
      req(raw_rds_object(), input$rdsEmbedding1, input$rdsEmbedding2, input$rdsBarcode, input$rdsTotalCounts, input$rdsDefaultCluster, input$rdsAssaySlot)
      
      obj <- raw_rds_object()
      tryCatch({
        
        # 1. Construct obs_data from object metadata and all reductions
        obs_df <- extract_full_metadata(obj)
        
        barcode_col_name <- if (input$rdsBarcode == "rownames") {
          obs_df[[".__barcode"]] <- rownames(obs_df)
          ".__barcode"
        } else {
          input$rdsBarcode
        }
        
        obs_data <- data.table::as.data.table(obs_df)
        
        # 2. Construct gene_table from object feature metadata (rowData or meta.features)
        gene_table <- if (inherits(obj, "Seurat")) {
          # For Seurat, feature metadata is in the 'meta.features' slot of an assay
          assay_name <- DefaultAssay(obj)
          if (!is.null(assay_name) && assay_name %in% names(obj@assays)) {
            data.table::as.data.table(obj[[assay_name]]@meta.features, keep.rownames = "var_names")
          } else {
            data.table::data.table() # Return empty table if no assay found
          }
        } else { # Assumes SCE
          # For SingleCellExperiment, use rowData()
          data.table::as.data.table(rowData(obj), keep.rownames = "var_names")
        }
        
        # Fallback: If no metadata was extracted, just use the gene names
        if (is.null(gene_table) || nrow(gene_table) == 0 || !"var_names" %in% names(gene_table)) {
          gene_names <- rownames(obj)
          gene_table <- data.table::data.table(var_names = gene_names)
          showNotification("No feature metadata found, using gene names only.", type = "warning", duration = 5)
        }
        
        # 3. Construct a dynamic config list
        config <- list(
          obs = list(columns = list(
            embedding1 = input$rdsEmbedding1,
            embedding2 = input$rdsEmbedding2,
            barcode = barcode_col_name,
            total_counts = input$rdsTotalCounts,
            celltype = input$rdsDefaultCluster
          )),
          quant = list(
            barcode = barcode_col_name,
            slot = input$rdsAssaySlot # Store the selected slot
          )
        )
        
        # 4. Create the unified gene data retrieval function
        get_gene_data_rds <- function(gene_name) {
          all_gene_names <- rownames(obj)
          if (!gene_name %in% all_gene_names) {
            warning(paste("Gene not found in object:", gene_name)); return(NULL)
          }
          
          if (inherits(obj, "Seurat")) {
            assay_to_use <- DefaultAssay(obj)
            # Use the slot selected by the user
            counts_vec <- GetAssayData(obj, slot = config$quant$slot, assay = assay_to_use)[gene_name, ]
          } else { # Assumes SCE
            # Use the assay selected by the user
            assay_to_use <- config$quant$slot
            if (!assay_to_use %in% assayNames(obj)) {
              stop(paste("Selected assay not found:", assay_to_use))
            }
            counts_vec <- assay(obj, assay_to_use)[gene_name, ]
          }
          
          # Create a data.table with cell barcodes and expression values
          data_list <- list(
            colnames(obj),
            as.numeric(counts_vec)
          )
          names(data_list) <- c(barcode_col_name, "expression")
          return(as.data.table(data_list))
        }
        
        # 5. Store all results
        data_out$config <- config
        data_out$obs_data <- obs_data
        data_out$gene_table <- gene_table
        data_out$get_gene_data <- get_gene_data_rds
        data_out$data_source_name <- temp_data_source_name()
        
        showNotification("Data loaded successfully from object.", type = "message", duration = 5)
        
      }, error = function(e) {
        showModal(modalDialog(title = "Error Processing Object", tags$pre(e$message), easyClose = TRUE))
      })
    })
    
    # ---- Render Data Previews ----
    output$obsDataPreview <- DT::renderDataTable({
      req(data_out$obs_data)
      DT::datatable(
        head(data_out$obs_data, 10),
        options = list(pageLength = 10, scrollX = TRUE, info = FALSE, searching = FALSE),
        rownames = FALSE,
        caption = "Preview of cell metadata (obs)."
      )
    })
    
    output$varDataPreview <- DT::renderDataTable({
      req(data_out$gene_table)
      DT::datatable(
        head(data_out$gene_table, 10),
        options = list(pageLength = 10, scrollX = TRUE, info = FALSE, searching = FALSE),
        rownames = FALSE,
        caption = "Preview of gene/variable metadata (var)."
      )
    })
    
    return(data_out)
  })
}