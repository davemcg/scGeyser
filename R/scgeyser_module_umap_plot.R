#' UI for UMAP Plot Controls
#' @description Creates the UI for the UMAP plot controls in the sidebar.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList h4 h5 numericInput checkboxInput selectInput hr downloadButton
#'
umapPlotControlsUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("2. Plotting Options"),
    
    h5("Shared Settings"),
    numericInput(ns("pointSize"), "Base Point Size:", value = 2.0, min = 0.1, max = 10, step = 0.5),
    
    # New Subsampling Section
    h5("Subsampling"),
    checkboxInput(ns("enableSubsample"), "Enable Subsampling", value = TRUE),
    numericInput(ns("subsampleSize"), "Number of cells:", value = 100000, min = 100, step = 100),
    selectInput(ns("subsampleColumn"), "Stratify by:", choices = NULL),
    
    hr(),
    
    h5("Download"),
    downloadButton(ns("downloadGenePlot"), "Gene Plot"),
    downloadButton(ns("downloadMetadataPlot"), "Metadata Plot")
  )
}

#' UI for UMAP Plot Output Area
#' @description Creates the UI for the UMAP plot output in the main panel.
#' @param id Namespace ID.
#' @importFrom shiny NS tagList fluidRow column h4 plotOutput hr tags selectInput numericInput checkboxInput actionButton icon
#' @importFrom DT dataTableOutput
#'
umapPlotOutputUI <- function(id) {
  ns <- NS(id)
  tagList(
    # --- ROW 1: PLOTS ---
    fluidRow(
      column(
        width = 6,
        h4("Gene Expression"),
        plotOutput(ns("genePlot"), height = "600px")
      ),
      column(
        width = 6,
        h4("Categorical Metadata"),
        # The 'click' argument enables click events for this plot
        plotOutput(ns("metadataPlot"), height = "600px", click = ns("plot_click"))
      )
    ),
    
    hr(),
    
    # --- ROW 2: CONTROLS ---
    fluidRow(
      # Left controls
      column(
        width = 6,
        tags$strong("Selected Gene:"),
        verbatimTextOutput(ns("selectedGeneDisplay"), placeholder = TRUE)
      ),
      # Right controls
      column(
        width = 6,
        fluidRow(
          column(7,
                 selectInput(ns("metadataColumn"), "Select Metadata Column", choices = NULL)
          ),
          column(5,
                 # This input allows the user to change the number of nearest cells to find
                 numericInput(ns("numNeighbors"), "Nearest Cells (click)", value = 5, min = 1, max = 50)
          )
        ),
        fluidRow(column(3, checkboxInput(ns("focusOnClick"), "Click-to-focus", value = TRUE)),
                 column(3,                 
                        checkboxInput(ns("showLabels"), "Show category labels", value = TRUE))
        )
      )
    ),
    
    # --- ROW 3: ALIGNED BUTTONS ---
    fluidRow(
      column(6,
             actionButton(ns("generateGenePlot"), "Generate Gene Plot", class = "btn-primary w-100", icon = icon("paint-brush"))
      ),
      column(6,
             actionButton(ns("generateMetaPlot"), "Generate Metadata Plot", class = "btn-primary w-100", icon = icon("paint-brush"))
      )
    ),
    
    hr(),
    
    # --- ROW 4: TABLES ---
    fluidRow(
      column(
        width = 6,
        h4("Gene Table"),
        geneTableUI(ns("gene_table_module"))
      ),
      column(
        width = 6,
        h4("Metadata for Nearest Cells (click plot to update)"),
        # This table will display the metadata for the clicked cells
        DT::dataTableOutput(ns("clickInfoTable"))
      )
    )
  )
}


#' Server for the UMAP Plot Module
#' @description Handles the logic for generating and displaying UMAP plots for gene expression and metadata.
#' @param id Namespace ID.
#' @param loaded_data A reactive list from the dataLoaderServer module.
#' @importFrom shiny moduleServer NS reactiveVal reactive eventReactive observeEvent req updateSelectInput isolate nearPoints renderText validate need renderPlot downloadHandler
#' @importFrom DT renderDataTable datatable
#' @importFrom stats median runif
#' @importFrom dplyr group_by summarise slice_sample ungroup bind_rows filter mutate left_join
#' @importFrom data.table copy fifelse
#' @importFrom pals trubetskoy okabe alphabet2 alphabet glasbey viridis
#' @importFrom scattermore geom_scattermore
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot annotate aes scale_color_manual scale_alpha_identity geom_point guides ggtitle scale_colour_gradient2 scale_colour_gradientn labs
#'
umapPlotServer <- function(id, loaded_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #-- Reactive State ----
    clicked_cells_data <- reactiveVal(NULL)
    focused_category <- reactiveVal(NULL)
    highlighted_labels_data <- reactiveVal(NULL)
    redraw_meta_trigger <- reactiveVal()
    redraw_gene_trigger <- reactiveVal()
    gene_to_plot <- reactiveVal(NULL) # Holds the gene selected from the table
    
    #-- Module Wiring ----
    selected_gene_from_table <- geneTableServer("gene_table_module", reactive(loaded_data$gene_table))
    
    # When a gene is clicked in the table, update our reactive value
    observeEvent(selected_gene_from_table(), {
      gene_to_plot(selected_gene_from_table())
    }, ignoreNULL = TRUE)
    
    #-- UI Updates & State Management ----
    observeEvent(loaded_data$obs_data, {
      req(loaded_data$obs_data, loaded_data$config, loaded_data$gene_table)
      
      categorical_cols <- names(loaded_data$obs_data)[sapply(loaded_data$obs_data, function(x) is.character(x) || is.factor(x))]
      selected_col <- if (!is.null(loaded_data$config$obs$columns$celltype) && loaded_data$config$obs$columns$celltype %in% categorical_cols) {
        loaded_data$config$obs$columns$celltype
      } else {
        categorical_cols[1]
      }
      updateSelectInput(session, "metadataColumn", choices = categorical_cols, selected = selected_col)
      # Update the new subsampling column selector as well
      updateSelectInput(session, "subsampleColumn", choices = categorical_cols, selected = selected_col)
      
      # Clear any selected gene when new data loads
      gene_to_plot(NULL)
    })
    
    # Button triggers
    observeEvent(input$generateMetaPlot, {
      focused_category(NULL)
      clicked_cells_data(NULL)
      highlighted_labels_data(NULL)
      redraw_meta_trigger(runif(1))
    })
    
    observeEvent(input$generateGenePlot, {
      redraw_gene_trigger(runif(1))
    })
    
    #-- Core Reactive Logic ----
    
    # Create a new reactive for the (potentially) subsampled data
    subsampled_obs_data_r <- reactive({
      req(loaded_data$obs_data)
      
      # If subsampling is disabled, return the full dataset
      if (!isTRUE(input$enableSubsample)) {
        return(loaded_data$obs_data)
      }
      
      # If enabled, require the necessary inputs
      req(input$subsampleSize, input$subsampleColumn)
      
      df <- loaded_data$obs_data
      total_rows <- nrow(df)
      sample_size <- as.integer(input$subsampleSize)
      
      # If the total number of rows is less than or equal to the desired sample size, return the whole dataframe.
      if (total_rows <= sample_size) {
        return(df)
      }
      
      # Perform proportional stratified sampling
      set.seed(123) # for reproducibility
      df %>%
        dplyr::group_by(.data[[input$subsampleColumn]]) %>%
        dplyr::slice_sample(prop = sample_size / total_rows) %>%
        dplyr::ungroup() |> 
        data.table::data.table()
    })
    
    labels_df_r <- eventReactive(c(subsampled_obs_data_r(), input$metadataColumn), {
      req(subsampled_obs_data_r(), input$metadataColumn, loaded_data$config)
      embed1_col <- loaded_data$config$obs$columns$embedding1
      embed2_col <- loaded_data$config$obs$columns$embedding2
      metadata_col <- input$metadataColumn
      subsampled_obs_data_r() %>%
        group_by(.data[[metadata_col]]) %>%
        dplyr::summarise(
          !!embed1_col := median(.data[[embed1_col]]),
          !!embed2_col := median(.data[[embed2_col]])
        )
    })
    
    metadata_plot_r <- eventReactive(redraw_meta_trigger(), {
      req(loaded_data$config, subsampled_obs_data_r(), input$metadataColumn)
      core_plot_theme <- get_core_plot_theme()
      embed1_col <- loaded_data$config$obs$columns$embedding1
      embed2_col <- loaded_data$config$obs$columns$embedding2
      metadata_col <- input$metadataColumn
      
      p_meta <- ggplot(subsampled_obs_data_r(), aes(x = .data[[embed1_col]], y = .data[[embed2_col]], color = .data[[metadata_col]])) +
        scale_color_manual(values = unname(rep(c(pals::trubetskoy(), pals::okabe(), pals::alphabet2(), pals::alphabet(), pals::glasbey()), 100)))
      
      plot_df <- copy(subsampled_obs_data_r())
      current_focus <- focused_category()
      
      if (!is.null(current_focus) && isTRUE(input$focusOnClick)) {
        num_points <- nrow(subsampled_obs_data_r())
        background_alpha <- min(1.0, 10000 / num_points)
        plot_df[, plot_alpha := fifelse(get(metadata_col) %in% current_focus, 1.0, background_alpha)]
      } else {
        plot_df[, plot_alpha := 1.0]
      }
      
      p_meta <- p_meta +
        scattermore::geom_scattermore(data = plot_df[order(plot_alpha)], aes(alpha = plot_alpha), pointsize = as.numeric(input$pointSize), pixels = c(800,800)) +
        scale_alpha_identity()
      
      if (isTRUE(input$showLabels)) {
        all_labels_df <- labels_df_r()
        highlighted_labels_df <- highlighted_labels_data()
        if (!is.null(highlighted_labels_df)) {
          highlighted_names <- highlighted_labels_df[[metadata_col]]
          regular_labels_df <- all_labels_df %>% filter(!.data[[metadata_col]] %in% highlighted_names)
        } else { regular_labels_df <- all_labels_df }
        p_meta <- p_meta + ggrepel::geom_text_repel(data = regular_labels_df, aes(label = .data[[metadata_col]]), bg.color = 'white', max.overlaps = Inf, box.padding = 0.1, segment.size = 0.1, size = 1)
      }
      
      if (!is.null(clicked_cells_data()) && isTRUE(input$focusOnClick)) {
        median_point_data <- clicked_cells_data() %>% dplyr::summarise(!!embed1_col := median(.data[[embed1_col]], na.rm = TRUE), !!embed2_col := median(.data[[embed2_col]], na.rm = TRUE))
        p_meta <- p_meta + geom_point(data = median_point_data, aes(x = .data[[embed1_col]], y = .data[[embed2_col]]), color = "white",  size = 0.4, shape = 4, stroke = 0.8) +
          geom_point(data = median_point_data, aes(x = .data[[embed1_col]], y = .data[[embed2_col]]), color = "black",  size = 0.4, shape = 4, stroke = 0.4)
      }
      
      if (!is.null(highlighted_labels_data()) && isTRUE(input$focusOnClick)) {
        p_meta <- p_meta + ggrepel::geom_text_repel(data = highlighted_labels_data(), aes(label = .data[[metadata_col]]), bg.color = 'white', bg.r = 0.15, fontface = 'bold', max.overlaps = Inf, box.padding = 0.2, segment.size = 0.2, size = 1.5)
      }
      
      caption_text <- paste("scGeyser", loaded_data$data_source_name)
      p_meta <- p_meta + guides(color = "none") + ggtitle(metadata_col) + core_plot_theme + labs(caption = caption_text)
      return(p_meta)
    })
    
    gene_plot_r <- eventReactive(redraw_gene_trigger(), {
      req(loaded_data$config, loaded_data$obs_data)
      core_plot_theme <- get_core_plot_theme()
      
      selected_gene_val <- gene_to_plot()
      
      if (is.null(selected_gene_val) || selected_gene_val == "") {
        return(ggplot() + annotate("text", x=0, y=0, label="Select a gene from the table and click 'Generate Gene Plot'.", size=3) + cowplot::theme_cowplot())
      }
      
      gene_data <- loaded_data$get_gene_data(selected_gene_val)
      
      if (is.null(gene_data)) {
        return(ggplot() + annotate("text", x=0, y=0, label=paste("Could not load data for:", selected_gene_val), size=3) + cowplot::theme_cowplot())
      }
      
      # Join gene data with the FULL observation data first
      full_gene_data <- loaded_data$obs_data %>%
        left_join(gene_data, by = setNames(loaded_data$config$quant$barcode, loaded_data$config$obs$columns$barcode))
      
      # Handle different expression data types
      if (!is.null(loaded_data$config$quant$slot) && loaded_data$config$quant$slot == "counts") {
        full_gene_data <- full_gene_data %>%
          mutate(scount = (expression / .data[[loaded_data$config$obs$columns$total_counts]]) * 1e4) %>%
          mutate(gene_count = ifelse(is.na(scount), 0, scount))
      } else { # For 'data', 'scale.data', or MGIF data
        full_gene_data <- full_gene_data %>%
          mutate(gene_count = ifelse(is.na(expression), 0, expression))
      }
      
      # Apply special subsampling: keep all expressing cells, subsample non-expressing
      plot_data <- if (!isTRUE(input$enableSubsample)) {
        full_gene_data
      } else {
        expressing_cells <- full_gene_data %>% filter(gene_count > 0)
        non_expressing_cells <- full_gene_data %>% filter(gene_count == 0)
        
        n_expressing <- nrow(expressing_cells)
        target_total_size <- as.integer(input$subsampleSize)
        n_non_expressing_to_sample <- max(50000,
                                          max(0, target_total_size - n_expressing)) # ensure you have at least 50000
        
        subsampled_non_expressing <- if (n_non_expressing_to_sample < nrow(non_expressing_cells)) {
          set.seed(123) # for reproducibility
          non_expressing_cells %>% dplyr::slice_sample(n = n_non_expressing_to_sample)
        } else {
          non_expressing_cells
        }
        
        dplyr::bind_rows(subsampled_non_expressing, expressing_cells)
      }
      
      embed1_col <- loaded_data$config$obs$columns$embedding1
      embed2_col <- loaded_data$config$obs$columns$embedding2
      
      # Re-filter for expressing cells from the final `plot_data` to draw them on top
      expressing_cells <- plot_data %>% filter(gene_count > 0)
      
      # --- Dynamic Plot Aesthetics ---
      # Set up color scale and legend title based on data type
      slot_type <- loaded_data$config$quant$slot
      if (!is.null(slot_type) && slot_type == "scale.data") {
        # Handle scaled data which can be negative/positive
        max_abs_val <- max(abs(expressing_cells$gene_count), na.rm = TRUE)
        if (is.infinite(max_abs_val)) max_abs_val <- 1 # Fallback
        
        color_aes <- aes(color = gene_count)
        color_scale <- scale_colour_gradient2(
          low = "blue", mid = "grey80", high = "red",
          midpoint = 0, limits = c(-max_abs_val, max_abs_val),
          name = "Scaled Expression"
        )
      } else {
        # Handle counts or log-normalized data (non-negative)
        color_aes <- aes(color = log2(gene_count + 1))
        color_scale <- scale_colour_gradientn(
          colours = pals::viridis(256),
          name = "log2(Expression + 1)"
        )
      }
      
      caption_text <- paste("scGeyser", loaded_data$data_source_name)
      p_gene <- ggplot(plot_data, aes(x = .data[[embed1_col]], y = .data[[embed2_col]])) +
        scattermore::geom_scattermore(color = 'grey80', pointsize = as.numeric(input$pointSize), pixels = c(800,800)) +
        scattermore::geom_scattermore(data = expressing_cells, pointsize = as.numeric(input$pointSize), mapping = color_aes, pixels = c(800,800)) +
        color_scale +
        ggtitle(selected_gene_val) +
        core_plot_theme +
        labs(caption = caption_text)
      
      return(p_gene)
    })
    
    #-- Plot Interaction (Click) ----
    observeEvent(input$plot_click, {
      req(subsampled_obs_data_r(), loaded_data$config, input$numNeighbors, input$metadataColumn, input$plot_click)
      
      embed1_col <- loaded_data$config$obs$columns$embedding1
      embed2_col <- loaded_data$config$obs$columns$embedding2
      
      # Find nearest cells in the subsampled data that is currently plotted
      nearest_cells <- nearPoints(
        df = subsampled_obs_data_r(),
        coordinfo = input$plot_click,
        xvar = embed1_col,
        yvar = embed2_col,
        maxpoints = input$numNeighbors,
        threshold = Inf
      )
      
      clicked_cells_data(nearest_cells)
      
      if (isTRUE(input$focusOnClick)) {
        metadata_col <- isolate(input$metadataColumn)
        categories_to_highlight <- unique(nearest_cells[[metadata_col]])
        
        labels_to_show <- labels_df_r() %>% filter(.data[[metadata_col]] %in% categories_to_highlight)
        highlighted_labels_data(labels_to_show)
        focused_category(categories_to_highlight)
        
        redraw_meta_trigger(runif(1))
      }
    })
    
    #-- Render Outputs ----
    output$selectedGeneDisplay <- renderText({
      validate(need(gene_to_plot(), "Select gene in table below"))
      gene_to_plot()
    })
    
    output$metadataPlot <- renderPlot({ metadata_plot_r() }, res = 300)
    output$genePlot <- renderPlot({ gene_plot_r() }, res = 300)
    
    output$clickInfoTable <- DT::renderDataTable({
      req(clicked_cells_data())
      DT::datatable(
        clicked_cells_data(), 
        options = list(pageLength = 25, scrollX = TRUE), 
        rownames = FALSE, 
        selection = 'none'
      )
    })
    
    #-- Download Handlers ----
    output$downloadMetadataPlot <- downloadHandler(
      filename = function() {
        paste0("umap_plot_metadata_", input$metadataColumn, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggsave(file, plot = metadata_plot_r(), width = 10, height = 10, dpi = 300)
      }
    )
    
    output$downloadGenePlot <- downloadHandler(
      filename = function() {
        req(gene_to_plot())
        paste0("umap_plot_gene_", gene_to_plot(), "_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggsave(file, plot = gene_plot_r(), width = 10, height = 10, dpi = 300)
      }
    )
  })
}