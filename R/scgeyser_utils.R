#' HELPER FUNCTIONS | R/utils.R ----
#' Define Core Plotting Theme
#'
#' Creates a consistent ggplot theme for all plots to reduce code duplication.
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 theme element_blank element_text unit element_rect
#' @return A list of ggplot theme elements.
#' 
get_core_plot_theme <- function() {
  list(
    cowplot::theme_cowplot(font_size = 7), 
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 2),
      legend.key.width = unit(0.07, 'cm'),
      legend.key.height = unit(0.15, 'cm'),
      plot.title = element_text(hjust = 0.5, size = 4, face = "bold"),
      axis.title = element_text(size = 4), 
      strip.background = element_rect(fill="grey90", linetype = "blank"), 
      strip.text = element_text(face="bold", size = 4),
      plot.caption = element_text(size = 3),
      legend.position = c(0, 0)
    )
  )
}

#' @title Summarize Grouped Data
#' @description Calculates counts and ratios for nested groups within a dataframe.
#' @param data The input data.frame or data.table.
#' @param ... One or more unquoted column names to group by. The ratio is
#'   calculated relative to the count of the first grouping variable.
#' @param threshold A numeric value (0-1). Rows with a `Ratio` below this
#'   value will be filtered out. Set to `NULL` to disable.
#' @importFrom rlang ensyms
#' @importFrom dplyr group_by summarise n left_join mutate arrange desc filter
#' @return A summarized and arranged tibble.

summarize_grouped_data <- function(data, ..., threshold = 0.01) {
  Sum <- Ratio <- Count <- NULL
  groupings <- ensyms(...)
  missing <- setdiff(as.character(groupings), colnames(data))
  if (length(missing) >= 1) {
    stop(paste(missing, "is not a column name."))
  }
  
  if (length(groupings) < 1) {
    stop("At least one grouping variable must be provided.")
  }
  
  sum_vals <- data %>% 
    group_by(!!groupings[[1]]) %>% 
    summarise(Sum = n(), .groups = 'drop')
  
  result <- data %>% 
    group_by(!!!groupings) %>% 
    summarise(Count = n(), .groups = 'drop') %>%
    left_join(sum_vals, by = as.character(groupings[[1]])) %>%
    mutate(Ratio = Count / Sum) %>%
    arrange(!!groupings[[1]], desc(Ratio))
  
  if (!is.null(threshold) && threshold > 0) {
    result <- result %>% filter(Ratio > threshold)
  }
  
  return(result)
}