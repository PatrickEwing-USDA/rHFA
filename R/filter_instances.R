#' @title count_instances
#'
#' @description Count the number of times each category of `group_col` appears
#' per `rep_col`.
#'
#' @param data A *data.frame*.
#' @param group_col A *character* or *numeric* column in `data` for primary grouping.
#' @param rep_col A *character* or *numeric* column in `data` indicating replicates of
#' `group_col`.
#' @param na.rm A *logical* value. If TRUE, remove rows containing missing values.
#'
#' @return A *vector* with the names of levels in `group_col` and the unique counts.
#'
#' @seealso `filter_instances()`, `autofilter_instances()`
#' 
#' @importFrom stats na.omit
#' @export

count_instances <- function(data, group_col, rep_col, na.rm = FALSE) {
  if (na.rm) {
    data <- droplevels(na.omit(data))
  }
  unique_data <- unique(data[, c(group_col, rep_col), drop = FALSE])
  instance_counts <- tapply(unique_data[[rep_col]], unique_data[[group_col]], length)
  return(instance_counts)
}

#' @title filter_instances
#'
#' @description Keep categories of `group_col` that occur in at least
#' `min_times` of `rep_col` categories.
#'
#' @param data A *data.frame*.
#' @param group_col A *character* or *numeric* column in `data` for primary grouping.
#' @param rep_col A *character* or *numeric* column in `data` indicating replicates of
#' `group_col`.
#' @param min_times A *numeric* value for the minimum number of instances, e.g., minimum site-years.
#' @param na.rm A *logical* value. If TRUE, remove rows containing missing values.
#'
#' @return A *data.frame*, a filtered version of `data`.
#'
#' @seealso `count_instances()`, `autofilter_instances()`
#' @export
filter_instances <- function(data, group_col, rep_col, min_times = 2, na.rm = FALSE) {
  instance_counts <- count_instances(data, group_col, rep_col, na.rm)
  valid_groups <- names(instance_counts[instance_counts >= min_times])
  filtered_data <- droplevels(data[data[[group_col]] %in% valid_groups, , drop = FALSE])
  return(filtered_data)
}

#' @title autofilter_instances
#'
#' @description Filters geno-years, site-years, and genos per site-year to
#' `min_times` instances.
#'
#' @param data A *data.frame*.
#' @param site A *character* or *numeric* column in `data` representing site.
#' @param year A *character* or *numeric* column in `data` representing year.
#' @param geno A *character* or *numeric* column in `data` representing genotype.
#' @param min_times A *numeric* value for the minimum instances, e.g., minimum site-years.
#' @param na.rm A *logical* value. If TRUE, remove rows containing missing values.
#' @param max_cycles A *numeric* value for the maximum number of iterations for filtering.
#'
#' @return A *data.frame*, a filtered version of `data`.
#' @export
autofilter_instances <- function(data, site, year, geno, min_times = 2, na.rm = TRUE, max_cycles = 999) {
  rn <- rownames(data)
  if (is.null(rn)) {
    rn <- seq_len(nrow(data))
    rownames(data) <- rn
  }
  data$siteyear <- paste(data[[site]], data[[year]], sep = "___")
  prev_nr <- nrow(data)
  cycles <- 0
  repeat {
    data <- data |>
      filter_instances(geno, year, min_times, na.rm) |>
      filter_instances(site, year, min_times, na.rm) |>
      filter_instances("siteyear", geno, min_times, na.rm)
    new_nr <- nrow(data)
    if (new_nr == prev_nr || cycles >= max_cycles) break
    prev_nr <- new_nr
    cycles <- cycles + 1
  }
  keep_rn <- rn %in% rownames(data)
  out <- data[keep_rn, , drop = FALSE]
  if (cycles >= max_cycles) {
    message("Stopped filtering after ", cycles, " cycles. Consider increasing `max_cycles`.")
  } else {
    message("Filtered geno-years, site-years, and genos per site-year to at least ", min_times, " instances over ", cycles, ifelse(cycles > 1, " cycles.", " cycle."))
  }
  return(out)
}
