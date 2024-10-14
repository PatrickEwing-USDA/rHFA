#' @title id_home
#'
#' @description Identify the home site based on highest relative phenotype value.
#'
#' @param data data.frame of performance data by site, year, and variety, for example.
#' @param site character, column name indicating spatial units
#' @param year character, column name indicating temporal units
#' @param geno character, column name indicating genetic units
#' @param pheno character, column name indicating the phenotype (ex. yield)
#' @param method character, method of finding highest relative performing location. Defaults to 'blup'. See details.
#' @param verbose logical (TRUE). Print messages?
#' @param ... other arguments passed to `method`.
#'
#' @details The home site for a given variety is defined as the location where
#' variety performs best across years, relative to other varieties. It's calculated
#' by:
#'
#' 1. Calculate relative phenotype within site-year by scaling (mean = 0) and
#' scaling (st. dev. = 1) phenotype across varieties within site-years.
#' 2. Find the expected relative phenotype for each variety within
#' within sites across years.
#' 3. For each variety, identify the site with the highest expected relative phenotype value.
#'
#' @returns the data.frame `data` with three new columns:
#' - `'rel_<pheno>'`: numeric, the relative phenotype value of each variety within site-year
#' - `is_home`: logical of whether a site is a variety's home site.
#' - `<method>_rel_<pheno>`: numeric, the expected phenotype of each variety at that location
#' as derived by `method`.
#'
#' @seealso `tapply()`, `lme4::lmer()`, `id_top_pheno()`
#'
#' @export

id_home <- function(data, site, year, geno, pheno, method = c('blup', 'mean', 'median'), verbose = FALSE, ...) {
  
  rel_colname <- paste0('rel_', pheno)
  
  # Preserve row order
  if (is.null(rownames(data))) {
    rownames(data) <- seq_len(nrow(data))
  }
  rn <- rownames(data)
  data$rn <- rn
  
  # Center and scale performance within site-year
  site_year <- paste(data[, site], data[, year], sep = '_')
  
  .scale_pheno <- function(x, pheno) {
    # Create the relative phenotype column
    colname <- paste0('rel_', pheno)
    x[, colname] <- scale(x[, pheno])
    return(x)
  }
  
  # Apply scaling to each site-year group
  data_split <- split(data, site_year)
  data_list <- lapply(data_split, .scale_pheno, pheno)
  data <- do.call(rbind, data_list)
  
  # Find highest relative yield for each genotype using method
  data <- id_top_pheno(data, site, geno, rel_colname, method, verbose, ...)
  
  # Restore original rownames and order
  rownames(data) <- data$rn
  data$rn <- NULL
  data <- data[rn, ]
  
  return(data)
}
