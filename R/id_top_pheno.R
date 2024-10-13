#' @title id_top_pheno
#'
#' @description id the level of x[, site] where the mean relative value of
#' x[, pheno] is highest. By default, uses shrinkage estimates of mean
#' relative x[, pheno] (via lme4) if any site occurs in at least two
#' years. Otherwise (or optionally) uses means.
#'
#' @param x data.frame or matrix for a single unit, ex a genotype. Coerced to data.frame.
#' @param site character of grouping column, ex field site locations
#' @param geno character, column name indicating genetic units
#' @param pheno character of response, ex yield, biomass
#' @param method character, method of finding highest relative performing location. Defaults to 'blup'. See details.
#' @param verbose logical (TRUE). Print messages?
#' @param ... other arguments passed to method.
#'
#' @details method is passed to tapply(), so must return a single value per genotype-site combination.
#'
#' @returns x, with an additional logical column "is_home", which indicates
#' whether the value of x[, site] is the home site, and a column "method_rel_pheno", the
#' data used to identify home site.
#'
#' @importFrom nlme lme
#' @importFrom stats fitted
#' 
#' @export

id_top_pheno <- function(x, site, geno, pheno, method=c('blup', 'mean', 'median'), verbose=TRUE, ...) {
  
  if (!(method %in% c('blup', 'mean', 'median', 'quantile'))) {
    warning("You are trying to identify home site using an untested method.
            Try sticking to 'blup', 'mean', or 'median'.")
  }
  
  # remove traces of old calculations
  is_old_est <- paste0(method, '_rel_') |>
    grepl(colnames(x))
  x <- x[, !is_old_est]
  
  missing_obs <- is.na(x[, pheno])
  x <- x[!missing_obs, ]
  
  x[, c(geno, site)] <- x[, c(geno, site)] |>
    lapply(factor) |>
    lapply(droplevels)
  
  # calculate mean pheno within site.
  if (method == 'blup') {

    est_col <- paste0('blup_', pheno)
    x[, 'genosite'] <- paste(x[, geno], x[, site], sep=':')
    
    site_eff =
      paste(pheno, 1, sep='~') |>
      formula() |>
      lme(data=x, random <- ~ 1 | genosite) |>
      fitted()
    x[, est_col] <- site_eff
    
    site_eff <- x[, c(geno, site, est_col)]
    
  } else {
    # Calculate
    est_col <- paste(method,
                     pheno,
                     sep='_')
    site_eff <- tapply(x[, pheno],
                       x[, c(geno, site)],
                       method,
                       ...)   # faster than aggregate
    
    # Munge
    sites <- sapply(colnames(site_eff), rep, nrow(site_eff)) |>
      c()
    genos <- rep(rownames(site_eff), ncol(site_eff))
    vals <- c(site_eff)
    site_eff <- data.frame(genos, sites, vals)
    colnames(site_eff) <- c(geno, site, est_col)
  }
  
  # identify site with max mean
  .max_site <- function(x) {
    x[, site] <- as.character(x[, site])
    out <- x[, est_col] |>
      which.max()
    out <- x[out, site]
    out
  }
  max_site <- tapply(site_eff,
                     site_eff[,geno],
                     .max_site) |>
    data.frame(max_site=_)
  
  # merge, preserving column names
  cn <- colnames(x)
  
  x <- merge(x, max_site, by.x=geno, by.y='row.names', sort=FALSE) #|>
  if (method != 'blup') {
    x <- merge(x, site_eff, by=c(geno, site), sort=FALSE)
  }
  
  x <- x[, c(cn, 'max_site')]
  
  x$is_home <- x[, site] == x[, 'max_site']
  x$max_site <- NULL
  x$genosite <- NULL
  
  return(x)
}