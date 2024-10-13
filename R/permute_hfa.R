# See permute_hfa_internals.R for many of the .functions.

#' @title permute_hfa
#'
#' @description test the magnitude and significance of the home field advantage
#' versus what would be expected by chance.
#'
#' @param data data.frame
#' @param level parameter at which to apply the home field advantage. Options are
#' 'population' (default), 'genotype', 'site', and 'year'.
#' @param site column name containing spatial environmental information
#' @param year column name containing temporal environmental information
#' @param geno column name containing genotype or variety information
#' @param pheno column name containing phenotype or performance information
#' @param popn column name differentiating sub-populations of genotypes within the data set.
#' @param times how many times do you want to permute (in addition to observed data)?
#' @param method *character* for estimating relative performance at each location. Default is
#' 'blup'. Other arguments that name functions (ex. `median`, `mean`, `quantile`) can be used.
#' @param parallel *logical* use parallel processing? Currently only on UNIX.
#' @param seed *numeric* for reproducible permutations
#' @param ... added arguments passed to `method`.
#'
#' @details This function is the only one that's needed for an HFA analysis. It identifies a home
#' site and calculates the home field advantage metric at the level chosen based on observed data.
#' It then permutes data within site-year and repeats this process to estimate the difference
#' between observed and expected (permuted) HFA - and tests whether this difference is zero.
#'
#' Note that observed HFA may differ slightly from "manual" calculations (ex. via `lm()`)
#' as it uses sparse matrix operations for efficiency. See `Matrix::qr()` and `Matrix::qrR()`.
#'
#' @return A list with two items: `'home_field'` is a `data.frame` that
#' summarizes the results. `perms` is a matrix (or list of matrices) of permutations.
#'
#' @seealso `id_home()`
#' @importFrom parallel detectCores mclapply
#' @importFrom stats setNames na.omit
#' 
#' @export

permute_hfa <- function(data,
                        level = c('population', 'genotype', 'year', 'site'),
                        site = NA,
                        year = NA,
                        geno = NA,
                        pheno = NA,
                        popn = NA,
                        times = 99,
                        method = 'blup',
                        parallel = TRUE,
                        seed = NULL,
                        ...) {
  
  ncpu <- ifelse(parallel, detectCores(), 1)
  level <- match.arg(level)
  
  # new column names
  rel_pheno <- paste0('rel_', pheno)
  
  # formula for calculating HFA effect
  ff <- switch(level,
               'population' = ' ~ site*year + geno + is_home',  # overall HFA
               'genotype' = ' ~ site*year + geno + geno:is_home',
               'year' = '~ site*year + geno + year:is_home',  # is this right?
               'site' = '~ year*site + geno + site:is_home')  # is this right? HFA for each genotype
  
  ff <- formula(gsub('geno', geno, gsub('year', year, gsub('site', site, ff))))
  
  # select and format data into list of dataframes for each population
  #   adjusting the column selection inside permute_hfa
  dd <- if (is.na(popn)) {
    na.omit(data[, c(site, year, geno, pheno)])
  } else {
    na.omit(data[, c(site, year, geno, pheno, popn)])
  }
  
  if (is.na(popn)) {
    dd <- list(dd)
  } else {
    dd <- split(dd, dd[, popn])
  }
  
  # calculate relative yields
  dd <- mclapply(dd, function(x) {
    # calc rel_pheno and is_home
    x <- id_home(x, site, year, geno, pheno, method, verbose = FALSE, ...)
    return(x)
  }, mc.cores = ncpu)
  
  # permute HFA within each population
  results <- lapply(dd, function(x) {
    # Set up structured permutations
    sets <- generate_sets(x, site, year, times, seed)
    
    # Permute HFA
    coef_permute <- do.call(cbind, mclapply(sets, function(ss) {
      # ID home site
      x[, c(pheno, rel_pheno)] <- x[ss, c(pheno, rel_pheno)]  # permute phenotypes within site-year
      x <- id_top_pheno(x, site = site, geno = geno, pheno = rel_pheno, method = method, verbose = FALSE, ...)
      
      # calculate HFA using Matrix and qr decomposition
      home_coef <- calculate_hfa(x, ff, pheno, geno, site, year)
      return(home_coef)
    }, mc.cores = ncpu))
    
    colnames(coef_permute) <- c('observed', paste0('perm', 1:(ncol(coef_permute) - 1)))
    
    # Calculate p-values and effects
    test <- cbind(
      intervals = calculate_intervals(coef_permute),
      p_val = .two_tailed(coef_permute)
    )
    
    results <- list(results = test, perms = coef_permute)
    return(results)
  })
  
  # format results
  if (is.null(names(results))) {
    names(results) <- paste0('popn', seq_len(length(results)))
  }
  
  if (level == 'population') {
    test_results <- rownames(do.call(rbind, lapply(results, function(x) x$results)), names(results))
    perms <- rownames(do.call(rbind, lapply(results, function(x) x$perms)), names(results))
  } else {
    lvl <- switch(level, 'genotype' = geno, 'year' = year, 'site' = site)
    test_results <- do.call(rbind, lapply(names(results), function(x) {
      res <- data.frame(popn = x, lvl = rownames(results[[x]]$results), results[[x]]$results, stringsAsFactors = FALSE)
      colnames(res) <- gsub('lvl', lvl, colnames(res))
      return(res)
    }))
    perms <- lapply(results, function(x) x$perms)
  }
  
  # finish formula to include the phenotype
  ff <- formula(paste(pheno, '~', as.character(ff)[2]))
  
  out <- list(home_field = test_results, data = do.call(rbind, dd), perms = perms, formula = ff)
  attr(out, 'Print_Warning') <- "A print method is on its way. For now, you probably want to call obj$home_field."
  
  return(out)
}
