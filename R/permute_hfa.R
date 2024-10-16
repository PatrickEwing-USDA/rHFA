#' @title .quick_resid, .quick_coef <internal>
#'
#' @description calculate residuals (coefficients) of a linear model using Matrix, about an
#' order of magnitude faster than residuals(lm.fit()) for decently large, highly
#' sparse model matrices.
#'
#' @param ff formula
#' @param data data.frame containing formula
#'
#' @return A vector of residuals (coefficients).
#'
#' @importFrom Matrix Matrix qr qr.resid qr.coef
#' @importFrom stats model.frame model.matrix

.quick_resid = function(ff, data) {
  
  mf = model.frame(ff, data)
  X = model.matrix(ff, mf) |>
    Matrix()
  mX = qr(X) |>
    suppressWarnings()
  Y = mf[, 1]
  
  out = qr.resid(mX, Y) |>
    suppressWarnings()
  
  return(out)
}

.quick_coef = function(ff, data) {
  
  mf = model.frame(ff, data)
  X = model.matrix(ff, mf) |>
    Matrix()
  mX = qr(X) |>
    suppressWarnings()
  Y = mf[, 1]
  
  out = qr.coef(mX, Y) |>
    suppressWarnings()
  
  return(out)
}


#' @title .generate_sets
#'
#' @description generate permutation sets structured by site-year
#'
#' @param x data.frame
#' @param site column identifying spatial environment
#' @param year column identifying temporal environment
#' @param times number of sets to produce.
#' @param seed optional random seed for reproducibility.
#'
#' @return a data.frame of permutation sets - row numbers used to reorder x. The first column is the original data.
#'
#' @importFrom permute shuffleSet how

.generate_sets <- function(x, site, year, times, seed = NULL) {
  
  control <- as.factor(paste(x[, site], x[, year], sep = '_'))
  N <- nrow(x)
  
  set.seed(seed)
  ss <- shuffleSet(N, times, control = how(blocks = control)) # permutations
  ss <- rbind(seq_len(N), ss) # include original data
  ss <- as.data.frame(t(ss))
  
  return(ss)
}


#' @title .calculate_hfa
#'
#' @description calculate the home field advantage based on the formula, ff. Uses efficient
#' packages for fastest permutation results.
#'
#' @param x data.frame
#' @param ff formula to use. part_pheno ~ geno */+ is_home
#' @param site character, column name indicating spatial units
#' @param year character, column name indicating temporal units
#' @param geno character, column name indicating genetic units
#' @param pheno character, column name indicating the phenotype (ex. yield)
#'
#' @return a data.frame of home_coefficients
#' 
#' @importFrom stats formula

.calculate_hfa <- function(x,
                           ff,
                           pheno,
                           geno = NA,
                           site = NA,
                           year = NA) {
  # Convert formula to character and reconstruct
  ff <- as.character(ff)
  ff <- formula(paste(pheno, ff[2], sep = '~'))
  
  # Calculate coefficients using the formula
  home_coef <- .quick_coef(ff, x)
  
  # Extract is_home coefficients
  is_home <- grepl('is_homeTRUE', names(home_coef))
  home_coef <- home_coef[is_home]
  
  # Format the coefficient names for better readability
  names(home_coef) <- gsub(':is_homeTRUE', '', names(home_coef))
  names(home_coef) <- gsub(geno, '', names(home_coef))
  names(home_coef) <- gsub(year, '', names(home_coef))
  names(home_coef) <- gsub(site, '', names(home_coef))
  names(home_coef) <- gsub('TRUE', '', names(home_coef))
  
  return(home_coef)
}

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
    sets <- .generate_sets(x, site, year, times, seed)
    
    # Permute HFA
    coef_permute <- do.call(cbind, mclapply(sets, function(ss) {
      # ID home site
      x[, c(pheno, rel_pheno)] <- x[ss, c(pheno, rel_pheno)]  # permute phenotypes within site-year
      x <- id_top_pheno(x, site = site, geno = geno, pheno = rel_pheno, method = method, verbose = FALSE, ...)
      
      # calculate HFA using Matrix and qr decomposition
      home_coef <- .calculate_hfa(x, ff, pheno, geno, site, year)
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
