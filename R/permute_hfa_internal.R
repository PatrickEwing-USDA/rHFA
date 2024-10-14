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


#' @title generate_sets
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
#' @export

generate_sets <- function(x, site, year, times, seed = NULL) {
  
  control <- as.factor(paste(x[, site], x[, year], sep = '_'))
  N <- nrow(x)
  
  set.seed(seed)
  ss <- shuffleSet(N, times, control = how(blocks = control)) # permutations
  ss <- rbind(seq_len(N), ss) # include original data
  ss <- as.data.frame(t(ss))
  
  return(ss)
}


#' @title calculate_hfa
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
#' @export

calculate_hfa <- function(x,
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
