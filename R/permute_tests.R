#' @title calculate_intervals
#'
#' @description calculate median and 90% confidence intervals of difference of permutation and observed data
#'
#' @param x permutation results. Rows are values, columns are permutations. The first column is observed.
#'
#' @return a matrix of columns "median", "q05" (5th percentile), "q95" (95th percentile) of (observed - permutation)
#' 
#' @importFrom stats median quantile
#' @export

calculate_intervals <- function(x) {
  
  # Subtract observed (first column) from each permutation result and multiply by -1
  difference <- sweep(x, 1, x[, 1, drop = FALSE], '-') * -1
  
  # Calculate the observed value, median, 5th percentile, and 95th percentile
  out <- cbind(
    HFA = x[, 1, drop=TRUE],
    expected = apply(x, 1, median, na.rm=TRUE),
    diff_q50 = apply(difference, 1, median, na.rm = TRUE),
    diff_q05 = apply(difference, 1, quantile, 0.05, na.rm = TRUE),
    diff_q95 = apply(difference, 1, quantile, 0.95, na.rm = TRUE)
  )
  
  return(out)
}


#' @title .two_tailed <internal>
#'
#' @description perform two-tailed test on permutation values
#'
#' @param x matrix. Each column is the result of a different permutation. The first column is the original data.
#'
#' @return a vector of p-values for each row in x.

.two_tailed <- function(x) {
  
  # Perform the two-tailed test: count how many times each permutation is greater than or equal to the observed value
  alpha <- rowSums(sweep(x, 1, x[, 1], '>=')) / ncol(x)
  
  # Calculate the p-values
  p_val <- apply(cbind(alpha, 1 - alpha), 1, min) * 2
  
  return(p_val)
}

#' @title .one_tailed <internal>
#'
#' @description perform a one-tailed test comparing permuted to observed results
#'
#' @param x *matrix* where the first column is observations and the remaining
#' are permutations
#' @param direction *string* comparison formula, defaults to '<='
#'
#' @details when `direction` = '<=', tests whether observations are *greater*
#' than expected. For HFA analysis, this would test for specialization.
#'
#' @returns a *vector* of the probability observations are different from
#' permutations in the `direction` indicated.
#'
#' @seealso `specialist_test()`, `generalist_test()`

.one_tailed <- function(x, direction = c('<=', '>')) {
  direction <- match.arg(direction)
  
  # Calculate the number of permutations meeting the condition
  alpha <- rowSums(sweep(x, 1, x[, 1], direction))
  
  # Calculate the proportion of permutations that meet the condition
  alpha <- alpha / ncol(x)
  
  # Calculate p-values
  p_val <- 1 - alpha
  
  return(p_val)
}


#' @title specialist_test
#'
#' @description test whether a line is likely a specialist.
#'
#' @param permute_results *list*: the result from a call to `permute_hfa()`
#'
#' @details Runs a one-tailed test of whether observed home field advantage
#' is *greater* than expected based on chance (permutations). Run `permute_hfa()`
#' first and feed the output into this function.
#'
#' @returns a *list*, modifying `permute_results` to include the test
#' in `permute_results$home_field['p_specialist']`. This is the *p*-value,
#' the probability of being called a specialist due to chance.
#'
#' @seealso `permute_hfa()`, `generalist_test()`, `.one_tailed()`
#'
#' @export

specialist_test <- function(permute_results) {
  perms <- permute_results$perms
  
  # Compare permutations to observations (observations on right of inequality)
  test <- do.call(c, lapply(perms, .one_tailed, direction = '<='))
  
  permute_results$home_field$p_specialist <- test
  
  return(permute_results)
}


#' @title generalist_test
#'
#' @description test whether a line is likely a generalist.
#'
#' @param permute_results *list*: the result from a call to `permute_hfa()`
#'
#' @details Runs a one-tailed test of whether observed home field advantage
#' is *less* than expected based on chance (permutations). Run `permute_hfa()`
#' first and feed the output into this function.
#'
#' @returns a *list*, modifying `permute_results` to include the test
#' in `permute_results$home_field['p_generalist']`. This is the *p*-value,
#' the probability of being called a generalist due to chance.
#'
#' @seealso `permute_hfa()`, `specialist_test()`, `.one_tailed()`
#'
#' @export

generalist_test <- function(permute_results) {
  perms <- permute_results$perms
  
  # Compare permutations to observations (observations on right of inequality)
  test <- do.call(c, lapply(perms, .one_tailed, direction = '>'))
  
  permute_results$home_field$p_generalist <- test
  
  return(permute_results)
}

