#' @title get_ss
#'
#' @description Calculate the sum of squares for a model and present a summary of
#' the predictors, their sum of squares, percentage variance explained, F-values,
#' and p-values in a data frame. It can work with both an ANOVA model or fit a model
#' using `car::Anova`.
#'
#' @param model A model object or an ANOVA object.
#' @param test.statistic A character string indicating the test statistic to use. Defaults to 'F'.
#' @param ... Additional arguments passed to `Anova()`.
#'
#' @return A data frame with the statistical measures, including sum of squares, percentage variance,
#' F-values, and p-values.
#' @importFrom car Anova
#' @export

get_ss <- function(model, test.statistic = 'F', ...) {
  
  # Check if the input is an anova object
  if ("anova" %in% class(model)) {
    a <- model
  } else {
    # Run car::Anova() if model is not an anova object
    message("Running car::Anova(model)")
    
    if (!requireNamespace("car", quietly = TRUE)) {
      stop("Package 'car' is required but not installed. Please install 'car' to use this function.")
    }
    
    a <- car::Anova(model, ...)
  }
  
  # Calculate percentage variance
  pVar <- round(a[, 1] / sum(a[, 1]) * 100, 2)
  
  # Prepare output data frame
  out <- data.frame(
    SumSq = a[, 1],
    PercentVar = pVar,
    df = a[, 2],
    F_value = round(a[, 3], 4),
    p_value = signif(a[, 4], 3)
  )
  
  # Set row and column names
  rownames(out) <- rownames(a)
  colnames(out) <- c("Sum Sq", "Percent Var", "df", "F value", "Pr(>F)")
  
  # Set 'is_home' p-value to NA if applicable
  if ("is_home" %in% rownames(out)) {
    out["is_home", "p_value"] <- NA
  }
  
  return(out)
}
