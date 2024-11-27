#' @title print.rHFA
#' 
#' @description Print method for `class("rHFA")`
#' 
#' @param x object of `class("rHFA")`. Output from `permute_hfa()`.
#' 
#' @export

print.rHFA <- function(x, ...) {
  ff <- x$formula |>
    as.character()
  ff <- paste(ff[2], ff[1], ff[3])
  # ff <- paste(ff[2], ff[3], sep=ff[1])
  out <- x$home_field
  is_numeric <- sapply(out, is.numeric)
  out[is_numeric] <- lapply(out[is_numeric], signif, 3)
  
  nperms <- ncol(x$perms)-1
  nperms <- paste('Results of', nperms, 'permutations:')
  
  cat('Formula:\n', ff)
  cat('\n\n\n')
  cat(nperms, '\n\n')
  print(out)
  invisible(x)
}
