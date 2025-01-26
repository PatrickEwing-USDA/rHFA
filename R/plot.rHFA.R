#' @title plot.rHFA
#' 
#' @description Default plotting function to visualize permutation results of HFA analysis
#' 
#' @param x object of class `rHFA`, i.e., the output of `permute_hfa`
#' @param color_vector factor vector of length `nrow(x$data)` specifying grouping for colors, if desired.
#' 
#' Infers the level on which HFA was calculated and plots observed HFA over the range of permuted
#' HFA. Permuted values are shown by caterpillars in gray showing median, 90%, and 95% confidence intervals.
#' 
#' Colored points are observed.
#' 
#' The output is a ggplot grob, so you can modify it accordingly. 
#' 
#' @returns
#' A ggplot object.
#' 
#' @seealso `permute_hfa()`, `ggplot()`
#' 
#'@importFrom ggplot2 ggplot stat_summary geom_point aes labs scale_color_brewer scale_color_manual theme_minimal theme coord_flip element_blank
#'@importFrom stats median quantile
#'
#'@export


.reshaper <- function(y, x) {
  dd = x[[y]]
  nc <- ncol(dd)
  nr <- nrow(dd)
  
  id <- rownames(dd) |>
    gsub(':is_home', '', x=_) |>
    rep(nc)
  permutation <- lapply(colnames(dd), rep, nr) |>
    unlist()
  value <- matrix(dd, ncol=1)
  dd <- data.frame(
    group=y,
    id=id,
    value=value,
    permutation=permutation
  )
}

#' @title plot_data
#' 
#' @description extract data used for plot.rHFA
#' 
#' @param x object of class `rHFA`, i.e., the output of `permute_hfa`
#' @param expected logical. Add expected `pheno` to the obs? Defaults to FALSE
#' 
#' Reformat `x` into easy-to-use data.frames for plotting or modifying output
#' of rHFA plotting functions.
#' 
#' Expected pheno is the expected yield, etc. of each level (ex genotype). 
#' 
#' @returns list of two `data.frame`s:
#' - "obs" are the observations, from `x$home_field`.
#' - "perms" are the permutations in tall format, from `x$perms`
#' 
#' @seealso `permute_hfa`
#' 
#' @export


plot_data <- function(x, expected=FALSE) {
  stopifnot('rHFA' %in% class(x))
  obs <- x$home_field[1:3]
  perms <- x$perms
  
  if (!is.list(perms)) perms = list(popn1=perms)
  
  perms <- lapply(names(perms), .reshaper, perms) |>
    do.call(rbind, args=_)
  colnames(perms)[1:3] = colnames(obs)
  
  out <- list(obs=obs,
              perms=perms)
  return(out)
}


#' @title plot.rHFA
#' 
#' @description Default plotting function to visualize permutation results of HFA analysis
#' 
#' @param x object of class `rHFA`, i.e., the output of `permute_hfa`
#' @param color_vector factor vector of length `nrow(x$data)` specifying grouping for colors, if desired.
#' 
#' Infers the level on which HFA was calculated and plots observed HFA over the range of permuted
#' HFA. Permuted values are shown by caterpillars in gray showing median, 90%, and 95% confidence intervals.
#' 
#' Colored points are observed.
#' 
#' The output is a ggplot grob, so you can modify it accordingly. 
#' 
#' @returns
#' A ggplot object.
#' 
#' @seealso `permute_hfa()`, `ggplot()`, `plot_data()`
#' 
#'@importFrom ggplot2 ggplot stat_summary geom_point aes labs scale_color_brewer scale_color_manual theme_minimal theme coord_flip element_blank
#'@importFrom stats median quantile
#'
#'@export
plot.rHFA <- function(x, grouping=NULL) {
  perms <- plot_data(x)$perms
  
  CI90 = function(x) {
    data.frame(y=median(x),
               ymin=quantile(x, .05),
               ymax=quantile(x, .95))
  }
  CI95 = function(x) {
    data.frame(y=median(x),
               ymin=quantile(x, 0.025),
               ymax=quantile(x, 0.975))
  }
  
  perms$grouping = perms$permutation
  if (!is.null(grouping)) {
    is_observed = perm$permutation == 'observed'
    perm[is_observed, 'grouping'] = grouping
  }
  
  xname = names(perms)[2]
  yname = names(perms)[3]
  
  plt = ggplot(perms,
               aes(x=.data[[xname]],
                   y=.data[[yname]])) +
    stat_summary(geom='pointrange',
                 fun.data=CI90,
                 color='darkgray',
                 linewidth=1.5) +
    stat_summary(geom='pointrange',
                 fun.data=CI95,
                 color='darkgray',
                 linewidth=.5) +
    geom_point(data=subset(perms, permutation=='observed'),
               mapping=aes(color=grouping),
               size=3) +
    labs(color=NULL) +
    theme_minimal()
  
  if (x$level %in% c('genotype', 'site')){
    plt = plt +
      coord_flip()
  }
  if (x$level == 'population') {
    plt <- plt + 
      labs(x='Population')
  }
  
  return(plt)
}