#' Illinois Corn Yield Data
#' 
#' A cleaned-up version of data from the University of Illinois Extension
#' corn yield trials from 1999-2018. Companies entered elite hybrids into this
#' trial, often targetted toward subregions of Illinois. The study designers
#' chose a few of them each year to plant in all locations as "checks". Only
#' these checks are included. 
#' 
#' The data is highly unbalanced and sparse: Hybrids and sites both turn over
#' frequently and check hybrids in one year may be experimental in another. If a
#' hybrid was a check in any year, it is included.
#' 
#' Data was also processed gently for consistency of site and company names.
#' 
#' @format `corn`: 
#' A data frame with 2,051 rows and 5 columns:
#' \describe{
#'   \item{YEAR}{Year of observation}
#'   \item{COMPANY}{Company that entered the hybrid}
#'   \item{HYBRID}{Hybrid name. Approximately analogous to genotype}
#'   \item{SITE}{Town in Illinois where the trial was conducted}
#'   \item{YIELD}{Grain yield in standard US bushels per acre, standardized to
#'   15.5\% moisture. One standard bushel of corn grain weighs 56 pounds.}
#' }
#' 
#' @source <https://vt.cropsci.illinois.edu/corn/>
"corn"