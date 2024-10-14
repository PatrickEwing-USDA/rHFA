#' @title .calculate_temporal_hfa <internal>
#' 
#' @description calculate hfa for each year
#' 
#' @param x data.frame
#' @param ff formula
#' @param year column in data.frame `x` representing year.
#' 
#' @details the formula should contain a term, <year>_num:is_home, and also specify a response (e.g. pheno)
#' 
#' @return a data.frame with the coefficients plus relevant 
#' 
#' @importFrom stats p.adjust
#' @importFrom utils type.convert

.calculate_temporal_hfa <- function(x, ff, year) {
  y <- as.character(ff)[2]
  
  # get the coefficients
  home_coef <- summary(lm(ff, x))$coef
  
  # format
  home_coef <- home_coef[grepl(':is_homeTRUE', rownames(home_coef)), ]
  
  rownames(home_coef) <- gsub(':is_homeTRUE', '', rownames(home_coef))
  rownames(home_coef) <- gsub(year, '', rownames(home_coef))
  
  # return values
  home_coef <- data.frame(year = rownames(home_coef),
                          year_num = type.convert(rownames(home_coef)),
                          home_coef,
                          p.adj = p.adjust(home_coef[, 'Pr(>|t|)']),
                          stringsAsFactors = FALSE)
  colnames(home_coef) <- gsub('year', year, colnames(home_coef))
  colnames(home_coef) <- gsub('Pr...t..', 'p.value', colnames(home_coef))
  colnames(home_coef) <- gsub('Std..Error', 'Std.Error', colnames(home_coef))
  
  return(home_coef)
}

#' @title temporal_hfa
#' 
#' @description Calculate trend in home field advantage across time. The formula is:
#' 
#' `pheno ~ year + site + geno + year:site + year:is_home`. 
#' 
#' The year:is_home coefficient is the home field advantage in each year. 
#' 
#' @param data data.frame
#' @param site column in `data` indicating spatial environment
#' @param year column in `data` indicating temporal environment
#' @param geno column in `data` indicating genotype
#' @param pheno column in `data` indicating phenotype
#' @param popn column in `data` indicating groups of genotypes
#' @param blup_home identify home site using BLUP?
#' @param formula_override formula override, default is NA.
#' @param parallel use parallel processing? Default is TRUE.
#' 
#' @return A list with ANOVA results, model, temporal HFA coefficients, and the used formulas.
#' 
#' @importFrom parallel detectCores mclapply
#' @importFrom car Anova
#' @importFrom stats formula lm
#' 
#' @export

temporal_hfa <- function(data, site, year, geno, pheno, popn = NA, 
                         blup_home = TRUE, 
                         formula_override = NA,
                         parallel = TRUE) {
  
  if (grepl('mingw', version$os)) parallel <- FALSE
  ncpu <- ifelse(parallel, detectCores(), 1)
  
  # set up dataframe
  if (is.na(popn)) {
    dd <- list(data)
  } else {
    dd <- split(data, data[, popn])
  }
  
  # formula
  if (is.na(formula_override)) {
    ff <- gsub('pheno', pheno, "pheno ~ year + site + geno + year:site + year:is_home")
    ff <- gsub('site', site, ff)
    ff <- gsub('year', year, ff)
    ff <- gsub('geno', geno, ff)
    ff <- formula(ff)
  } else {
    ff <- formula_override
  }
  
  # has home already been identified?
  id_home_site <- !any(grepl('is_home', colnames(data)))
  
  if (id_home_site) {
    dd <- mclapply(dd, function(x) {
      id_home(x, site, year, geno, pheno, blup = blup_home, verbose = FALSE)
    }, mc.cores = ncpu)
  }
  
  # calculate annual hfa
  out <- do.call(rbind, mclapply(dd, .calculate_temporal_hfa, ff, year, mc.cores = ncpu))
  
  # population, if specified
  if (!is.na(popn)) {
    popn_id <- sapply(strsplit(rownames(out), '\\.'), '[', 1)
    
    out <- data.frame(
      popn = popn_id,
      out,
      stringsAsFactors = FALSE
    )
    colnames(out) <- gsub('popn', popn, colnames(out))
  }
  
  # year/temporal
  if (is.factor(data[, year])) {
    out[, year] <- factor(out[, year], levels = levels(data[, year]))
  }
  
  # run anova
  year_num <- paste0(year, '_num')
  ff_aov <- paste('Estimate', year_num, sep = '~')
  if (!is.na(popn)) {
    ff_aov <- paste(ff_aov, '*', popn)
  }
  
  out_lm <- lm(formula(ff_aov), data = out)
  out_aov <- Anova(out_lm)
  
  out <- list(anova = out_aov,
              model = out_lm,
              temporal_hfa = out,
              formulas = c(hfa = ff, anova = ff_aov))
  
  return(out)
}
