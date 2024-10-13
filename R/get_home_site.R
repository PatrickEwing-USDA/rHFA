#' @title get_home_site
#'
#' @description Extract home location information for each genotype based on the results of `id_home()`.
#'
#' @param data A data.frame from a call to `id_home()`.
#' @param geno A column in `data` with genotype.
#' @param site A column in `data` with location environment.
#'
#' @return A `data.frame` with columns `geno`, `site`, and `home_years`, which indicates the number
#' of years a genotype appeared at a home site.
#' @export
#' 
get_home_site <- function(data, geno, site) {
  
  # Check if 'is_home' column is present in data
  if (!("is_home" %in% colnames(data))) {
    stop("Did you call `id_home()` on `data`?")
  }
  
  # Filter rows where 'is_home' is TRUE
  is_home <- data[, 'is_home']
  out <- data[is_home, c(geno, site)]
  
  # Count occurrences of each genotype
  counts <- tapply(out[, geno], out[, geno], length)
  
  # Add 'home_years' column indicating number of home years
  out <- unique(out)
  out[, 'home_years'] <- counts[out[, geno]]
  
  return(out)
}

#' @title home_distance_spatial
#'
#' @description Calculate the distance between home sites for genotypes based on spatial locations.
#'
#' @param data A data.frame from a call to `id_home()`.
#' @param locations A `SpatialPointsDataFrame` (recommended) or a data.frame with columns `site`, `lat`, and `long`.
#' @param geno Genotype column in `data` and in `locations`.
#' @param site Spatial environment column in `data`.
#' @param lat Latitude or y-coordinate column in `locations`.
#' @param long Longitude or x-coordinate column in `locations`.
#' @param great_circle A logical indicating whether coordinates are on an ellipse (set TRUE for lat/long).
#'
#' @return A distance matrix. Row and column names represent genotypes, and values represent distances
#' between the home sites of genotype pairs.
#' @seealso `spDists` from the `"sp"` package.
#' 
#' @importFrom sp spDists
#' @export

home_distance_spatial <- function(data, locations, geno, site, lat = NA, long = NA, great_circle = TRUE) {
  
  if (!("is_home" %in% colnames(data))) {
    stop("Did you call `id_home()` on `data`?")
  }
  
  # Check that site exists in locations
  if (!(site %in% colnames(locations))) {
    stop(paste(site, "is not in the columns of `locations`."))
  }
  
  # Validate that all site levels in data exist in locations
  unmatched_sites <- data[, site] %in% locations[, site] == FALSE
  if (any(unmatched_sites)) stop("Some sites in `data` are not found in `locations`.")
  
  # Calculate distances using spatial data
  if (inherits(locations, "SpatialPoints")) {
    distances <- sp::spDists(locations)
  } else {
    distances <- sp::spDists(as.matrix(locations[, c(long, lat)]), longlat = great_circle)
  }
  
  colnames(distances) <- locations[[site]]
  rownames(distances) <- locations[[site]]
  
  # Get home site data
  home_data <- get_home_site(data, geno, site)
  
  # Create the genotype distance matrix
  genotypes <- as.character(home_data[, geno])
  home_sites <- as.character(home_data[, site])
  names(home_sites) <- genotypes
  n <- length(genotypes)
  
  dist_matrix <- matrix(nrow = n, ncol = n, dimnames = list(genotypes, genotypes))
  
  for (i in genotypes) {
    for (j in genotypes) {
      dist_matrix[i, j] <- distances[home_sites[i], home_sites[j]]
    }
  }
  
  return(dist_matrix)
}

#' @title home_distance_environmental
#'
#' @description Calculate the environmental distance between home sites of genotypes.
#'
#' @param data A data.frame from a call to `id_home()`.
#' @param locations A data.frame with a column `site`.
#' @param geno Genotype column in `data`.
#' @param site Site environment column in `data` and `locations`.
#' @param vars A vector of columns in `locations` with environmental measurements.
#' @param scale A logical indicating whether `locations[, vars]` should be centered and scaled. Defaults to TRUE.
#' @param ... Additional arguments passed to `dist()`.
#'
#' @return A distance matrix where rows and columns represent genotypes, and values represent distances
#' between the home sites of the genotype pairs.
#' 
#' @importFrom stats dist
#' @export

home_distance_environmental <- function(data, locations, geno, site, vars = c(), scale = TRUE, ...) {
  
  vars <- unlist(vars)  # Ensure vars is a character vector
  
  # Check for 'is_home' column in data
  if (!("is_home" %in% colnames(data))) {
    stop("Did you call `id_home()` on `data`?")
  }
  
  # Validate the presence of the site column in locations
  if (!(site %in% colnames(locations))) {
    stop(paste(site, "is not in the columns of `locations`."))
  }
  
  # Validate that all site levels in data exist in locations
  unmatched_sites <- data[, site] %in% locations[, site] == FALSE
  if (any(unmatched_sites)) stop("Some sites in `data` are not found in `locations`.")
  
  # Prepare environmental data for distance calculation
  env_data <- locations[, vars, drop = FALSE]
  rownames(env_data) <- locations[, site]
  
  if (scale) {
    env_data <- scale(env_data)  # Center and scale environmental variables
  }
  
  if (any(is.na(env_data))) {
    warning("NAs found in `locations[, vars]`.")
  }
  
  # Calculate environmental distances
  env_distances <- as.matrix(dist(env_data, ...))
  
  # Get home site data
  home_data <- get_home_site(data, geno, site)
  
  # Create the genotype distance matrix
  genotypes <- as.character(home_data[, geno])
  home_sites <- as.character(home_data[, site])
  names(home_sites) <- genotypes
  n <- length(genotypes)
  
  dist_matrix <- matrix(nrow = n, ncol = n, dimnames = list(genotypes, genotypes))
  
  for (i in genotypes) {
    for (j in genotypes) {
      dist_matrix[i, j] <- env_distances[home_sites[i], home_sites[j]]
    }
  }
  
  return(dist_matrix)
}

