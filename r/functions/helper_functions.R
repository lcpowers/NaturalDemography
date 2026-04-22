## Various project functions

ellipse_area <- function(a1,a2,pm){
  
  area = pi*(a1/2)*(a2/2) * (100-pm)/100 /100 # mm^2 to cm^2
  return(area)
  
}

w_mean <- function(v, w) {
  if (all(is.na(v))) return(NA_real_)
  # use weights only where v is not NA
  ww <- w[!is.na(v)]
  vv <- v[!is.na(v)]
  sum(vv * ww) / sum(ww)
}

# Weighted median helper (handles coverage fractions)
w_median <- function(v, w) {
  if (all(is.na(v))) return(NA_real_)
  o <- order(v)
  v <- v[o]; w <- w[o]
  w <- w / sum(w)
  cumw <- cumsum(w)
  v[which.min(abs(cumw - 0.5))]
}

# Coverage-weighted circular mean for aspect (degrees 0–360)
w_circ_mean_deg <- function(values, coverage_fraction) {
  v <- values
  w <- coverage_fraction
  
  # Remove missing pairs of values and weights
  keep <- !is.na(v) & !is.na(w)
  if (!any(keep)) return(NA_real_)
  
  # Subset to valid entries only
  v <- v[keep]
  w <- w[keep]
  
  # Convert degrees to radians
  rad <- v * pi / 180
  
  # Compute weighted vector components
  x <- sum(cos(rad) * w)
  y <- sum(sin(rad) * w)
  
  # Convert back to degrees and wrap into [0, 360)
  mean_angle <- atan2(y, x) %% (2 * pi)
  mean_angle * 180 / pi
}

# Helper function for coverage-weighted range (max - min)
w_range <- function(values, coverage_fraction) {
  v <- values[!is.na(values) & !is.na(coverage_fraction)]
  if(length(v) == 0) return(NA_real_)
  max(v) - min(v)
}

# Potential insolation via hillshade integration over growing-season sun positions
# lat_dd: site latitude in decimal degrees
# Returns SpatRaster [0,1]; higher = more solar exposure
p_insol <- function(dem_spat, lat_dd) {
  trn <- terra::terrain(dem_spat, v = c("slope","aspect"), unit = "radians", neighbors = 8)

  # Solar noon elevations for summer solstice and equinox
  elevs <- c(90 - lat_dd + 23.45, 90 - lat_dd)
  elev_wts <- c(0.6, 0.4)  # summer weighted higher (longer days)
  azims <- seq(0, 315, by = 45)

  hs_accum <- NULL
  for(ei in seq_along(elevs)) {
    hs_az <- NULL
    for(az in azims) {
      hs <- terra::shade(trn[["slope"]], trn[["aspect"]], angle = elevs[ei], direction = az)
      hs_az <- if(is.null(hs_az)) hs else (hs_az + hs)
    }
    hs_az <- hs_az / length(azims)
    hs_accum <- if(is.null(hs_accum)) hs_az * elev_wts[ei] else hs_accum + hs_az * elev_wts[ei]
  }
  hs_accum
}

# Convert odds ratios to predicted probabilities
predict_p <- function(or, p0) {
  odds0 <- p0 / (1 - p0)
  odds1 <- odds0 * or
  p1 <- odds1 / (1 + odds1)
  return(p1)
}  
