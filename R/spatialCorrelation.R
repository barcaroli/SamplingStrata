# ==============================================================================
# File:      spatialCorrelation.R
# Functions: vgmModelNames(), spatialCorr() [internal]
#
# Purpose:
#   Variogram models of the prediction errors used by the "spatial" method.
#   The covariance of the prediction errors of two distinct units i, j is
#     Cov(e_i, e_j) = (1 - nugget) * s_i * s_j * rho(d_ij)
#   where rho is the correlation function of the variogram model fitted on the
#   residuals (r = range):
#     exponential  rho(d) = exp(-kappa * d / r)
#     spherical    rho(d) = 1 - 1.5 (d / r) + 0.5 (d / r)^3 for d < r, 0 for d >= r
#     gaussian     rho(d) = exp(-kappa * (d / r)^2)
#   kappa is not used by the spherical model. The same functions are coded in
#   C++ in buildStrataDFSpatial.cpp (used by the genetic algorithm); the R
#   version is used by aggrStrataSpatial() (K-means initial solutions).
# ==============================================================================

# Canonical names of the variogram models: full names or gstat abbreviations
# ("Exp", "Sph", "Gau"), case-insensitive. A dataframe (one row per domain)
# is converted column by column, keeping its row names.
vgmModelNames <- function(x) {
  if (is.data.frame(x)) {
    x[] <- lapply(x, vgmModelNames)
    return(x)
  }
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x) || length(x) == 0 || anyNA(x))
    stop("'vgm_model' must be \"exponential\" (\"Exp\"), \"spherical\" (\"Sph\") ",
         "or \"gaussian\" (\"Gau\")", call. = FALSE)
  codes <- c(exponential = "exponential", exp = "exponential",
             spherical = "spherical", sph = "spherical",
             gaussian = "gaussian", gau = "gaussian")
  out <- unname(codes[tolower(x)])
  if (anyNA(out))
    stop("'vgm_model' must be \"exponential\" (\"Exp\"), \"spherical\" (\"Sph\") ",
         "or \"gaussian\" (\"Gau\"); found: ",
         paste(unique(x[is.na(out)]), collapse = ", "), call. = FALSE)
  out
}

# Correlation of the prediction errors at distance(s) d (range > 0)
spatialCorr <- function(d, range, kappa, vgm_model = "exponential") {
  switch(vgm_model,
         exponential = exp(-kappa * d / range),
         spherical = {
           h <- d / range
           ifelse(h < 1, 1 - 1.5 * h + 0.5 * h^3, 0)
         },
         gaussian = exp(-kappa * (d / range)^2),
         stop("unknown variogram model: ", vgm_model, call. = FALSE))
}
