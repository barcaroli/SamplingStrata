# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

buildStrataDFSpatial <- function(dataset, fitting = as.numeric( c(1.0)), range = as.numeric( c(0.0)), kappa = 3.0, progress = FALSE, verbose = FALSE) {
    .Call(`_SamplingStrata_buildStrataDFSpatial`, dataset, fitting, range, kappa, progress, verbose)
}

