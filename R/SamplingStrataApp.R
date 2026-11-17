#' Launch the SamplingStrata Shiny application
#'
#' Starts an interactive application that guides the user through the whole
#' SamplingStrata workflow: data loading, variable selection and
#' categorisation, frame and strata construction, precision constraints,
#' optional take-all strata, optimisation, evaluation, sample selection, and
#' export of a reproducible R script.
#'
#' The application requires the packages \pkg{shiny}, \pkg{shinydashboard},
#' \pkg{shinyjs}, \pkg{DT}, \pkg{ggplot2} and \pkg{corrplot}.
#'
#' @param launch.browser Logical; if \code{TRUE} (default) the app is opened in
#'   the default web browser.
#' @param ... Further arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return Called for its side effect (starts the app); returns invisibly
#'   when the app is closed.
#'
#' @examples
#' if (interactive()) {
#'   SamplingStrataApp()
#' }
#'
#' @export
SamplingStrataApp <- function(launch.browser = TRUE, ...) {
  needed <- c("shiny", "shinydashboard", "shinyjs", "DT", "ggplot2", "corrplot")
  missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1),
                                 quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    if (interactive() &&
        isTRUE(utils::askYesNo(paste0(
          "The following packages are needed by the app but not installed: ",
          paste(missing_pkgs, collapse = ", "), ". Install them from CRAN?")))) {
      utils::install.packages(missing_pkgs)
      missing_pkgs <- missing_pkgs[!vapply(missing_pkgs, requireNamespace,
                                           logical(1), quietly = TRUE)]
    }
    if (length(missing_pkgs) > 0)
      stop("Please install the required package(s): ",
           paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  appDir <- system.file("shiny", "SamplingStrataApp", package = "SamplingStrata")
  if (!nzchar(appDir))
    stop("Could not find the app directory; try reinstalling SamplingStrata.",
         call. = FALSE)
  shiny::runApp(appDir, launch.browser = launch.browser, ...)
}
