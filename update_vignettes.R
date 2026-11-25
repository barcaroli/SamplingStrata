# ------------------------------------------------------------------------------
# Pre-compiled vignettes
# The *.Rmd.orig files are the SOURCES of the vignettes: their R chunks are
# executed here. knitr::knit() turns each of them into a *.Rmd file in which
# code, printed results and figures are already written as static content, so
# that R CMD build / R CMD check (and CRAN) do not rerun the optimisations.
# Edit only the .orig files, then (after reinstalling the package) run this
# script from the package root folder.
# The working directory is set to 'vignettes' so that the figures are written
# in vignettes/fig/ (fig.path = "fig/" in each .orig), where the generated .Rmd
# files look for them. The folder is not called 'figure' because R CMD check
# reports a directory with that name as a leftover from knitr.
# ------------------------------------------------------------------------------
# Each vignette is knitted in a new, clean R process (xfun::Rscript_call):
# the result does not depend on the packages and objects of the current
# session, and the version of SamplingStrata installed on disk is used even if
# an older one is still loaded here.
knit_clean <- function(input, output)
  xfun::Rscript_call(knitr::knit, list(input = input, output = output))
owd <- setwd("vignettes")
knit_clean("SamplingStrata.Rmd.orig", output = "SamplingStrata.Rmd")
# knitr::knit("drake_continuous.Rmd.orig", output = "drake_continuous.Rmd")
knit_clean("models.Rmd.orig", output = "models.Rmd")
knit_clean("spatial.Rmd.orig", output = "spatial.Rmd")
setwd(owd)
