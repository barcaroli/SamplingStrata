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
owd <- setwd("vignettes")
knitr::knit("SamplingStrata.Rmd.orig", output = "SamplingStrata.Rmd")
# knitr::knit("drake_continuous.Rmd.orig", output = "drake_continuous.Rmd")
knitr::knit("models.Rmd.orig", output = "models.Rmd")
knitr::knit("spatial.Rmd.orig", output = "spatial.Rmd")
setwd(owd)
