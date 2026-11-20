## Adds the argument mut_factor (default 0.05) to optimStrata() and to the
## functions that implement the "continuous" and "spatial" methods.
## Run from the root folder of the SamplingStrata source package:
##   setwd("D:/Google Drive/Sampling/SamplingStrata 1.6")
##   source("apply_mut_factor.R")
## then regenerate the documentation with devtools::document().

edit_file <- function(path, fixed = list(), regex = list()) {
  raw <- readChar(path, file.info(path)$size, useBytes = TRUE)
  crlf <- grepl("\r\n", raw, fixed = TRUE)
  s <- gsub("\r\n", "\n", raw, fixed = TRUE)
  for (r in fixed) {
    n <- lengths(regmatches(s, gregexpr(r[1], s, fixed = TRUE)))
    if (n != as.integer(r[3])) stop(path, ": expected ", r[3], " match(es), found ", n,
                                    " for: ", substr(r[1], 1, 60))
    s <- gsub(r[1], r[2], s, fixed = TRUE)
  }
  for (r in regex) {
    n <- lengths(regmatches(s, gregexpr(r[1], s, perl = TRUE)))
    if (n != as.integer(r[3])) stop(path, ": expected ", r[3], " match(es), found ", n,
                                    " for regex: ", r[1])
    s <- gsub(r[1], r[2], s, perl = TRUE)
  }
  if (crlf) s <- gsub("\n", "\r\n", s, fixed = TRUE)
  writeChar(s, path, eos = NULL, useBytes = TRUE)
  cat("updated", path, "\n")
}

doc <- paste0(
"#' @param mut_factor Mutation factor (used only by the \"continuous\" and \"spatial\" methods): a mutated\n",
"#'   gene (a cut point, encoded in [0, 1] on the range of the stratification variable) is multiplied by\n",
"#'   (1 + mut_factor) or (1 - mut_factor), with equal probability. Small values allow a fine local\n",
"#'   search around the current boundaries. Default is 0.05 (in previous versions it was fixed to 0.5).\n")
add_doc <- function(path) {
  raw <- readChar(path, file.info(path)$size, useBytes = TRUE)
  crlf <- grepl("\r\n", raw, fixed = TRUE)
  s <- gsub("\r\n", "\n", raw, fixed = TRUE)
  if (grepl("@param mut_factor", s, fixed = TRUE)) stop(path, ": mut_factor already documented")
  i <- regexpr("#' @param mut_chance", s, fixed = TRUE)
  rest <- substring(s, i + 5)
  j <- i + 5 + regexpr("#' @param", rest, fixed = TRUE) - 2
  s <- paste0(substr(s, 1, j), doc, substring(s, j + 1))
  if (crlf) s <- gsub("\n", "\r\n", s, fixed = TRUE)
  writeChar(s, path, eos = NULL, useBytes = TRUE)
}

edit_file("R/optimStrata.R", fixed = list(
  c("mut_chance = NA, elitism_rate = 0.2, suggestions = NULL,",
    "mut_chance = NA, mut_factor = 0.05, elitism_rate = 0.2, suggestions = NULL,", 1),
  c("mut_chance = mut_chance, elitism_rate = elitism_rate, \n                                highvalue = 1e+08, suggestions = suggestions, realAllocation = TRUE, \n                                writeFiles = writeFiles, showPlot = showPlot, parallel = parallel, \n                                cores = cores)",
    "mut_chance = mut_chance, mut_factor = mut_factor,\n                                elitism_rate = elitism_rate, \n                                highvalue = 1e+08, suggestions = suggestions, realAllocation = TRUE, \n                                writeFiles = writeFiles, showPlot = showPlot, parallel = parallel, \n                                cores = cores)", 1),
  c("iter = iter, pops = pops, mut_chance = mut_chance, \n",
    "iter = iter, pops = pops, mut_chance = mut_chance, mut_factor = mut_factor,\n", 1)))
add_doc("R/optimStrata.R")

edit_file("R/optimizeStrata2.R",
  fixed = list(c("minnumstr = 2, iter = 50, pops = 20, mut_chance = NA, elitism_rate = 0.2,",
                 "minnumstr = 2, iter = 50, pops = 20, mut_chance = NA, mut_factor = 0.05, elitism_rate = 0.2,", 1)),
  regex = list(c("(?s)(strataGenalg2\\(errors = erro\\[\\[i\\]\\],.*?showPlot)\\)",
                 "\\1, mut_factor = mut_factor)", 3)))
add_doc("R/optimizeStrata2.R")

edit_file("R/optimizeStrataSpatial.R",
  fixed = list(c("                                   mut_chance = NA, \n",
                 "                                   mut_chance = NA, \n                                   mut_factor = 0.05, \n", 1)),
  regex = list(c("(?s)(strataGenalgSpatial\\(errors = erro\\[\\[i\\]\\],.*?showPlot)\\)",
                 "\\1, mut_factor = mut_factor)", 3)))
add_doc("R/optimizeStrataSpatial.R")

for (f in c("R/strataGenalg2.R", "R/strataGenalgSpatial.R"))
  edit_file(f, fixed = list(
    c("    writeFiles, \n    showPlot\n)", "    writeFiles, \n    showPlot,\n    mut_factor = 0.05\n)", 1),
    c("mutationFactor=0.5,", "mutationFactor=mut_factor,", 1)))

edit_file("NEWS.Rmd", fixed = list(c(
  "* The function *prepareSuggestion* that prepares suggestion for the optimization with \"continuous\" method has been completely redefined and improved.\n",
  "* The function *prepareSuggestion* that prepares suggestion for the optimization with \"continuous\" method has been completely redefined and improved.\n\n* New argument *mut_factor* of *optimStrata* (methods \"continuous\" and \"spatial\"): amplitude of the mutation of the cut points in the genetic algorithm, previously fixed to 0.5. The new default, 0.05, allows a finer local search of the boundaries.\n", 1)))

cat("\nDone. Now run devtools::document() to update the files in man/ (and rebuild NEWS.md from NEWS.Rmd).\n")
