unlink("C:/Users/Giulio/AppData/Local/R/win-library/4.5/00LOCK-SamplingStrata", recursive = TRUE)
Sys.setenv(PATH = paste("C:/rtools44/usr/bin", "C:/rtools44/mingw64/bin", Sys.getenv("PATH"), sep = ";"))
devtools::install_github("barcaroli/SamplingStrata", build_vignettes = FALSE)
