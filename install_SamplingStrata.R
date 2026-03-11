Sys.setenv(PATH = paste("C:/rtools44/usr/bin", "C:/rtools44/mingw64/bin", Sys.getenv("PATH"), sep = ";"))
devtools::install_github("barcaroli/SamplingStrata", build_vignettes = TRUE)
