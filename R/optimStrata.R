optimStrata <- function (method = c("atomic", "continuous", "spatial"), framesamp, 
                         framecens = NULL, model = NULL, nStrata = NA, errors, alldomains = TRUE, 
                         dom = NULL, strcens = FALSE, minnumstr = 2, iter = 50, pops = 20, 
                         mut_chance = NA, elitism_rate = 0.2, suggestions = NULL, 
                         writeFiles = FALSE, showPlot = TRUE, parallel = TRUE, cores = NA, 
                         fitting = NA, range = NA, kappa = NA) 
{
  check_consecutive_domainvalue <- function(df, df_name) {
    colnames(df) <- toupper(colnames(df))
    domain_ids <- sort(unique(as.integer(as.character(df$DOMAINVALUE))))
    expected_domains <- seq_len(length(domain_ids))
    if (length(domain_ids) == 0 || any(is.na(domain_ids)) || 
        !identical(domain_ids, expected_domains)) {
      stop(paste0("'", df_name, "$DOMAINVALUE' must contain consecutive integer values from 1 to N with no gaps. Found: ", 
                  paste(domain_ids, collapse = ", ")))
    }
  }
  check_spatial_frame <- function(framesamp, errors, nStrata, minnumstr, df_name = "framesamp") {
    if (is.null(framesamp)) {
      stop(paste0("The spatial sampling frame (", df_name, ") dataframe is missing"))
    }
    if (is.null(errors)) {
      stop("The 'precision constraints' (errors) dataframe is missing")
    }

    cols <- names(framesamp)
    y_cols <- grep("^Y[0-9]+$", cols, value = TRUE)
    x_cols <- grep("^X[0-9]+$", cols, value = TRUE)
    var_cols <- grep("^var[0-9]+$", cols, value = TRUE)
    lon_col <- intersect(c("lon", "LON"), cols)[1]
    lat_col <- intersect(c("lat", "LAT"), cols)[1]

    if (length(y_cols) == 0) {
      stop(paste0("No Y* variables found in '", df_name, "'. Use buildFrameSpatial() before optimStrata(method = 'spatial')."))
    }
    if (length(x_cols) == 0) {
      stop(paste0("No X* variables found in '", df_name, "'. Use buildFrameSpatial() before optimStrata(method = 'spatial')."))
    }
    if (length(var_cols) != length(y_cols)) {
      stop(paste0("The number of var* columns in '", df_name, "' must match the number of Y* variables."))
    }
    if (is.na(lon_col) || is.na(lat_col)) {
      stop(paste0("Coordinates lon/lat must be present in '", df_name, "'."))
    }
    if (!("domainvalue" %in% cols)) {
      stop(paste0("Column domainvalue must be present in '", df_name, "'."))
    }

    numeric_cols <- c(x_cols, y_cols, var_cols, lon_col, lat_col)
    not_numeric <- numeric_cols[!vapply(framesamp[numeric_cols], is.numeric, logical(1))]
    if (length(not_numeric) > 0) {
      stop(paste0("These spatial frame columns must be numeric: ", paste(not_numeric, collapse = ", ")))
    }

    non_finite <- vapply(framesamp[numeric_cols], function(x) sum(!is.finite(x)), integer(1))
    if (any(non_finite > 0)) {
      bad <- paste(names(non_finite)[non_finite > 0], non_finite[non_finite > 0], sep = "=", collapse = ", ")
      stop(paste0("Non-finite values found in spatial frame columns: ", bad,
                  ". Remove or impute them before optimStrata(method = 'spatial')."))
    }

    non_positive_var <- vapply(framesamp[var_cols], function(x) sum(x <= 0), integer(1))
    if (any(non_positive_var > 0)) {
      bad <- paste(names(non_positive_var)[non_positive_var > 0], non_positive_var[non_positive_var > 0], sep = "=", collapse = ", ")
      stop(paste0("Non-positive variance values found in spatial frame columns: ", bad,
                  ". Filter records with var* <= 0 before optimStrata(method = 'spatial')."))
    }

    domain_size <- table(framesamp$domainvalue)
    small_domains <- domain_size[domain_size < minnumstr]
    if (length(small_domains) > 0) {
      bad <- paste(names(small_domains), as.integer(small_domains), sep = "=", collapse = ", ")
      stop(paste0("Some domains have fewer units than minnumstr=", minnumstr, ": ", bad,
                  ". Accumulate small domains or lower minnumstr."))
    }

    if (!all(is.na(nStrata))) {
      domains <- sort(unique(as.integer(as.character(framesamp$domainvalue))))
      nstrata_check <- nStrata
      if (length(nstrata_check) == 1) {
        nstrata_check <- rep(nstrata_check, length(domains))
      }
      if (length(nstrata_check) != length(domains)) {
        stop("Length of nStrata must be 1 or equal to the number of domains for method = 'spatial'.")
      }
      names(nstrata_check) <- domains
      max_feasible <- floor(as.integer(domain_size[as.character(domains)]) / minnumstr)
      too_many <- nstrata_check > max_feasible
      if (any(too_many, na.rm = TRUE)) {
        bad <- paste(domains[too_many], "requested", nstrata_check[too_many],
                     "max", max_feasible[too_many], sep = "=", collapse = ", ")
        stop(paste0("nStrata is not feasible for some domains with minnumstr=", minnumstr, ": ", bad))
      }
    }

    invisible(TRUE)
  }
  if (!(method %in% c("atomic", "continuous", "spatial"))) 
    stop("Method should be one in 'atomic','continuous','spatial'")
  if (any(nStrata < 2)) 
    stop("There is at least one element in nStrata that is < 2")
  if (alldomains == TRUE & !is.null(dom)) 
    stop("Processing of all domains set TRUE, but a given domain has been indicated")
  if (method == "atomic") {
    if (is.null(errors)) 
      stop("The 'precision constraints' (errors) dataframe is missing")
    if (is.null(framesamp)) 
      stop("The 'sampling frame' (framesamp) dataframe is missing")
    if (!is.null(framesamp)) 
      checkInput(errors, sampframe = framesamp)
    check_consecutive_domainvalue(framesamp, "framesamp")
    check_consecutive_domainvalue(errors, "errors")
    if (!is.null(framecens)) 
      checkInput(errors, sampframe = framecens)
    nvarX <- length(grep("X", colnames(framesamp)))
    for (i in (1:nvarX)) {
      st <- paste0("if (!is.numeric(framesamp$X", i, ")) stop('Stratification variable(s) must be numeric - Transform them to numeric from factor or character before run optimization')")
      eval(parse(text = st))
    }
    strata <- buildStrataDF(framesamp, model = model, progress = FALSE)
    if (!is.null(framecens)) {
      cens <- buildStrataDF(framecens, model = model, 
                            progress = FALSE)
      strcens <- TRUE
    }
    if (is.null(framecens)) {
      cens <- NULL
      strcens <- FALSE
    }
    if (!is.na(nStrata[1])) {
      initialStrata <- nStrata
      addStrataFactor <- 0
    }
    if (is.na(nStrata[1])) {
      initialStrata <- NA
      addStrataFactor <- 0
    }
    solut <- optimizeStrata(errors = errors, strata = strata, 
                            cens = cens, strcens = strcens, alldomains = alldomains, 
                            dom = dom, initialStrata = initialStrata, addStrataFactor = addStrataFactor, 
                            minnumstr = minnumstr, iter = iter, pops = pops, 
                            mut_chance = mut_chance, elitism_rate = elitism_rate, 
                            highvalue = 1e+08, suggestions = suggestions, realAllocation = TRUE, 
                            writeFiles = writeFiles, showPlot = showPlot, parallel = parallel, 
                            cores = cores)
    newstrata <- updateStrata(strata, solut)
    framenew <- updateFrame(frame = framesamp, newstrata = newstrata)
    if (!is.null(framecens)) {
      colnames(framecens) <- toupper(colnames(framecens))
      framenew$STRATUM <- as.character(framenew$STRATUM)
      framecens$LABEL <- 99999
      framecens$STRATUM <- "99999"
      framenew <- rbind(framenew, framecens)
      framenew$STRATUM <- as.numeric(framenew$LABEL)
      nvarX <- length(grep("X", colnames(framecens)))
      for (i in c(1:nvarX)) {
        st <- paste("framecens$X", i, " <- NULL", sep = "")
        eval(parse(text = st))
      }
      framecens$X1 <- 99999
      cens <- buildStrataDF(framecens, progress = FALSE)
      cens$X1 <- NULL
      cens$SOLUZ <- cens$N
      cens$CENS <- 1
      solut$aggr_strata <- rbind(solut$aggr_strata, cens)
    }
    solution <- list(indices = solut$indices, framenew = framenew, 
                     aggr_strata = solut$aggr_strata)
  }
  if (method == "continuous") {
    if (is.null(errors)) 
      stop("The 'precision constraints' (errors) dataframe is missing")
    if (is.null(framesamp)) 
      stop("The 'sampling frame' (framesamp) dataframe is missing")
    checkInput(errors, sampframe = framesamp)
    check_consecutive_domainvalue(framesamp, "framesamp")
    check_consecutive_domainvalue(errors, "errors")
    if (!is.null(framecens)) 
      checkInput(errors, sampframe = framecens)
    if (!is.null(framecens)) 
      strcens <- TRUE
    if (!is.na(fitting[1])) 
      stop("Fitting value(s) not required with this method")
    if (!is.na(range[1])) 
      stop("Range value(s) not required with this method")
    if (!is.na(kappa)) 
      stop("Kappa value not required with this method")
    solution <- optimizeStrata2(errors = errors, framesamp = framesamp, 
                                framecens = framecens, strcens = strcens, model = model, 
                                alldomains = alldomains, dom = dom, nStrata = nStrata, 
                                minnumstr = minnumstr, iter = iter, pops = pops, 
                                mut_chance = mut_chance, elitism_rate = elitism_rate, 
                                highvalue = 1e+08, suggestions = suggestions, realAllocation = TRUE, 
                                writeFiles = writeFiles, showPlot = showPlot, parallel = parallel, 
                                cores = cores)
  }
  if (method == "spatial") {
    check_spatial_frame(framesamp, errors, nStrata, minnumstr)
    checkInput(errors, sampframe = framesamp)
    check_consecutive_domainvalue(framesamp, "framesamp")
    check_consecutive_domainvalue(errors, "errors")
    nvarY <- length(grep("Y", names(framesamp)))
    if (is.na(fitting[1])) 
      stop("Fitting values of spatial models must be given")
    if (is.na(range[1])) 
      stop("Range values of spatial models must be given")
    if (is.na(kappa)) 
      kappa <- 3
    if (nvarY != length(as.numeric(fitting))) 
      stop("Fitting values must be equal to the number of Y's")
    if (nvarY != length(as.numeric(range))) 
      stop("Range values must be equal to the number of Y's")
    nvars <- length(grep("var", names(framesamp)))
    if (nvarY != nvars) 
      stop("Variances in the 'framesamp' dataframe must be given (one for each Y)")
    for (i in (1:nvars)) {
      stmt <- paste("if (min(framesamp$var", i, ") < 0) stop('Variance var", 
                    i, " of variable Y", i, " has negative values in framesamp')", 
                    sep = "")
      eval(parse(text = stmt))
    }
    if (sum(grep("lon", colnames(framesamp))) == 0 | sum(grep("lat", 
                                                              colnames(framesamp))) == 0) 
      stop("Coordinates (lon and lat) must be given in 'framesamp' dataframe")
    if (!is.null(framecens)) {
      checkInput(errors, sampframe = framecens)
      strcens <- TRUE
      nvars <- length(grep("var", names(framecens)))
      if (nvarY != nvars) 
        stop("Variances in the 'framecens' dataframe must be given (one for each Y)")
      for (i in (1:nvars)) {
        stmt <- paste("if (min(framecens$var", i, ") < 0) stop('Variance var", 
                      i, " of variable Y", i, " has negative values in framecens')", 
                      sep = "")
        eval(parse(text = stmt))
      }
      if (sum(grep("lon", colnames(framecens))) == 0 | 
          sum(grep("lat", colnames(framecens))) == 0) 
        stop("Coordinates (lon and lat) must be given in 'framecens' dataframe")
    }
    solution <- optimizeStrataSpatial(errors = errors, framesamp = framesamp, 
                                      framecens = framecens, strcens = strcens, alldomains = alldomains, 
                                      dom = dom, nStrata = nStrata, minnumstr = minnumstr, 
                                      iter = iter, pops = pops, mut_chance = mut_chance, 
                                      elitism_rate = elitism_rate, highvalue = 1e+08, 
                                      suggestions = suggestions, realAllocation = TRUE, 
                                      writeFiles = writeFiles, showPlot = showPlot, parallel = parallel, 
                                      cores = cores, fitting = fitting, range = range, 
                                      kappa = kappa)
  }
  return(solution)
}
