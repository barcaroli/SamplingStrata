# ==============================================================================
# File:      optimStrata.R
# Functions: optimStrata() [exported]
#
# Purpose:
#   Main entry point for the optimisation of the stratification. After checking
#   the inputs (consecutive domain codes, structure of the spatial frame), it
#   calls optimizeStrata() (method = "atomic"), optimizeStrata2() (method =
#   "continuous") or optimizeStrataSpatial() (method = "spatial") and returns
#   their results in a common format.
# ==============================================================================

#' @title Optimization of the stratification of a sampling frame given a sample survey
#'
#' @description
#' Wrapper function to call the different optimization functions: (i) optimizeStrata (method = "atomic"); (ii) optimizeStrata2 (method = "continuous"); (iii) optimizeStrataSpatial (method = "spatial").
#' For continuity reasons, these functions are still available to be used standalone.
#'
#' @param method This parameter allows to choose the method to be applied in the optimization step: (i) optimizeStrata (method = "atomic"); (ii) optimizeStrata (method = "continuous"); (iii) optimizeStrata (method = "spatial")
#' @param errors This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
#'   related to target variables of the survey.
#' @param framesamp This is the dataframe containing the information related to the sampling frame.
#' @param framecens This the dataframe containing the units to be selected in any case. It has same structure than "framesamp" dataframe.
#' @param nStrata Indicates the number of strata to be obtained in the final solution.
#' @param model In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
#'   to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
#'   target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters
#'   beta, sig2, gamma (> 1 in case of heteroscedasticity).
#'   Default is NULL.
#' @param alldomains Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
#'   then a value must be given to parameter 'dom'.
#' @param dom Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval
#'   (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
#' @param strcens Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
#' @param minnumstr Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
#' @param iter Indicates the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
#' @param pops The dimension of each generations in terms of individuals. Default is 20.
#' @param mut_chance Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
#'   of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
#'   but a slower convergence, while low values permit a faster convergence, but the final solution can be
#'   distant from the optimal one. Default is NA, in correspondence of which it is computed as 1/(vars+1)
#'   where vars is the length of elements in the solution.
#' @param elitism_rate This  parameter indicates the rate of better solutions that must be preserved from one generation
#'   to another. Default is 0.2 (20\%).
#' @param suggestions Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL.
#' @param writeFiles Indicates if the various dataframes and plots produced during the execution have to be written in the working
#'   directory.
#'   Default is FALSE.
#' @param showPlot Indicates if the plot showing the trend in the value of the objective function has to be shown or not. In parallel = TRUE, this defaults to FALSE
#'   Default is TRUE.
#' @param parallel Should the analysis be run in parallel.
#'   Default is TRUE.
#' @param cores If the analysis is run in parallel, how many cores should be used.
#'   If not specified n-1 of total available cores are used OR if number of domains < (n-1) cores, then number of cores equal to number of domains are used.
#' @param fitting Fitting of the model(s) (in terms of R squared).
#'   It can be vector with as many elements as the number of target variables Ys, or a
#'   dataframe with as many rows as the different domains in the frame.
#' @param range Maximum range for spatial autocorrelation.
#'   It can be vector with as many elements as the number of target variables Ys, or a
#'   dataframe with as many rows as the different domains in the frame.
#' @param kappa Factor used in evaluating spatial autocorrelation.
#'
#' @return
#' List containing
#' (1) the vector of the solution,
#' (2) the optimal aggregated strata,
#' (3) the total sampling frame with the label of aggregated strata
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' library(SamplingStrata)
#' ############################
#' # Example of "atomic" method
#' ############################
#' data(swissmunicipalities)
#' swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
#' frame <- buildFrameDF(df = swissmunicipalities,
#'                       id = "id",
#'                       domainvalue = "REG",
#'                       X = c("POPTOT","HApoly"),
#'                       Y = c("Surfacesbois", "Airind"))
#' ndom <- length(unique(frame$domainvalue))
#' cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
#'                          CV1 = rep(0.1,ndom),
#'                          CV2 = rep(0.1,ndom),
#'                          domainvalue = c(1:ndom)))
#' strata <- buildStrataDF(frame)
#' kmean <- KmeansSolution(strata,cv,maxclusters=30)
#' nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
#'                  FUN=function(x) length(unique(x)))
#' sugg <- prepareSuggestion(kmean = kmean,
#'                           frame = frame,
#'                           nStrata = nstrat)
#' evalSuggestion(frame=frame,sugg=sugg,errors=cv)
#' solution <- optimStrata(method ="atomic",
#'                         framesamp = frame,
#'                         errors = cv,
#'                         nStrata = nstrat,
#'                         suggestions = kmean,
#'                         iter = 50,
#'                         pops = 10)
#' outstrata <- solution$aggr_strata
#' framenew <- solution$framenew
#' s <- selectSample(framenew, outstrata)
#' ################################
#' # Example of "continuous" method
#' ################################
#' kmean <- KmeansSolution2(frame = frame,
#'                          errors = cv,
#'                          maxclusters = 10)
#' nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
#'                  FUN=function(x) length(unique(x)))
#' sugg <- prepareSuggestion(kmean = kmean,
#'                           frame = frame,
#'                           nStrata = nstrat)
#' attr(sugg, "seeds")   
#' evalSuggestion(frame, errors = cv, sugg, which = 1)
#' evalSuggestion(frame, errors = cv, sugg, which = 2)
#' solution <- optimStrata(method = "continuous",
#'                         framesamp = frame,
#'                         errors = cv,
#'                         nStrata = nstrat,
#'                         iter = 50,
#'                         pops = 10,
#'                         suggestions = sugg)
#' framenew <- solution$framenew
#' outstrata <- solution$aggr_strata
#' s <- selectSample(framenew,outstrata)
#' #############################
#' # Example of "spatial" method
#' #############################
#' library(sp)
#' data("meuse")
#' data("meuse.grid")
#' meuse.grid$id <- c(1:nrow(meuse.grid))
#' coordinates(meuse) <- c('x','y')
#' coordinates(meuse.grid) <- c('x','y')
#' library(gstat)
#' library(automap)
#' v <- variogram(lead ~ dist + soil, data = meuse)
#' fit.vgm.lead <- autofitVariogram(lead ~ dist + soil,
#'                                  meuse,
#'                                  model = "Exp")
#' plot(v, fit.vgm.lead$var_model)
#' lead.kr <- krige(lead ~ dist + soil, meuse, meuse.grid,
#'                  model = fit.vgm.lead$var_model)
#' lead.pred <- ifelse(lead.kr[1]$var1.pred < 0,0, lead.kr[1]$var1.pred)
#' lead.var <- ifelse(lead.kr[2]$var1.var < 0, 0, lead.kr[2]$var1.var)
#' df <- as.data.frame(list(dom = rep(1,nrow(meuse.grid)),
#'                          lead.pred = lead.pred,
#'                          lead.var = lead.var,
#'                          lon = meuse.grid$x,
#'                          lat = meuse.grid$y,
#'                          id = c(1:nrow(meuse.grid))))
#' frame <-buildFrameSpatial(df = df,
#'                           id = "id",
#'                           X = c("lead.pred"),
#'                           Y = c("lead.pred"),
#'                           variance = c("lead.var"),
#'                           lon = "lon",
#'                           lat = "lat",
#'                           domainvalue = "dom")
#' cv <- as.data.frame(list(DOM = rep("DOM1",1),
#'                          CV1 = rep(0.05,1),
#'                          domainvalue = c(1:1) ))
#' solution <- optimStrata(method = "spatial",
#'                         errors = cv,
#'                         framesamp = frame,
#'                         iter = 25,
#'                         pops = 10,
#'                         nStrata = 5,
#'                         fitting = 1,
#'                         kappa = 1,
#'                         range = fit.vgm.lead$var_model$range[2])
#' framenew <- solution$framenew
#' outstrata <- solution$aggr_strata
#' frameres <- SpatialPixelsDataFrame(points = framenew[c("LON","LAT")],
#'                                    data = framenew)
#' frameres$LABEL <- as.factor(frameres$LABEL)
#' spplot(frameres,c("LABEL"), col.regions=bpy.colors(5))
#' s <- selectSample(framenew,outstrata)
#' }
#'
#' @keywords survey
#' @export
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
  if (!is.character(method) || !all(method %in% c("atomic", "continuous", "spatial"))) 
    stop("Method should be one in 'atomic','continuous','spatial'")
  method <- match.arg(method)
  if (any(nStrata < 2, na.rm = TRUE)) 
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
    fitting_val1 <- if (is.data.frame(fitting)) fitting[1, 1] else fitting[1]
    range_val1   <- if (is.data.frame(range))   range[1, 1]   else range[1]
    if (is.na(fitting_val1))
      stop("Fitting values of spatial models must be given")
    if (is.na(range_val1))
      stop("Range values of spatial models must be given")
    if (is.na(kappa))
      kappa <- 3
    fitting_len <- if (is.data.frame(fitting)) ncol(fitting) else length(as.numeric(fitting))
    range_len   <- if (is.data.frame(range))   ncol(range)   else length(as.numeric(range))
    if (nvarY != fitting_len)
      stop("Fitting values must be equal to the number of Y's")
    if (nvarY != range_len)
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
