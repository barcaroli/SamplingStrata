# ==============================================================================
# File:      optimizeStrata2.R
# Functions: optimizeStrata2() [exported]
#
# Purpose:
#   Optimisation of the stratification with the "continuous" method (only
#   continuous X variables): for each domain the genetic algorithm
#   (strataGenalg2()) searches for the cut points of the X variables that
#   minimise the cost of the Bethel allocation under the precision constraints.
#   Returns the solution, the aggregated strata and the frame labelled with the
#   new strata. Called by optimStrata(method = "continuous").
# ==============================================================================

#' @title Best stratification of a sampling frame for multipurpose surveys (only with continuous stratification variables)
#'
#' @description
#' This function runs a set of other functions to optimise the stratification of a sampling frame, only when stratification variables are of the continuous type. It differentiates from 'optimizeStrata' that accepts both continuous and categorical variables
#'
#' @param errors This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
#'   related to target variables of the survey.
#' @param framesamp This is the (mandatory) dataframe containing the information related to the sampling frame.
#' @param framecens This the (optional) dataframe containing the units to be selected in any case. It has same structure than "frame" dataframe.
#' @param strcens Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
#' @param model In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
#'   to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
#'   target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters
#'   beta, sig2, gamma (> 1 in case of heteroscedasticity).
#'   Default is NULL.
#' @param alldomains Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
#'   then a value must be given to parameter 'dom'.
#' @param dom Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval
#'   (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
#' @param nStrata Indicates the number of strata for each variable.
#' @param minnumstr Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
#' @param iter Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
#' @param pops The dimension of each generations in terms of individuals. Default is 20.
#' @param mut_chance Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
#'   of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
#'   but a slower convergence, while low values permit a faster convergence, but the final solution can be
#'   distant from the optimal one. Default is NA, in correspondence of which it is computed as 1/(vars+1)
#'   where vars is the length of elements in the solution.
#' @param elitism_rate This  parameter indicates the rate of better solutions that must be preserved from one generation
#'   to another. Default is 0.2 (20\%).
#' @param highvalue Parameter for genetic algorithm. In should not be changed
#' @param suggestions Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL.
#' @param realAllocation If FALSE, the allocation is based on INTEGER values; if TRUE, the allocation is based on REAL values.
#'   Default is TRUE.
#' @param writeFiles Indicates if the various dataframes and plots produced during the execution have to be written in the working
#'   directory.
#'   Default is FALSE.
#' @param showPlot Indicates if the plot showing the trend in the value of the objective function has to be shown or not. In parallel = TRUE, this defaults to FALSE
#'   Default is TRUE.
#' @param parallel Should the analysis be run in parallel.
#'   Default is TRUE.
#' @param cores If the analysis is run in parallel, how many cores should be used.
#'   If not specified n-1 of total available cores are used OR if number of domains < (n-1) cores, then number of cores equal to number of domains are used.
#'
#' @return
#' A list containing
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
#' ################################
#' # Example of "continuous" method
#' ################################
#' kmean <- KmeansSolution2(frame = frame,
#'                          errors = cv,
#'                          maxclusters = 10)
#' nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
#'                  FUN=function(x) length(unique(x)))
#' solution <- optimizeStrata2(framesamp = frame,
#'                         errors = cv,
#'                         nStrata = nstrat,
#'                         iter = 50,
#'                         pops = 10)
#' framenew <- solution$framenew
#' outstrata <- solution$aggr_strata
#' s <- selectSample(framenew,outstrata)
#' }
#'
#' @keywords survey
#' @export
optimizeStrata2 <- function (errors, framesamp, framecens = NULL, strcens = FALSE, 
                             model = NULL, alldomains = TRUE, dom = NULL, nStrata = c(5), 
                             minnumstr = 2, iter = 50, pops = 20, mut_chance = NA, elitism_rate = 0.2, 
                             highvalue = 1e+08, suggestions = NULL, realAllocation = TRUE, 
                             writeFiles = FALSE, showPlot = TRUE, parallel = TRUE, cores = NA) 
{
  if (strcens == FALSE) {
    cens = NULL
    censi = NULL
  }
  frame <- framesamp
  if (alldomains == TRUE) 
    dom <- NULL
  colnames(frame) <- toupper(colnames(frame))
  ndom <- length(unique(frame$DOMAINVALUE))
  # a single number of strata is valid for every domain
  if (length(nStrata) == 1) nStrata <- rep(nStrata, ndom)
  nvarX <- length(grep("X", names(frame)))
  if (alldomains == FALSE) {
    frame <- frame[frame$DOMAINVALUE == dom, ]
  }
  if (strcens == TRUE & is.null(framecens)) 
    stop("No data in the cens dataframe")
  if (strcens == TRUE) {
    censi <- NULL
    if (alldomains == FALSE) {
      colnames(framecens) <- toupper(colnames(framecens))
      framecens <- framecens[framecens$DOMAINVALUE == 
                               dom, ]
      if (nrow(framecens) == 0) {
        censi <- NULL
        cens <- NULL
      }
      if (nrow(framecens) > 0) {
        censi <- framecens
      }
    }
    if (nrow(framecens) > 0) {
      colnames(framecens) <- toupper(colnames(framecens))
      framecensold <- framecens
      framecens$X1 <- max(nStrata) + 1
      nvarX <- length(grep("X", names(framecens)))
      if (nvarX > 1) {
        for (i in (2:nvarX)) {
          eval(parse(text = paste("framecens$X", i, 
                                  " <- NULL", sep = "")))
        }
      }
      cens <- buildStrataDF(framecens, progress = FALSE, 
                            verbose = FALSE)
      cens$CENS <- 1
      censtot <- cens
    }
  }
  if (writeFiles == TRUE) {
    dire <- getwd()
    direnew <- paste(dire, "/output", sep = "")
    # the folder is no longer deleted: files of a previous run are overwritten
    if (!dir.exists(direnew)) 
      dir.create(direnew)
  }
  if (parallel == TRUE & ndom == 1) 
    parallel <- FALSE
  if (parallel == FALSE & !missing(cores) & !is.na(cores)) {
    cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
    cores <- 1
    Sys.sleep(0.314)
  }
  if (alldomains == FALSE & (parallel == TRUE | !missing(cores) | 
                             !is.na(cores))) {
    cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
    cores <- 1
    Sys.sleep(0.314)
  }
  colnames(errors) <- toupper(colnames(errors))
  colnames(frame) <- toupper(colnames(frame))
  errors$DOMAINVALUE <- as.factor(errors$DOMAINVALUE)
  erro <- split(errors, list(errors$DOMAINVALUE))
  stcamp <- split(frame, list(frame$DOMAINVALUE))
  domain_ids <- unique(frame$DOMAINVALUE)
  if (!is.null(suggestions)) {
    for (i in (unique(frame$DOMAINVALUE))) {
      nvalues <- nrow(suggestions[suggestions$domainvalue == 
                                    i, ])
      if (nvalues != nvarX * (nStrata[i] - 1)) 
        stop("Number of values in suggestions not compatible with nStrata")
    }
    suggestdom <- split(suggestions, list(suggestions$domainvalue))
  }
  if (strcens == TRUE & !is.null(cens) > 0) {
    colnames(cens) <- toupper(colnames(cens))
    stcens <- vector("list", length(domain_ids))
    names(stcens) <- as.character(domain_ids)
    for (domain_id in domain_ids) {
      stcens[[as.character(domain_id)]] <- cens[cens$DOM1 == domain_id, , drop = FALSE]
    }
  }
  if (alldomains == TRUE) {
    if (ndom > length(nStrata)) 
      stop("'Number of strata (nStrata) ' vector lenght (=", 
           length(nStrata), ") \nis not compatible with number of domains (=", 
           ndom, ")\nSet nStrata with a number of elements equal to the number of domains")
    vettsol <- NULL
    outstrata <- NULL
    if (parallel) {
      if (missing(cores) | is.na(cores)) {
        cores <- parallel::detectCores() - 1
        if (ndom < cores) 
          cores <- ndom
      }
      if (cores < 2) 
        stop("\nOnly one core available: no parallel processing possible. \n               \nPlease change parameter parallel to FALSE and run again")
      cat("\n *** Starting parallel optimization for ", 
          ndom, " domains using ", cores, " cores\n")
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, ls(), envir = environment())
      showPlot <- FALSE
      par_ga_sol = pblapply(cl = cl, X = domain_ids, 
                            FUN = function(i) {
                              erro[[i]] <- erro[[i]][, -ncol(errors)]
                              flagcens <- strcens
                              cens_i <- NULL
                              if (strcens == TRUE && exists("stcens")) {
                                cens_i <- stcens[[as.character(i)]]
                              }
                              if (strcens == TRUE) {
                                if (!is.null(cens_i) && nrow(cens_i) > 0) {
                                  cens <- cens_i
                                  flagcens = TRUE
                                }
                              }
                              if (strcens == TRUE) {
                                if (is.null(cens_i) || nrow(cens_i) == 0) {
                                  cens <- NULL
                                  flagcens = FALSE
                                }
                              }
                              if (strcens == FALSE) {
                                cens <- NULL
                                flagcens = FALSE
                              }
                              if (!is.null(suggestions) & alldomains == 
                                  TRUE) {
                                suggest <- matrix(0, nrow = 2, ncol = (nStrata[i] - 
                                                                         1) * nvarX)
                                suggest[1, ] <- suggestdom[[i]]$suggestions1
                                suggest[2, ] <- suggestdom[[i]]$suggestions2
                              }
                              if (!is.null(suggestions) & alldomains == 
                                  FALSE) {
                                suggest <- matrix(0, nrow = 2, ncol = (nStrata[i] - 
                                                                         1) * nvarX)
                                suggest[1, ] <- suggestdom[[i]]$suggestions1
                                suggest[2, ] <- suggestdom[[i]]$suggestions2
                              }
                              if (is.null(suggestions)) {
                                suggest <- NULL
                              }
                              if (nrow(stcamp[[i]]) > 0) {
                                solut <- strataGenalg2(errors = erro[[i]], 
                                                       frame = stcamp[[i]], cens = cens, strcens = flagcens, 
                                                       model, ncuts = (nStrata[i] - 1), dominio = i, 
                                                       minnumstr, iter, pops, mut_chance, elitism_rate, 
                                                       suggestions = suggest, realAllocation, 
                                                       writeFiles, showPlot)
                                rbga.results <- solut[[3]]
                                if (nrow(stcamp[[i]]) == 1) {
                                  solut <- list(c(1), stcamp[[i]][c(1:grep("DOMAINVALUE", 
                                                                           colnames(stcamp[[i]])))])
                                  error <- data.frame(erro[[i]], stringsAsFactors = TRUE)
                                  strat <- data.frame(solut[[2]], stringsAsFactors = TRUE)
                                  solut[[2]]$SOLUZ <- sum(bethel_cpp(strat, 
                                                                     error, realAllocation = T))
                                  if (solut[[2]]$SOLUZ > solut[[2]]$N) 
                                    solut[[2]]$SOLUZ <- solut[[2]]$N
                                }
                                vettsol <- solut[[1]]
                                if (length(outstrata) > 0) 
                                  colnames(outstrata) <- toupper(colnames(outstrata))
                                colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
                                outstrata <- solut[[2]]
                                list(vettsol = vettsol, outstrata = outstrata, 
                                     rbga.results = rbga.results)
                              }
                            })
      vettsol <- do.call(rbind, lapply(par_ga_sol, `[[`, 
                                       1))
      outstrata <- do.call(rbind, lapply(par_ga_sol, `[[`, 
                                         2))
      results <- do.call(rbind, lapply(par_ga_sol, `[[`, 
                                       3))
      for (i in domain_ids) {
        rbga.object <- par_ga_sol[[i]]$rbga.results
        max <- max(rbga.object$best, rbga.object$mean)
        min <- min(rbga.object$best, rbga.object$mean)
        if (writeFiles == TRUE) {
          stmt <- paste("png(filename = file.path(direnew, 'plotdom", 
                        i, ".png'),height=5, width=7, units='in', res=144)", 
                        sep = "")
          eval(parse(text = stmt))
        }
        plot(rbga.object$best, type = "l", main = "", 
             ylim = c(min, max), xlab = "Iteration (Generation)", 
             ylab = "Best (black lower line) and mean (red upper line) evaluation value")
        lines(rbga.object$mean, col = "red")
        title(paste("Domain #", i, " - Sample cost", 
                    round(min(rbga.object$best), 2)), col.main = "red")
        if (writeFiles == TRUE) 
          dev.off()
      }
    }
    else {
      for (i in domain_ids) {
        cat("\n-------------------------------------------------------")
        cat("\n *** Domain : ", i)
        if (nrow(stcamp[[i]]) > 0) 
          cat("\n Number of strata : ", nrow(stcamp[[i]]))
        erro[[i]] <- erro[[i]][, -ncol(errors)]
        flagcens <- strcens
        cens_i <- NULL
        if (strcens == TRUE && exists("stcens")) {
          cens_i <- stcens[[as.character(i)]]
        }
        if (strcens == TRUE) {
          if (!is.null(cens_i) && nrow(cens_i) > 0) {
            cens <- cens_i
            flagcens = TRUE
          }
        }
        if (strcens == TRUE) {
          if (is.null(cens_i) || nrow(cens_i) == 0) {
            cens <- NULL
            flagcens = FALSE
          }
        }
        if (strcens == FALSE) {
          cens <- NULL
          flagcens = FALSE
        }
        if (!is.null(suggestions) & alldomains == TRUE) {
          suggest <- matrix(0, nrow = 2, ncol = (nStrata[i] - 
                                                   1) * nvarX)
          suggest[1, ] <- suggestdom[[i]]$suggestions1
          suggest[2, ] <- suggestdom[[i]]$suggestions2
        }
        if (!is.null(suggestions) & alldomains == FALSE) {
          suggest <- matrix(0, nrow = 2, ncol = (nStrata[i] - 
                                                   1) * nvarX)
          suggest[1, ] <- suggestions$suggestions1
          suggest[2, ] <- suggestions$suggestions2
        }
        if (is.null(suggestions)) {
          suggest <- NULL
        }
        if (nrow(stcamp[[i]]) > 0) {
          solut <- strataGenalg2(errors = erro[[i]], 
                                 frame = stcamp[[i]], cens = cens, strcens = flagcens, 
                                 model, ncuts = (nStrata[i] - 1), dominio = i, 
                                 minnumstr, iter, pops, mut_chance, elitism_rate, 
                                 suggestions = suggest, realAllocation, writeFiles, 
                                 showPlot)
          if (nrow(stcamp[[i]]) == 1) {
            solut <- list(c(1), stcamp[[i]][c(1:grep("DOMAINVALUE", 
                                                     colnames(stcamp[[i]])))])
            error <- data.frame(erro[[i]], stringsAsFactors = TRUE)
            strat <- data.frame(solut[[2]], stringsAsFactors = TRUE)
            solut[[2]]$SOLUZ <- sum(bethel_cpp(strat, 
                                               error, realAllocation = T))
            if (solut[[2]]$SOLUZ > solut[[2]]$N) 
              solut[[2]]$SOLUZ <- solut[[2]]$N
          }
          vettsol <- rbind(vettsol, solut[[1]])
          if (length(outstrata) > 0) 
            colnames(outstrata) <- toupper(colnames(outstrata))
          colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
          outstrata <- rbind(outstrata, solut[[2]])
          rbga.object <- solut[[3]]
          max <- max(rbga.object$best, rbga.object$mean)
          min <- min(rbga.object$best, rbga.object$mean)
          if (writeFiles == TRUE) {
            stmt <- paste("png(filename = file.path(direnew, 'plotdom", 
                          i, ".png'),height=5, width=7, units='in', res=144)", 
                          sep = "")
            eval(parse(text = stmt))
          }
          plot(rbga.object$best, type = "l", main = "", 
               ylim = c(min, max), xlab = "Iteration (Generation)", 
               ylab = "Best (black lower line) and mean (red upper line) evaluation value")
          lines(rbga.object$mean, col = "red")
          title(paste("Domain #", i, " - Sample cost", 
                      round(min(rbga.object$best), 2)), col.main = "red")
          if (writeFiles == TRUE) 
            dev.off()
        }
      }
    }
  }
  if (alldomains == FALSE) {
    i <- dom
    erro[[i]] <- erro[[i]][, -ncol(errors)]
    flagcens <- strcens
    if (strcens == TRUE) {
      flagcens <- TRUE
      if (!is.null(cens)) 
        colnames(cens) <- toupper(colnames(cens))
      censi <- cens[cens$DOM1 == i, ]
      if (nrow(censi) == 0) {
        flagcens <- FALSE
      }
    }
    if (!is.null(suggestions)) {
      # the two chromosomes of the selected domain (suggestions may contain
      # all the domains, as returned by prepareSuggestion())
      sugdom <- suggestions[suggestions$domainvalue == i, , drop = FALSE]
      suggest <- matrix(0, nrow = 2, ncol = (nStrata[i] - 1) * nvarX)
      suggest[1, ] <- sugdom$suggestions1
      suggest[2, ] <- sugdom$suggestions2
    }
    if (is.null(suggestions)) {
      suggest <- NULL
    }
    solut <- strataGenalg2(errors = erro[[i]], frame = frame, 
                           cens = cens, strcens = flagcens, model, ncuts = (nStrata[i] - 1), dominio = i, minnumstr, iter, pops, mut_chance, 
                           elitism_rate, suggestions = suggest, realAllocation, 
                           writeFiles, showPlot)
    vettsol <- solut[[1]]
    outstrata <- solut[[2]]
    rbga.object <- solut[[3]]
    if (writeFiles == TRUE) {
      stmt <- paste("png(filename = file.path(direnew, 'plotdom", 
                    dom, ".png'),height=5, width=7, units='in', res=144)", 
                    sep = "")
      eval(parse(text = stmt))
    }
    max <- max(rbga.object$best, rbga.object$mean)
    min <- min(rbga.object$best, rbga.object$mean)
    plot(rbga.object$best, type = "l", main = "", ylim = c(min, 
                                                           max), xlab = "Iteration (Generation)", ylab = "Best (black lower line) and mean (red upper line) evaluation value")
    lines(rbga.object$mean, col = "red")
    title(paste("Domain # ", dom, " - Sample cost", round(min(rbga.object$best), 
                                                          2)), col.main = "red")
    if (writeFiles == TRUE) 
      dev.off()
  }
  colnames(outstrata) <- toupper(colnames(outstrata))
  dimens <- sum(round(outstrata$SOLUZ))
  cat("\n *** Sample size : ", dimens)
  cat("\n *** Number of strata : ", nrow(outstrata))
  cat("\n---------------------------")
  if (writeFiles == TRUE) {
    write.table(outstrata, file = file.path(direnew, "outstrata.txt"), 
                sep = "\t", row.names = FALSE, col.names = TRUE, 
                quote = FALSE)
    cat("\n...written output to ", direnew, "/outstrata.txt\n")
  }
  vettsoldf <- as.data.frame(vettsol, stringsAsFactors = TRUE)
  colnames(vettsoldf) <- c("ID", "LABEL")
  vettsoldf$STRATO <- vettsoldf$LABEL
  framenew <- merge(frame, vettsoldf, by = c("ID"))
  if (strcens == TRUE) {
    framecens$STRATO <- nStrata[i] + 1
    framecens$LABEL <- nStrata[i] + 1
    if (alldomains == FALSE) {
      framecens <- framecens[framecens$DOMAINVALUE == 
                               dom, ]
      framecensold <- framecensold[framecensold$DOMAINVALUE == 
                                     dom, ]
      framecens$STRATO <- nStrata + 1
      framecens$LABEL <- nStrata + 1
    }
    for (i in (1:nvarX)) {
      eval(parse(text = paste("framecens$X", i, " <- framecensold$X", 
                              i, sep = "")))
    }
    framecens$STRATO <- as.character(max(as.numeric(outstrata$STRATO)) + 
                                       1)
    framecens$LABEL <- max(as.numeric(outstrata$STRATO)) + 
      1
    censtot$STRATO <- as.character(max(as.numeric(outstrata$STRATO)) + 
                                     1)
    colnames(framecens) <- toupper(colnames(framecens))
    framenew <- rbind(framenew, framecens)
    censtot$SOLUZ <- censtot$N
    outstrata <- rbind(outstrata, censtot)
  }
  solution <- list(indices = vettsol, aggr_strata = outstrata, 
                   framenew = framenew)
  if (writeFiles == TRUE) {
    setwd(dire)
  }
  return(solution)
}
