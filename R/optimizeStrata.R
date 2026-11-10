# ==============================================================================
# File:      optimizeStrata.R
# Functions: optimizeStrata() [exported]
#
# Purpose:
#   Optimisation of the stratification with the "atomic" method: the atomic
#   strata (cross-classification of categorical X variables) of each domain are
#   aggregated by a genetic algorithm (strataGenalg()) searching for the
#   partition that minimises the cost of the Bethel allocation under the
#   precision constraints. Can run the domains in parallel. Called by
#   optimStrata(method = "atomic").
# ==============================================================================

#' @title Best stratification of a sampling frame for multipurpose surveys
#'
#' @description
#' This function runs a set of other functions to optimise the stratification of a sampling frame
#'
#' @param errors This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
#'   related to target variables of the survey.
#' @param strata This is the (mandatory) dataframe containing the information related to "atomic" strata, i.e. the strata obtained by
#'   the Cartesian product of all auxiliary variables X's. Information concerns the identifiability of strata
#'   (values of X's) and variability of Y's (for each Y, mean and standard error in strata).
#' @param cens This the (optional) dataframe containing the takeall strata, those strata whose units must be selected in whatever sample. It has same structure than "strata" dataframe.
#' @param strcens Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
#' @param alldomains Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
#'   then a value must be given to parameter 'dom'.
#' @param dom Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval
#'   (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
#' @param initialStrata This is the initial limit on the number of strata in the different domains for each solution. Default is NA, and in this case it is set equal to the
#'   number of atomic strata in each domain.
#' @param addStrataFactor This parameter indicates the probability that at each mutation the number of strata may increase with
#'   respect to the current value. Default is 0.0.
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
#' A list containing (1) the vector of the solution and (2) the optimal aggregated strata
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
#' solution <- optimizeStrata(strata = strata,
#'                         errors = cv,
#'                         initialStrata = nstrat,
#'                         suggestions = kmean,
#'                         iter = 50,
#'                         pops = 10)
#' outstrata <- solution$aggr_strata
#' newstrata <- updateStrata(strata,solution)
#' framenew <- updateFrame(frame, newstrata)
#' s <- selectSample(framenew, outstrata)
#' }
#'
#' @keywords survey
#' @export
optimizeStrata <- 
  function (errors, strata, cens = NULL, strcens = FALSE, alldomains = TRUE, 
            dom = NULL, initialStrata = NA, addStrataFactor = 0, minnumstr = 2, 
            iter = 50, pops = 20, mut_chance = NA, elitism_rate = 0.2, 
            highvalue = 1e+08, suggestions = NULL, realAllocation = TRUE, 
            writeFiles = FALSE, showPlot = TRUE, parallel = TRUE, cores = NA) 
  {
    if (writeFiles == TRUE) {
      dire <- getwd()
      direnew <- paste(dire, "/output", sep = "")
      # the folder is no longer deleted: files of a previous run are overwritten
      if (!dir.exists(direnew)) 
        dir.create(direnew)
      #setwd(direnew)
    }
    
    if (parallel == FALSE & !missing(cores) & !is.na(cores)) {
      cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
      cores <- 1
      Sys.sleep(0.314)
    }
      
    if (is.na(initialStrata[1])) 
      initialStrata <- as.numeric(table(strata$DOM1))
    nstrata = initialStrata
    colnames(errors) <- toupper(colnames(errors))
    colnames(strata) <- toupper(colnames(strata))
    errors$DOMAINVALUE <- as.factor(errors$DOMAINVALUE)
    erro <- split(errors, list(errors$DOMAINVALUE))
    stcamp <- split(strata, list(strata$DOM1))
    if (!is.null(suggestions)) 
      suggestdom <- split(suggestions, list(suggestions$domainvalue))
    if (strcens == TRUE) {
      colnames(cens) <- toupper(colnames(cens))
      k <- length(levels(as.factor(strata$DOM1)))
      stcens <- NULL
      for (i in (1:k)) {
        stcens[[i]] <- cens[cens$DOM1 == i, ]
      }
    }
    ndom <- length(levels(as.factor(strata$DOM1)))
    if (alldomains == TRUE) {
      if (ndom > length(nstrata)) 
        stop("'initialStrata' vector lenght (=", length(nstrata), 
             ") \nis not compatible with number of domains (=", 
             ndom, ")\nSet initialStrata with a number of elements equal to the number of domains")
      vettsol <- NULL
      outstrata <- NULL
      if (parallel == TRUE & ndom == 1) parallel <- FALSE
      if (parallel) {
        if (missing(cores) | is.na(cores)) {
          cores <- parallel::detectCores() - 1
          if (ndom < cores) 
            cores <- ndom
        }
        if (cores < 2) 
          stop("\nOnly one core available: no parallel processing possible. 
               \nPlease change parameter parallel = FALSE and run again")
        cat("\n *** Starting parallel optimization for ", 
            ndom, " domains using ", cores, " cores\n")
        cl <- parallel::makePSOCKcluster(cores)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterExport(cl = cl, ls(), envir = environment())
        
        showPlot <- FALSE
        
        par_ga_sol = pblapply(
          cl = cl, X = unique(strata$DOM1), FUN = function(i)  {         
                                erro[[i]] <- erro[[i]][, -ncol(errors)]
                                cens <- NULL
                                flagcens <- strcens
                                if (strcens == TRUE) {
                                  if (nrow(stcens[[i]]) > 0) {
                                    cens <- stcens[[i]]
                                    flagcens = TRUE
                                  }
                                }
                                if (strcens == TRUE) {
                                  if (nrow(stcens[[i]]) == 0) {
                                    cens <- NULL
                                    flagcens = FALSE
                                  }
                                }
                                if (strcens == FALSE) {
                                  cens <- NULL
                                  flagcens = FALSE
                                }
                                if (!is.null(suggestions) & alldomains == TRUE) {
                                  suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
                                  suggest[1, ] <- suggestdom[[i]]$suggestions
                                }
                                if (!is.null(suggestions) & alldomains == FALSE) {
                                  suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
                                  suggest[1, ] <- suggestions$suggestions
                                }
                                if (is.null(suggestions)) {
                                  suggest <- NULL
                                }
                                if (nrow(stcamp[[i]]) > 0) {
                                  solut <- strataGenalg(errors = erro[[i]], 
                                                                         strata = stcamp[[i]], cens = cens, strcens = flagcens, 
                                                                         dominio = i, initialStrata = nstrata[i], 
                                                                         minnumstr, iter, pops, mut_chance, elitism_rate, 
                                                                         addStrataFactor, highvalue, suggestions = suggest, 
                                                                         realAllocation, writeFiles, showPlot)
                                  if (nrow(stcamp[[i]]) == 1) {
                                    solut <- list(c(1), stcamp[[i]][c(1:grep("DOM1", 
                                                                             colnames(stcamp[[i]])))])
                                    error <- data.frame(erro[[i]],stringsAsFactors = TRUE)
                                    strat <- data.frame(solut[[2]],stringsAsFactors = TRUE)
                                    solut[[2]]$SOLUZ <- sum(bethel_cpp(strat, error, 
                                                                   realAllocation = T, minnumstrat = minnumstr))
                                    if (solut[[2]]$SOLUZ > solut[[2]]$N) 
                                      solut[[2]]$SOLUZ <- solut[[2]]$N
                                  }
                                  vettsol <- solut[[1]]
                                  if (length(outstrata) > 0) 
                                    colnames(outstrata) <- toupper(colnames(outstrata))
                                  colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
                                  outstrata <- solut[[2]]
                                  rbga.results <- solut[[3]]
                                  list(vettsol = vettsol, outstrata = outstrata, 
                                       rbga.results = rbga.results)
                                }
                              }
        )
        
        vettsol <- do.call(c, lapply(par_ga_sol, `[[`, 1))
        outstrata <- do.call(rbind, lapply(par_ga_sol, `[[`, 2))
        
        for (i in (unique(strata$DOM1))) {
        # for (i in (seq_along(unique(strata$DOM1)))) {
          rbga.object <- par_ga_sol[[i]]$rbga.results
          max <- max(rbga.object$best, rbga.object$mean)
          min <- min(rbga.object$best, rbga.object$mean)
          if (writeFiles == TRUE) {
            stmt <- paste("png(filename = file.path(direnew, 'plotdom", i, ".png'),height=5, width=7, units='in', res=144)", sep = "")
            eval(parse(text = stmt))
          }  
          plot(rbga.object$best, type = "l", 
               main = "", ylim = c(min,max), 
               xlab = "Iteration (Generation)", 
               ylab = "Best (black lower line) and mean (red upper line) evaluation value")
          lines(rbga.object$mean, col = "red")
          title(paste("Domain #", i, " - Sample cost", 
                      round(min(rbga.object$best), 2)), 
                col.main = "red")
          if (writeFiles == TRUE) dev.off()
        }
      }
      else {
        for (i in (seq_along(unique(strata$DOM1)))) {
          cat("\n *** Domain : ", i, " ", as.character(errors$DOMAINVALUE[i]))
          cat("\n Number of strata : ", nrow(stcamp[[i]]))
          erro[[i]] <- erro[[i]][, -ncol(errors)]
          cens <- NULL
          flagcens <- strcens
          if (strcens == TRUE) {
            if (nrow(stcens[[i]]) > 0) {
              cens <- stcens[[i]]
              flagcens = TRUE
            }
          }
          if (strcens == TRUE) {
            if (nrow(stcens[[i]]) == 0) {
              cens <- NULL
              flagcens = FALSE
            }
          }
          if (strcens == FALSE) {
            cens <- NULL
            flagcens = FALSE
          }
          if (!is.null(suggestions) & alldomains == TRUE) {
            suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
            suggest[1, ] <- suggestdom[[i]]$suggestions
          }
          if (!is.null(suggestions) & alldomains == FALSE) {
            suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
            suggest[1, ] <- suggestions$suggestions
          }
          if (is.null(suggestions)) {
            suggest <- NULL
          }
          if (nrow(stcamp[[i]]) > 0) {
            solut <- SamplingStrata:::strataGenalg(errors = erro[[i]], 
                                  strata = stcamp[[i]], cens = cens, strcens = flagcens, 
                                  dominio = i, initialStrata = nstrata[i], 
                                  minnumstr, iter, pops, mut_chance, elitism_rate, 
                                  addStrataFactor, highvalue, suggestions = suggest, 
                                  realAllocation, writeFiles, showPlot)
            rbga.object <- solut[[3]]
            if (nrow(stcamp[[i]]) == 1) {
              solut <- list(c(1), stcamp[[i]][c(1:grep("DOM1", 
                                                       colnames(stcamp[[i]])))])
              error <- data.frame(erro[[i]],stringsAsFactors = TRUE)
              strat <- data.frame(solut[[2]],stringsAsFactors = TRUE)
              solut[[2]]$SOLUZ <- sum(SamplingStrata:::bethel_cpp(strat, error, 
                                             realAllocation = T, minnumstrat = minnumstr))
              if (solut[[2]]$SOLUZ > solut[[2]]$N) 
                solut[[2]]$SOLUZ <- solut[[2]]$N
            }
            vettsol <- c(vettsol, solut[[1]])
            if (length(outstrata) > 0) 
              colnames(outstrata) <- toupper(colnames(outstrata))
            colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
            outstrata <- rbind(outstrata, solut[[2]])

            max <- max(rbga.object$best, rbga.object$mean)
            min <- min(rbga.object$best, rbga.object$mean)
            if (writeFiles == TRUE) {
              stmt <- paste("png(filename = file.path(direnew, 'plotdom", i, ".png'),height=5, width=7, units='in', res=144)", sep = "")
              eval(parse(text = stmt))
            }  
            plot(rbga.object$best, type = "l", 
                 main = "", ylim = c(min,max), 
                 xlab = "Iteration (Generation)", 
                 ylab = "Best (black lower line) and mean (red upper line) evaluation value")
            lines(rbga.object$mean, col = "red")
            title(paste("Domain #", i, " - Sample cost", 
                        round(min(rbga.object$best), 2)), 
                  col.main = "red")
            if (writeFiles == TRUE) dev.off()
          }
        }
      }
    }
    if (alldomains == FALSE) {
      if (dom < 1 | dom > ndom) 
        stop("\nInvalid value of the indicated domain\n")
      i <- dom
      erro[[i]] <- erro[[i]][, -ncol(errors)]
      flagcens <- strcens
      if (strcens == TRUE) {
        flagcens <- TRUE
        colnames(cens) <- toupper(colnames(cens))
        censi <- cens[cens$DOM1 == i, ]
        if (nrow(censi) == 0) {
          flagcens <- FALSE
          censi <- NULL
        }
      }
      if (!is.null(suggestions) & alldomains == TRUE) {
        suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
        suggest[1, ] <- suggestdom[[i]]$suggestions
      }
      if (!is.null(suggestions) & alldomains == FALSE) {
        suggest <- matrix(0, nrow = 1, ncol = nrow(stcamp[[i]]))
        suggest[1, ] <- suggestdom[[dom]]$suggestions
      }
      if (is.null(suggestions)) {
        suggest <- NULL
      }
      solut <- strataGenalg(errors = erro[[i]], strata = stcamp[[i]], 
                            cens = censi, strcens = flagcens, dominio = i, initialStrata = nstrata[i], 
                            minnumstr, iter, pops, mut_chance, elitism_rate, 
                            addStrataFactor, highvalue, suggestions = suggest, 
                            realAllocation, writeFiles, showPlot)
      vettsol <- solut[[1]]
      outstrata <- solut[[2]]
      rbga.object <- solut[[3]]
      if (writeFiles == TRUE) {
        stmt <- paste("png(filename = file.path(direnew, 'plotdom",dom,".png'),height=5, width=7, units='in', res=144)", sep = "")
        eval(parse(text = stmt))
      }  
      max <- max(rbga.object$best, rbga.object$mean)
      min <- min(rbga.object$best, rbga.object$mean)
      plot(rbga.object$best, type = "l", 
           main = "", ylim = c(min,max), 
           xlab = "Iteration (Generation)", 
           ylab = "Best (black lower line) and mean (red upper line) evaluation value")
      lines(rbga.object$mean, col = "red")
      title(paste("Domain # ",dom," - Sample cost", 
                  round(min(rbga.object$best), 2)), 
            col.main = "red")
      if (writeFiles == TRUE) dev.off()
    }
    colnames(outstrata) <- toupper(colnames(outstrata))
    dimens <- sum(round(outstrata$SOLUZ))
    cat("\n *** Sample size : ", dimens)
    cat("\n *** Number of strata : ", nrow(outstrata))
    cat("\n---------------------------")
    if (writeFiles == TRUE) {
      write.table(outstrata, file = file.path(direnew, "outstrata.txt"), sep = "\t", 
                  row.names = FALSE, col.names = TRUE, quote = FALSE)
      cat("\n...written output to ", direnew,"/outstrata.txt\n")
    }
    solution <- list(indices = vettsol, aggr_strata = outstrata)
    if (writeFiles == TRUE) {
      setwd(dire)
    }
    return(solution)
  }
