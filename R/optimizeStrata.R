optimizeStrata <- 
  function (errors, strata, cens = NULL, strcens = FALSE, alldomains = TRUE, 
            dom = NULL, initialStrata = NA, addStrataFactor = 0, minnumstr = 2, 
            iter = 50, pops = 20, mut_chance = NA, elitism_rate = 0.2, 
            highvalue = 1e+08, suggestions = NULL, realAllocation = TRUE, 
            writeFiles = FALSE, showPlot = TRUE, parallel = TRUE, cores) 
  {
    if (writeFiles == TRUE) {
      dire <- getwd()
      direnew <- paste(dire, "/output", sep = "")
      if (dir.exists(direnew)) 
        unlink(direnew,recursive=TRUE)
      if (!dir.exists(direnew)) 
        dir.create(direnew)
      #setwd(direnew)
    }
    
    if(parallel == FALSE & !missing(cores)){
      cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
      cores <- 1
      Sys.sleep(0.314)
    }
      
    if (is.na(initialStrata)) 
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
        if (missing(cores)) {
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
          cl = cl, X = 1:ndom, FUN = function(i)  {         
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
                                    error <- data.frame(erro[[i]])
                                    strat <- data.frame(solut[[2]])
                                    solut[[2]]$SOLUZ <- sum(bethel(strat, error, 
                                                                   realAllocation = T))
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
        
        for (i in (1:ndom)) {
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
        for (i in 1:ndom) {
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
            solut <- strataGenalg(errors = erro[[i]], 
                                  strata = stcamp[[i]], cens = cens, strcens = flagcens, 
                                  dominio = i, initialStrata = nstrata[i], 
                                  minnumstr, iter, pops, mut_chance, elitism_rate, 
                                  addStrataFactor, highvalue, suggestions = suggest, 
                                  realAllocation, writeFiles, showPlot)
            if (nrow(stcamp[[i]]) == 1) {
              solut <- list(c(1), stcamp[[i]][c(1:grep("DOM1", 
                                                       colnames(stcamp[[i]])))])
              error <- data.frame(erro[[i]])
              strat <- data.frame(solut[[2]])
              solut[[2]]$SOLUZ <- sum(bethel(strat, error, 
                                             realAllocation = T))
              if (solut[[2]]$SOLUZ > solut[[2]]$N) 
                solut[[2]]$SOLUZ <- solut[[2]]$N
            }
            vettsol <- c(vettsol, solut[[1]])
            if (length(outstrata) > 0) 
              colnames(outstrata) <- toupper(colnames(outstrata))
            colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
            outstrata <- rbind(outstrata, solut[[2]])
            rbga.object <- solut[[3]]
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
