optimizeStrataSpatial <- 
  function (errors, 
            framesamp,
            framecens = NULL, 
            strcens = FALSE, 
            alldomains = TRUE, 
            dom = NULL, 
            nStrata = c(5), 
            fitting=c(1),
            range=c(0),
            kappa=3,
            minnumstr = 2, 
            iter = 50, 
            pops = 20, 
            mut_chance = NA, 
            # mutationFactor = 0.5,
            elitism_rate = 0.2, 
            highvalue = 1e+08, 
            suggestions = NULL, 
            realAllocation = TRUE, 
            writeFiles = FALSE, 
            showPlot = TRUE, 
            parallel = TRUE, 
            cores = NA) 
  { 
    if (strcens == FALSE) {
      cens=NULL
      censi=NULL
    }
    frame <- framesamp
    if (alldomains == TRUE) dom <- NULL
    colnames(frame) <- toupper(colnames(frame))
    nvarX <- length(grep("X", names(framesamp)))
    ndom <- length(levels(as.factor(frame$DOMAINVALUE)))
    if (alldomains == FALSE) {
      frame <- frame[frame$DOMAINVALUE == dom,]
    }
    if (strcens == TRUE & is.null(framecens))
      stop("No data in the cens dataframe")
    if (strcens == TRUE) {
      censi <- NULL
      if (alldomains == FALSE) {
        colnames(framecens) <- toupper(colnames(framecens))
        framecens <- framecens[framecens$DOMAINVALUE == dom,]
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
            eval(parse(text=paste("framecens$X",i," <- NULL",sep="")))
          }
        }  
        cens <- buildStrataDFSpatial(dataset=framecens,
                                     fitting=fitting,
                                     range=range,
                                     kappa=kappa,
                                     progress=FALSE,verbose=FALSE)
        cens$CENS <- 1
        censtot <- cens
      }
    }
    if (writeFiles == TRUE) {
      dire <- getwd()
      direnew <- paste(dire, "/output", sep = "")
      if (dir.exists(direnew)) 
        unlink(direnew,recursive=TRUE)
      if (!dir.exists(direnew)) 
        dir.create(direnew)
      #setwd(direnew)
    }
    if (parallel == TRUE & ndom == 1) parallel <- FALSE
    if(parallel == FALSE & !missing(cores) | !is.na(cores)){
      cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
      cores <- 1
      Sys.sleep(0.314)
    }
    if(alldomains == FALSE & (parallel == TRUE | !missing(cores) | !is.na(cores))){
      cat("Sequential optimization as parallel = FALSE, defaulting number of cores = 1")
      cores <- 1
      Sys.sleep(0.314)
    }
    # ncuts <- nStrata - 1
    # if (is.na(initialStrata)) 
    #   initialStrata <- as.numeric(table(strata$DOM1))
    # nstrata = initialStrata
    colnames(errors) <- toupper(colnames(errors))
    colnames(frame) <- toupper(colnames(frame))
    errors$DOMAINVALUE <- as.factor(errors$DOMAINVALUE)
    erro <- split(errors, list(errors$DOMAINVALUE))
    # stcamp <- split(strata, list(strata$DOM1))
    stcamp <- split(frame, list(frame$DOMAINVALUE))
    if (!is.null(suggestions)) {
      # nStrataSuggested <- nrow(sugg[sugg$domainvalue==1,])
      # if (nStrataSuggested != nStrata) stop("Number of strata in 'suggestions' is different from 'nStrata' value")
      for (i in (1:ndom)) {
        nvalues <- nrow(suggestions[suggestions$domainvalue==i,])
        if (nvalues != nvarX*(nStrata[i]-1)) stop("Number of values in suggestions not compatible with nStrata")
      }
      suggestdom <- split(suggestions, list(suggestions$domainvalue))
    }
    if (strcens == TRUE & !is.null(cens) > 0) {
      colnames(cens) <- toupper(colnames(cens))
      # k <- length(levels(as.factor(strata$DOM1)))
      k <- length(levels(as.factor(frame$DOMAINVALUE)))
      stcens <- NULL
      for (i in (1:k)) {
        stcens[[i]] <- cens[cens$DOM1 == i, ]
      }
    }
    # ndom <- length(levels(as.factor(strata$DOM1)))
    if (alldomains == TRUE) {
      if (ndom > length(nStrata))
        stop("'initialStrata' vector lenght (=", length(nStrata),
             ") \nis not compatible with number of domains (=",
             ndom, ")\nSet initialStrata with a number of elements equal to the number of domains")
      vettsol <- NULL
      outstrata <- NULL
      if (parallel) {
        if (missing(cores) | is.na(cores)) {
          cores <- parallel::detectCores() - 1
          if (ndom < cores) 
            cores <- ndom
        }
        if (cores < 2) 
          stop("\nOnly one core available: no parallel processing possible. 
               \nPlease change parameter parallel to FALSE and run again")
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
            # cens <- NULL
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
              suggest <- matrix(0, nrow = 2, ncol = (nStrata[i]-1)*nvarX)
              suggest[1, ] <- suggestdom[[i]]$suggestions1
              suggest[2, ] <- suggestdom[[i]]$suggestions2
            }
            if (!is.null(suggestions) & alldomains == FALSE) {
              suggest <- matrix(0, nrow = 2, ncol = (nStrata[i]-1)*nvarX)
              suggest[1, ] <- suggestdom[[i]]$suggestions1
              suggest[2, ] <- suggestdom[[i]]$suggestions2
            }
            if (is.null(suggestions)) {
              suggest <- NULL
            }
            if (nrow(stcamp[[i]]) > 0) {
              solut <- strataGenalgSpatial(errors = erro[[i]], 
                                            frame = stcamp[[i]],
                                            cens = cens, 
                                            strcens = flagcens, 
                                            fitting  = fitting,
                                            range = range,
                                            kappa = kappa,
                                            ncuts = (nStrata[i] - 1),
                                            dominio = i, 
                                            minnumstr, 
                                            iter, 
                                            pops, 
                                            mut_chance, 
                                            elitism_rate, 
                                            suggestions = suggest, 
                                            realAllocation, 
                                            writeFiles, 
                                            showPlot)
              if (nrow(stcamp[[i]]) == 1) {
                solut <- list(c(1), 
                              stcamp[[i]][c(1:grep("DOMAINVALUE", 
                                                   colnames(stcamp[[i]])))])
                error <- data.frame(erro[[i]],stringsAsFactors = TRUE)
                strat <- data.frame(solut[[2]],stringsAsFactors = TRUE)
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
        #       vettsol <- do.call(c, lapply(par_ga_sol, `[[`, 1)) 
        vettsol <- do.call(rbind, lapply(par_ga_sol, `[[`, 1))
        outstrata <- do.call(rbind, lapply(par_ga_sol, `[[`, 2))
        results <- do.call(rbind, lapply(par_ga_sol, `[[`, 3))
        
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
          # cens <- NULL
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
            suggest <- matrix(0, nrow = 2, ncol = (nStrata[i]-1)*nvarX)
            suggest[1, ] <- suggestdom[[i]]$suggestions1
            suggest[2, ] <- suggestdom[[i]]$suggestions2
          }
          if (!is.null(suggestions) & alldomains == FALSE) {
            suggest <- matrix(0, nrow = 2, ncol = (nStrata[i]-1)*nvarX)
            suggest[1, ] <- suggestdom[[i]]$suggestions1
            suggest[2, ] <- suggestdom[[i]]$suggestions2
          }
          if (is.null(suggestions)) {
            suggest <- NULL
          }
          if (nrow(stcamp[[i]]) > 0) {
            solut <- strataGenalgSpatial(errors = erro[[i]], 
                                          frame = stcamp[[i]],
                                          cens = cens, 
                                          strcens = flagcens, 
                                          fitting = fitting,
                                          range = range,
                                          kappa = kappa,
                                          ncuts = (nStrata[i] - 1),
                                          dominio = i, 
                                          minnumstr, 
                                          iter, 
                                          pops, 
                                          mut_chance, 
                                          elitism_rate, 
                                          suggestions = suggest, 
                                          realAllocation, 
                                          writeFiles, 
                                          showPlot)
            if (nrow(stcamp[[i]]) == 1) {
              solut <- list(c(1), stcamp[[i]][c(1:grep("DOMAINVALUE", 
                                                       colnames(stcamp[[i]])))])
              error <- data.frame(erro[[i]],stringsAsFactors = TRUE)
              strat <- data.frame(solut[[2]],stringsAsFactors = TRUE)
              solut[[2]]$SOLUZ <- sum(bethel(strat, error, 
                                             realAllocation = T))
              if (solut[[2]]$SOLUZ > solut[[2]]$N) 
                solut[[2]]$SOLUZ <- solut[[2]]$N
            }
            # vettsol <- c(vettsol, solut[[1]])
            vettsol <- rbind(vettsol, solut[[1]])
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
      # if (dom < 1 | dom > ndom) 
      #   stop("\nInvalid value of the indicated domain\n")
      i <- dom
      erro[[i]] <- erro[[i]][, -ncol(errors)]
      flagcens <- strcens
      if (strcens == TRUE) {
        flagcens <- TRUE
        if (!is.null(cens)) colnames(cens) <- toupper(colnames(cens))
        censi <- cens[cens$DOM1 == i, ]
        if (nrow(censi) == 0) {
          flagcens <- FALSE
        }
      }
      if (!is.null(suggestions) & alldomains == TRUE) {
        suggest <- matrix(0, nrow = 1, ncol = (nStrata[i]-1)*nvarX)
        suggest[1, ] <- suggestdom[[i]]$suggestions
      }
      if (!is.null(suggestions) & alldomains == FALSE) {
        suggest <- matrix(0, nrow = 1, ncol = (nStrata[i]-1)*nvarX)
        suggest[1, ] <- suggestions$suggestions
      }
      if (is.null(suggestions)) {
        suggest <- NULL
      }
      solut <- strataGenalgSpatial(errors = erro[[i]], 
                                    frame = frame,
                                    # frame = stcamp[[i]], 
                                    cens = cens, 
                                    strcens = flagcens, 
                                    fitting = fitting,
                                    range = range,
                                    kappa = kappa,
                                    ncuts = (nStrata[i] - 1),
                                    dominio = i, 
                                    minnumstr, 
                                    iter, 
                                    pops, 
                                    mut_chance, 
                                    elitism_rate, 
                                    suggestions = suggest, 
                                    realAllocation, 
                                    writeFiles, 
                                    showPlot)
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
    vettsoldf <- as.data.frame(vettsol,stringsAsFactors = TRUE)
    colnames(vettsoldf) <- c("ID","LABEL")
    vettsoldf$STRATO <- vettsoldf$LABEL
    framenew <- merge(frame,vettsoldf,by=c("ID"))
    # if (strcens == TRUE & !is.null(censi)) {
    if (strcens == TRUE) {
      if (alldomains == FALSE) {
        # colnames(framecens) <- toupper(colnames(framecens))
        # colnames(framecensold) <- toupper(colnames(framecensold))
        framecens <- framecens[framecens$DOMAINVALUE == dom,]
        framecensold <- framecensold[framecensold$DOMAINVALUE == dom,]
      }
      for (i in (1:nvarX)) {
        eval(parse(text=paste("framecens$X",i," <- framecensold$X",i,sep="")))
      }
      framecens$STRATO <- max(nStrata) + 1
      framecens$LABEL <- max(nStrata) + 1
      colnames(framecens) <- toupper(colnames(framecens))
      framenew <- rbind(framenew,framecens)
      censtot$SOLUZ <- censtot$N
      outstrata <- rbind(outstrata,censtot)
    }
    #-----------------------------------------    
    # new to tackle with erroneous allocation     
    # dataset <- framenew 
    # nX <- sum(grepl("X",colnames(frame)))
    # for(j in 1:nX){
    #   eval(parse(text=paste("frame$X",j," <- NULL",sep="")))
    # }
    # dataset$X1 <- framenew$LABEL
    # outstrata2 <- buildStrataDFSpatial(dataset,progress=FALSE,verbose=FALSE)
    # outstrata2$SOLUZ <- outstrata$SOLUZ
    #-----------------------------------------    
    solution <- list(indices = vettsol, aggr_strata = outstrata, framenew = framenew)
    if (writeFiles == TRUE) {
      setwd(dire)
    }
    return(solution)
  }
