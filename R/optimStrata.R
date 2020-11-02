#------------------------------------------
# optimStrata
#------------------------------------------
# wrapper function to call the different 
# optimization functions:
# 1. optimizeStrata (method = "atomic")
# 2. optimizeStrata (method = "continuous")
# 1. optimizeStrata (method = "spatial")
#------------------------------------------

optimStrata <- function(method=c("atomic","continuous","spatial"),
                        # common parameters
                        framesamp,
                        framecens=NULL,
                        model=NULL,
                        nStrata=NA,
                        errors,
                        alldomains=TRUE,
                        dom=NULL,
                        strcens=FALSE,
                        minnumstr=2,
                        iter=50,
                        pops=20,
                        mut_chance=NA,
                        elitism_rate=0.2,
                        suggestions=NULL,
                        writeFiles=FALSE,
                        showPlot=TRUE,
                        parallel=TRUE,
                        cores=NA,
                        # parameters only for optimizeStrataSpatial
                        fitting=NA,
                        range=NA,
                        kappa=NA) 
{
  # Control of method
  if(!(method %in% c("atomic","continuous","spatial"))) stop("Method should be one in 'atomic','continuous','spatial'")
  # if (missing(cores)) cores <- 1
  # Control of common parameters
  if (alldomains==TRUE & !is.null(dom)) stop("Processing of all domains set TRUE, but a given domain has been indicated")
  # Method 'atomic'
  if (method == "atomic") {
	  if (is.null(errors)) stop("The 'precision constraints' (errors) dataframe is missing")
    if (is.null(framesamp)) stop("The 'sampling frame' (framesamp) dataframe is missing")
    if (!is.null(framesamp)) checkInput(errors, sampframe=framesamp)
    if (!is.null(framecens)) checkInput(errors, sampframe=framecens)
#     checkInput(errors, strata)
# 	  if (!is.null(cens) & is.null(strcens)) stop("Takeall strata presence indicated (strcens=TRUE), but no related strata dataframe (cens) is given")
    nvarX <- length(grep("X",colnames(framesamp)))
    for (i in (1:nvarX)) {
      st <- paste0("if (!is.numeric(framesamp$X",i,")) stop('Stratification variable(s) must be numeric - Transform them to numeric from factor or character before run optimization')"
)
      eval(parse(text=st))
    }
    strata <- buildStrataDF(framesamp,model=model,progress = FALSE)
    if (!is.null(framecens)) {
      cens <- buildStrataDF(framecens,model=model,progress = FALSE)
      strcens <- TRUE
    }
    if (is.null(framecens)) {
      cens <- NULL
      strcens <- FALSE
    }
    if (!is.na(nStrata[1])) {
      initialStrata <- nStrata
      addStrataFactor <- 0.0
    }
    if (is.na(nStrata[1])) {
      initialStrata <- NA
      addStrataFactor <- 0.0
    }
	  solut <- optimizeStrata(
    	errors = errors, 
    	strata = strata, 
    	cens = cens, 
    	strcens = strcens,
    	alldomains = alldomains,
    	dom = dom,	
    	initialStrata = initialStrata, 
    	addStrataFactor = addStrataFactor, 
    	minnumstr = minnumstr, 
    	iter = iter, 
    	pops = pops, 
    	mut_chance = mut_chance, 
    	elitism_rate = elitism_rate,
    	highvalue = 1e+08, 
    	suggestions = suggestions,
    	realAllocation = TRUE,
    	writeFiles = writeFiles,
    	showPlot = showPlot, 
    	parallel = parallel,
    	cores = cores
	    )	
    newstrata <- updateStrata(strata, solut)
    framenew <- updateFrame(frame=framesamp,newstrata=newstrata)
    if (!is.null(framecens)) {
      colnames(framecens) <- toupper(colnames(framecens))
      framenew$STRATUM <- as.character(framenew$STRATUM)
      framecens$LABEL <- 99999
      framecens$STRATUM <- "99999"
      framenew <- rbind(framenew,framecens)
      framenew$STRATUM <- as.numeric(framenew$LABEL)
      nvarX <- length(grep("X",colnames(framecens)))
      for (i in c(1:nvarX)) {
        st <- paste("framecens$X",i," <- NULL",sep="")
        eval(parse(text=st))
      }
      framecens$X1 <- 99999
      cens <- buildStrataDF(framecens,progress=FALSE)
      cens$X1 <- NULL
      cens$SOLUZ <- cens$N
      cens$CENS <- 1
      solut$aggr_strata <- rbind(solut$aggr_strata,cens)
    }
    solution <- list(indices = solut$indices, framenew=framenew,aggr_strata=solut$aggr_strata)
  }
  # Method 'continuous'
  if (method == "continuous") {
	  if (is.null(errors)) stop("The 'precision constraints' (errors) dataframe is missing")
	  # if (!is.null(strata)) stop("Strata dataframe is not required with this method")
	  # if (!is.na(initialStrata)) stop("Initial number of strata is not required with this method")
	  # if (!is.na(addStrataFactor)) stop("'addStrataFactor' is not required with this method")
	  # if (!is.null(cens)) stop("Takeall strata dataframe is not required with this method")
    if (is.null(framesamp)) stop("The 'sampling frame' (framesamp) dataframe is missing")
    checkInput(errors, sampframe=framesamp)
    if (!is.null(framecens)) checkInput(errors, sampframe=framecens) 
    if (!is.null(framecens)) strcens <- TRUE
	  if (!is.na(fitting[1])) stop("Fitting value(s) not required with this method")
	  if (!is.na(range[1])) stop("Range value(s) not required with this method")
	  if (!is.na(kappa)) stop("Kappa value not required with this method")
	  solution <- optimizeStrata2(
      errors = errors, 
      framesamp = framesamp,
      framecens = framecens, 
      strcens = strcens, 
      model = model, 
      alldomains = alldomains, 
      dom = dom, 
      nStrata = nStrata, 
      minnumstr = minnumstr, 
      iter = iter, 
      pops = pops, 
      mut_chance = mut_chance, 
      elitism_rate = elitism_rate, 
      highvalue = 1e+08, 
      suggestions = suggestions, 
      realAllocation = TRUE, 
      writeFiles = writeFiles, 
      showPlot = showPlot, 
      parallel = parallel, 
      cores = cores
	)
  }
  # Method 'spatial'
  if (method=="spatial") {
    checkInput(errors, sampframe=framesamp)
    # if (!is.null(framecens)) checkInput(errors, sampframe=framecens)
	  nvarY <- length(grep("Y", names(framesamp)))
	  if (is.na(fitting[1])) stop("Fitting values of spatial models must be given")
	  if (is.na(range[1])) stop("Range values of spatial models must be given")
	  if (is.na(kappa)) kappa <- 3
	  if (nvarY != length(as.numeric(fitting))) stop("Fitting values must be equal to the number of Y's")
	  if (nvarY != length(as.numeric(range))) stop("Range values must be equal to the number of Y's")
	  nvars <- length(grep("var", names(framesamp)))
	  if (nvarY != nvars) stop("Variances in the 'framesamp' dataframe must be given (one for each Y)")
	  for (i in (1:nvars)) {
	    stmt <- paste("if (min(framesamp$var",i,") < 0) stop('Variance var",i," of variable Y",i," has negative values in framesamp')",sep="")
	  }
	  if (sum(grep("lon",colnames(framesamp))) == 0 | sum(grep("lat",colnames(framesamp))) == 0)  stop("Coordinates (lon and lat) must be given in 'framesamp' dataframe")
	  
	  if (!is.null(framecens)) {
	    checkInput(errors, sampframe=framecens)
	    strcens <- TRUE
	    nvars <- length(grep("var", names(framecens)))
	    if (nvarY != nvars) stop("Variances in the 'framecens' dataframe must be given (one for each Y)")
	    for (i in (1:nvars)) {
	      stmt <- paste("if (min(framecens$var",i,") < 0) stop('Variance var",i," of variable Y",i," has negative values in framecens')",sep="")
	      eval(parse(text=stmt))
	   }
	    if (sum(grep("lon",colnames(framecens))) == 0 | sum(grep("lat",colnames(framecens))) == 0)  stop("Coordinates (lon and lat) must be given in 'framecens' dataframe")
	  }
	  
	  solution <- optimizeStrataSpatial(
      errors = errors, 
      framesamp = framesamp,
      framecens = framecens, 
      strcens = strcens, 
      alldomains = alldomains, 
      dom = dom, 
      nStrata = nStrata, 
      minnumstr = minnumstr, 
      iter = iter, 
      pops = pops, 
      mut_chance = mut_chance, 
      elitism_rate = elitism_rate, 
      highvalue = 1e+08, 
      suggestions = suggestions, 
      realAllocation = TRUE, 
      writeFiles = writeFiles, 
      showPlot = showPlot, 
      parallel = parallel, 
      cores = cores,
      fitting = fitting,
      range = range,
      kappa = kappa
	  )
  }
  return(solution)
}