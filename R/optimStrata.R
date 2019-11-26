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
                        errors,
                        alldomains=TRUE,
                        dom=NULL,
                        minnumstr=2,
                        iter=50,
                        pops=20,
                        mut_chance=NA,
                        elitism_rate=0.2,
                        suggestions=NULL,
                        writeFiles=FALSE,
                        showPlot=TRUE,
                        parallel=TRUE,
                        cores,
                        # parameters only for optimizeStrata
                        strata=NULL,
                        initialStrata=NA,
                        addStrataFactor=0.0,
						            cens=NULL,
                        strcens=NULL,
                        # parameters only for optimizeStrata2 and optimizeStrataSpatial
                        framesamp=NULL,
                        framecens=NULL,
                        model=NULL,
                        nStrata=5,
                        # parameters only for optimizeStrataSpatial
                        fitting=1,
                        range=c(1),
                        kappa=3
                        ) 
{
  # Control of method
  if(!(method %in% c("atomic","continuous","spatial"))) stop("Method should be one in 'atomic','continuous','spatial'")
  # Control of common parameters
  if (alldomains==TRUE & is.null(dom)) stop("Processing of all domains set TRUE, but a given domain has been indicated")
  # Method 'atomic'
  if (method == "atomic") {
	  if (is.null(errors)) stop("The 'precision constraints' (errors) dataframe is missing")
    checkInput(errors, strata)
	  if (!is.null(cens) & is.null(strcens)) stop("Takeall strata presence indicated (strcens=TRUE), but no related strata dataframe (cens) is given")
	  optimizeStrata(
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
  }
  # Method 'continuous'
  if (method == "continuous") {
	  if (is.null(errors)) stop("The 'precision constraints' (errors) dataframe is missing")
	  if (!is.null(strata)) stop("Strata dataframe is not required with this method")
	  if (!is.na(initialStrata)) stop("Initial number of strata is not required with this method")
	  if (!is.null(addStrataFactor)) stop("'addStrataFactor' is not required with this method")
	  if (!is.null(cens)) stop("Takeall strata dataframe is not required with this method")
    checkInput(errors, frame=framesamp)
	  if (!is.null(fitting)) stop("Fitting value(s) not required with this method")
	  if (!is.null(range)) stop("Range value(s) not required with this method")
	  if (!is.null(range)) stop("Kappa value not required with this method")
	  if (!is.null(framecens)) checkInput(errors, frame=framecens) 
	  optimizeStrata2(
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
	  nvarY <- length(grep("Y", names(framesamp)))
	  if (is.null(fitting)) stop("Fitting values of spatial models must be given")
	  if (is.null(range)) stop("Range values of spatial models must be given")
	  if (nvarY != length(as.numeric(fitting))) stop("Fitting values must be equal to the number of Y's")
	  if (nvarY != length(as.numeric(range))) stop("Range values must be equal to the number of Y's")
	  nvars <- length(grep("var", names(framesamp)))
	  if (nvarY != nvars) stop("Variances in the frame dataframe must be given (one for each Y)")
	  if (grepl("lon",colnames(framesamp)) == FALSE | grepl("lat",colnames(framesamp)))  stop("Coordinates (lon and lat) must be given in the frame dataframe")
	  optimizeStrataSpatial(
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
      cores = cores,
      fitting = fitting,
      range = range,
      kappa = kappa
	  )
  }
}