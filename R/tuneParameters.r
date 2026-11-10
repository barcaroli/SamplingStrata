# ==============================================================================
# File:      tuneParameters.r
# Functions: tuneParameters() [exported]
#
# Purpose:
#   Runs a series of noptim optimisations with optimizeStrata() on one domain,
#   with the parameters of the genetic algorithm varying as given in vectors;
#   for each solution it selects nsampl samples from the frame and compares
#   the solutions (number of strata, cost, simulated CVs) to help choosing the
#   most suitable parameter values.
# ==============================================================================

# ----------------------------------------------------
# Function to choose suitable values of the parameters
# for the function "optimizeStrata"
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
#' @title Execution and compared evaluation of optimization runs
#'
#' @description
#' This function allows to execute a number of optimization runs,
#' varying in a controlled way the values of the parameters, in order to find their most suitable values.
#' by comparing the resulting solutions.
#' It can be applied only to a given domain per time.
#' Most parameters of this function are the same than those of the function 'optimizeStrata', but they are given
#' in a vectorial format. The length of each vector is given by the number of optimizations to be run:
#' it is therefore possible to define different combination of values of the parameters for each execution of 'optimizeStrata'.
#' After each optimization run, from the corrisponding optimized frame a given number of samples are drawn. For each of them, the
#' estimates of the target variables Y's are computed ("precision"), together with the associated coefficients of variations,
#' and the absolute differences between the values of the estimates and the true values in the population ("bias").
#' Information on the distribution of bias (differences) and precision (CV's) are outputted, and also boxplots
#' for each of them are produced, in order to permit a compared evaluation of the different solutions found in the different runs.
#' As the optimal solution is stored for each run, after the evaluation it is possible to use it directly, or as a
#' "suggestion" for a new optimization with more iterations (in order to improve it).
#'
#' @param noptim Number of optimization runs to be performed
#' @param nsampl Number of samples to be drawn from the optimized population frame after each optimization
#' @param frame The (mandatory) dataframe containing the sampling frame
#' @param errors This is the (mandatory) dataframe containing the precision levels expressed in terms of Coefficients of Variation
#'   that estimates on target variables Y's of the survey must comply
#' @param strata This is the (mandatory) dataframe containing the information related to "atomic" strata, i.e. the strata obtained by
#'   the Cartesian product of all auxiliary variables X's. Information concerns the identifiability of strata
#'   (values of X's) and variability of Y's (for each Y, mean and standard error in strata)
#' @param cens This the (optional) dataframe containing the takeall strata, those strata whose units must be selected in
#'   whatever sample. It has same structure than "strata" dataframe
#' @param strcens Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE
#' @param alldomains Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains. It must be left to its default (FALSE)
#' @param dom Indicates the domain on which the optimization runs must be performed. It is an integer value that has to be internal to the interval
#'   (1 <--> number of domains). It is mandatory, if not indicated, the default (1) is taken.
#' @param initialStrata This is the initial limit on the number of strata for each solution. Default is 3000.
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param addStrataFactor This parameter indicates the probability that at each mutation the number of strata may increase with
#'   respect to the current value. Default is 0.01 (1\%).
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param minnumstr Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param iter Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 20.
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param pops The dimension of each generations in terms of individuals. Default is 50.
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param mut_chance Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
#'   of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
#'   but a slower convergence, while low values permit a faster convergence, but the final solution can be
#'   distant from the optimal one. Default is 0.05.
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param elitism_rate This  parameter indicates the rate of better solutions that must be preserved from one generation
#'   to another. Default is 0.2 (20\%).
#'   This parameter has to be given in a vectorial format, whose length is given by the number of different optimisations
#'   ( = value of parameter 'noptim')
#' @param writeFiles Indicates if the various dataframes and plots produced during the execution have to be written in the working
#'   directory.
#'   Default is FALSE.
#'
#' @return
#' A dataframe containing for each iteration the number of strata, the cost of the solution and the values of the expected CV's
#'
#' @author
#' Giulio Barcaroli
#'
#' @examples
#' #
#' \dontrun{
#' #------------------------------------------------------------
#' # data setting
#' library(SamplingStrata)
#' data(swissstrata)
#' data(swisserrors)
#' data(swissframe)
#' # As this function can be applied only to a given domain per time,
#' # we select the first domain
#' frame <- swissframe[swissframe$domainvalue == 1,]
#' strata <- swissstrata[swissstrata$DOM1 == 1,]
#' errors <- swisserrors[swisserrors$domainvalue == 1,]
#' #------------------------------------------------------------
#' # parameters setting
#' noptim <- 10 # Number of runs
#' nsampl <- 100 # Number of samples to be drawn after each optimization
#' initialStrata <- ceiling(c(1:noptim)*0.1*(nrow(strata))) # Number of initial strata
#' addStrataFactor <- rep(0.01,noptim) # Rate for increasing initial strata
#' minnumstr <- rep(2,noptim) # Minimum number of units per stratum
#' iter <- rep(200,noptim) # Number of iterations for each optimization
#' pops <- rep(20,noptim) # Number of solutions for each iteration
#' mut_chance <- rep(0.004,noptim) # Mutation chance
#' elitism_rate <- rep(0.2,noptim) # Elitism rate
#' #------------------------------------------------------------
#' results <- tuneParameters (
#'   noptim,
#'   nsampl,
#'   frame,
#'   errors = errors,
#'   strata = strata,
#'   cens = NULL,
#'   strcens = FALSE,
#'   alldomains = FALSE,
#'   dom = 1,
#'   initialStrata,
#'   addStrataFactor,
#'   minnumstr,
#'   iter,
#'   pops,
#'   mut_chance,
#'   elitism_rate
#'   )
#' results
#' }
#'
#' @keywords survey
#' @export
tuneParameters <- function (noptim, nsampl, frame, errors, strata, 
    cens = NULL, strcens = FALSE, alldomains = FALSE, dom = 1, 
    initialStrata, addStrataFactor, minnumstr, iter, pops, mut_chance, 
    elitism_rate,writeFiles=FALSE) 
{
  if (writeFiles == TRUE) {
    dire <- getwd()
    direnew <- paste(dire,"/simulation",sep="")
    if(!dir.exists(direnew)) dir.create(direnew)
    setwd(direnew)
    on.exit(setwd(dire), add = TRUE)   # restore the working directory also on error
  }
  colnames(frame) <- toupper(colnames(frame))
  frame <- frame[frame$DOMAINVALUE == dom, ]
	checkInput(errors,strata,frame)
    if (length(initialStrata) != noptim) 
        stop("Number of optimisations does not equal values in initialStrata parameter")
    if (length(addStrataFactor) != noptim) 
        stop("Number of optimisations does not equal values in addStrataFactor parameter")
    if (length(minnumstr) != noptim) 
        stop("Number of optimisations does not equal values in minnumstr parameter")
    if (length(iter) != noptim) 
        stop("Number of optimisations does not equal values in iter parameter")
    if (length(pops) != noptim) 
        stop("Number of optimisations does not equal values in pops parameter")
    if (length(mut_chance) != noptim) 
        stop("Number of optimisations does not equal values in mut_chance parameter")
    if (length(elitism_rate) != noptim) 
        stop("Number of optimisations does not equal values in elitism_rate parameter")
    numY <- length(grep("Y", toupper(colnames(frame))))
    for (i in 1:numY) {
        stmt <- paste("Y", i, " <- sum(frame$Y", i, ")", sep = "")
        eval(parse(text = stmt))
    }
    estim <- array(0, c(noptim, nsampl, numY))
    nsimul <- c(1:noptim)
    nstrata <- rep(0, noptim)
    cost <- rep(0, noptim)
    for (i in 1:numY) {
        stmt <- paste("CV", i, " <- rep(0,noptim)", sep = "")
        eval(parse(text = stmt))
    }
    est <- NULL
    for (i in 1:numY) {
        if (i < numY) 
            stmt <- paste("est <- paste(est,'CV", i, "=CV", i, 
                ",',sep='')", sep = "")
        if (i == numY) 
            stmt <- paste("est <- paste(est,'CV", i, "=CV", i, 
                "',sep='')", sep = "")
        eval(parse(text = stmt))
    }
    stmt <- paste("simula <- list(nsimul=nsimul,nstrati=nstrata,cost=cost,", 
        est, ")", sep = "")
    eval(parse(text = stmt))
    simula <- as.data.frame(simula,stringsAsFactors = TRUE)
    nopt <- rep(0, (noptim * nsampl))
    nsamp <- rep(0, (noptim * nsampl))
    for (i in 1:numY) {
        stmt <- paste("diff", i, " <- rep(0,noptim*nsampl)", 
            sep = "")
        eval(parse(text = stmt))
    }
    est <- NULL
    for (i in 1:numY) {
        if (i < numY) 
            stmt <- paste("est <- paste(est,'diff", i, "=diff", 
                i, ",',sep='')", sep = "")
        if (i == numY) 
            stmt <- paste("est <- paste(est,'diff", i, "=diff", 
                i, "',sep='')", sep = "")
        eval(parse(text = stmt))
    }
    stmt <- paste("differ <- list(nopt=nopt,nsamp=nsamp,", est, 
        ")", sep = "")
    eval(parse(text = stmt))
    differ <- as.data.frame(differ,stringsAsFactors = TRUE)
    for (i in (1:noptim)) {
        solution <- optimizeStrata(errors = errors, strata = strata, 
            cens = NULL, strcens = FALSE, alldomains = FALSE, 
            dom = dom, initialStrata = initialStrata[i], addStrataFactor = addStrataFactor[i], 
            minnumstr = minnumstr[i], iter = iter[i], pops = pops[i], 
            mut_chance = mut_chance[i], elitism_rate = elitism_rate[i], 
            realAllocation = TRUE, highvalue = 1e+08, suggestions = NULL, writeFiles=FALSE)
        # (the files outstrata*.txt, solution*.txt, results*.txt and plot*.pdf
        # were removed/copied here, also in the working directory and also with
        # writeFiles = FALSE; optimizeStrata() is called with writeFiles = FALSE
        # and writes none of them, so these operations have been removed)
		outstrata <- solution$aggr_strata
        newstrata <- updateStrata(strata, solution)
        framenew <- updateFrame(frame, newstrata)
        simula$nstrati[i] <- nrow(outstrata)
        simula$cost[i] <- sum(ceiling(outstrata$SOLUZ))
        for (j in (1:nsampl)) {
            samp <- selectSample(framenew, outstrata, verbatim=FALSE)
            for (k in 1:numY) {
                stmt <- paste("estim[i,j,k] <- sum(samp$Y", k, 
                  " * samp$WEIGHTS)", sep = "")
                eval(parse(text = stmt))
                differ$nopt[(i - 1) * nsampl + j] <- i
                differ$nsamp[(i - 1) * nsampl + j] <- j
                stmt <- paste("differ$diff", k, "[(i-1)*nsampl+j]<- estim[i,j,k] - Y", 
                  k, " ", sep = "")
                eval(parse(text = stmt))
            }
        }
        for (k in 1:numY) {
            stmt <- paste("simula$CV", k, "[i] <- sd(estim[i,,k]) / mean(estim[i,,k])", 
                sep = "")
            eval(parse(text = stmt))
        }
        if (writeFiles == TRUE) {
          sink("Iteration.txt", append = TRUE)
          cat("\nIteration :", i, "\n")
          print(simula[i, ])
          sink()
        }
    }
    ind <- which(simula$cost == min(simula$cost))
    cat("\n--------------------")
    cat("\nBest optimization: run # ", ind)
    cat("\nRequired sample cost: ", simula$cost[ind])
    cat("\n--------------------")
    cat("\nValues of parameters")
    cat("\n--------------------")
    cat("\ninitialStrata: ", initialStrata[ind])
    cat("\naddStrataFactor: ", addStrataFactor[ind])
    cat("\nminnumstr: ", minnumstr[ind])
    cat("\niter: ", iter[ind])
    cat("\npops: ", pops[ind])
    cat("\nmut_chance: ", mut_chance[ind])
    cat("\nelitism_rate: ", elitism_rate[ind])
    cat("\n--------------------")
    res <- paste("results_", dom, ".csv", sep = "")
    if (writeFiles == TRUE) {
    write.table(simula, res, row.names = FALSE, col.names = TRUE, 
        quote = FALSE, sep = ",")
    dif <- paste("details_", dom, ".csv", sep = "")
    write.table(differ, dif, row.names = FALSE, col.names = TRUE, 
        quote = FALSE, sep = ",")
		
    stmt <- paste("pdf('strata_costs_", dom, ".pdf', width=14, height=10)", 
        sep = "")
    eval(parse(text = stmt))
	  split.screen(c(1, 2))
    screen(1)
    plot(simula$nsimul, simula$nstrati, xlim = c(1, nrow(simula)), 
        xlab = "Runs", ylab = "Number of strata", 
        bty = "n")
    title("Number of optimal strata / optimization run",cex.main=1.0,font.main=1)
    lines(simula$nsimul, simula$nstrati)
    screen(2)
    plot(simula$nsimul, simula$cost, xlim = c(1, nrow(simula)), xlab = "Runs", 
        ylab = "Sample cost", bty = "n")
    title("Solution cost / optimization run",cex.main=1.0,font.main=1)
    lines(simula$nsimul, simula$cost)
    close.screen(all.screens = TRUE)
    dev.off()
	
    stmt <- paste("pdf('estimates_", dom, ".pdf', width=14, height=10)", 
        sep = "")
    eval(parse(text = stmt))
	nopt <- differ$nopt
    k <- ceiling(numY/4)
    for (j in 1:k) {
        split.screen(c(2, 2))
        for (i in 1:4) {
            if (i + 4 * (j - 1) <= numY) {
                stmt <- paste("screen(", i, ")", sep = "")
                eval(parse(text = stmt))
                stmt <- paste("boxplot(diff", i + 4 * (j - 1), 
                  "~nopt,data=differ,xlab = 'Runs',col='orange')", 
                  sep = "")
                eval(parse(text = stmt))
                stmt <- paste("mtext(expression(Y", i + 
                  4 * (j - 1), "), side=3, adj=0, cex=1.0, line=1)", 
                  sep = "")
                eval(parse(text = stmt))
				title("Distribution of differences between true value and estimate",
					cex.main= 1.0,font.main=1)
            }
        }
        close.screen(all.screens = TRUE)
		dev.off()
    }
    }
    return(simula)
}
