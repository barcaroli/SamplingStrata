# It optimizes a binary chromosome using a genetic algorithm.
#
#----- Modification
# dataset          = sampling frame data
# errors           = precision constraints
#----- End modification
# string           = string to optimize
# popSize          = the population size
# iters            = number of generations
# mutationChance   = chance that a var in the string gets mutated
rbga2 <- function(dataset,
                  cens,
                  strcens,
                  model,
                  minnumstr,
                 errors,
                 ncuts,
                 stringMin=c(), 
                 stringMax=c(),
                 suggestions=NULL,
                 popSize=200, 
                 iters=100, 
                 mutationChance=NA,
                 mutationFactor=0.5,
                 elitism=NA,
                 monitorFunc=NULL, 
                 evalFunc=NULL,
                 showSettings=FALSE, 
                 verbose=FALSE,
                 # if TRUE, the elite is selected by integer sample size
                 # (attribute "ntotInt" returned by evalFunc), ties broken
                 # by the continuous one
                 eliteByInteger=getOption("SamplingStrata.eliteByInteger", TRUE)) {

    if (is.null(evalFunc)) {
        stop("A evaluation function must be provided. See the evalFunc parameter.");
    }
    
    vars = length(stringMin);
    if (is.na(mutationChance)) {
        mutationChance = 1/(vars+1);
    }
    if (is.na(elitism)) {
        elitism = floor(popSize/5)
    }
    
    # TODO: should do a variaty of sanity checks first
    if (verbose) cat("Testing the sanity of parameters...\n");
    if (length(stringMin) != length(stringMax)) {
        stop("The vectors stringMin and stringMax must be of equal length.");
    }
    if (popSize < 5) {
        stop("The population size must be at least 5.");
    }
    if (iters < 1) {
        stop("The number of iterations must be at least 1.");
    }
    if (!(elitism < popSize)) {
        stop("The population size must be greater than the elitism.");
    }
    
    if (showSettings) {
        if (verbose) cat("The start conditions:\n");
        result = list(stringMin=stringMin, stringMax=stringMax, suggestions=suggestions,
                      popSize=popSize, iters=iters,
                      elitism=elitism, mutationChance=mutationChance);
        class(result) = "rbga";
        
        cat(summary(result));
    } else {
        if (verbose) cat("Not showing GA settings...\n");
    }
    
    if (vars > 0) {
        if (!is.null(suggestions)) {
            if (verbose) cat("Adding suggestions to first population...\n");
            population = matrix(nrow=popSize, ncol=vars);
            # suggestionCount = dim(suggestions)[1]
            # for (i in 1:suggestionCount) {
            #     population[i,] = suggestions[i,]
            # }
            suggestionCount = 2
            population[1,] = suggestions[1,]
            population[2,] = suggestions[2,]
            if (verbose) cat("Filling others with random values in the given domains...\n");
            for (var in 1:vars) {
                population[(suggestionCount+1):popSize,var] = stringMin[var] +
                                   runif(popSize-suggestionCount)*(stringMax[var]-stringMin[var]);
            }
        } else {
            if (verbose) cat("Starting with random values in the given domains...\n");
            # start with an random population
            population = matrix(nrow=popSize, ncol=vars);
            # fill values
            for (var in 1:vars) {
                population[,var] = stringMin[var] +
                                   runif(popSize)*(stringMax[var]-stringMin[var]);
            }
        }
        
        # do iterations
        bestEvals = rep(NA, iters);
        meanEvals = rep(NA, iters);
        bestSolution = rep(NA, iters);
        bestStrata = rep(NA, iters);
        evalVals = rep(NA, popSize);
        evalInt = rep(NA, popSize);
        for (iter in 1:iters) {
            if (verbose) cat(paste("Starting iteration", iter, "\n"));

            # calculate each object
            if (verbose) cat("Calculating evaluation values... ");
            for (object in 1:popSize) {
                if (is.na(evalVals[object])) {
                         #---- Modification:
                      res <- evalFunc(dataset,
                                      cens,
                                      strcens,
                                      model,
                                      minnumstr,
                                      errors,
                                      population[object,],
                                      ncuts);
                      evalVals[object] = sum(res)
                      evalInt[object] = if (is.null(attr(res, "ntotInt"))) sum(res) else attr(res, "ntotInt")
                    #---- End modification:
                    if (verbose) cat(".");
                }
            }
            bestEvals[iter] = min(evalVals);
            meanEvals[iter] = mean(evalVals);
            # cat("\nBest solution: ",bestEvals[iter])
            if (verbose) cat(" done.\n");
            
            if (!is.null(monitorFunc)) {
                if (verbose) cat("Sending current state to rgba.monitor()...\n");
                # report on GA settings
                result = list(type="floats chromosome",
                              stringMin=stringMin, stringMax=stringMax,
                              popSize=popSize, iter=iter, iters=iters,
                              population=population, elitism=elitism, mutationChance=mutationChance,
                              evaluations=evalVals, best=bestEvals, mean=meanEvals);
                class(result) = "rbga";
                
                monitorFunc(result);
            }
            
            if (iter < iters) { # ok, must create the next generation
                if (verbose) cat("Creating next generation...\n");
              # cat("\n Iteration  ",iter)
                newPopulation = matrix(nrow=popSize, ncol=vars);
                newEvalVals = rep(NA, popSize);
                newEvalInt = rep(NA, popSize);
                
                if (verbose) cat("  sorting results...\n");
                if (eliteByInteger) {
                    # rank by integer sample size, then by continuous one
                    ord = order(evalInt, evalVals);
                    sortedEvaluations = list(x=evalVals[ord], ix=ord);
                } else {
                    sortedEvaluations = sort(evalVals, index=TRUE);
                }
                sortedInt = evalInt[sortedEvaluations$ix];
                sortedPopulation  = matrix(population[sortedEvaluations$ix,], ncol=vars);
                
                # save the best
                if (elitism > 0) {
                    if (verbose) cat("  applying elitism...\n");
                    newPopulation[1:elitism,] = sortedPopulation[1:elitism,];
                    newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
                    newEvalInt[1:elitism] = sortedInt[1:elitism]
                } # ok, save nothing
                
                # fill the rest by doing crossover
                if (vars > 1) {
                    if (verbose) cat("  applying crossover...\n");
                    for (child in (elitism+1):popSize) {
                        # ok, pick two random parents
                        parentProb = dnorm(1:popSize, mean = 0, sd = (popSize/3))
    			              parentIDs = sample(1:popSize, 2, prob = parentProb)
                        #parentIDs = sample(1:popSize, 2)
                        parents = sortedPopulation[parentIDs,];
                        crossOverPoint = sample(0:vars,1);
                        if (crossOverPoint == 0) {
                            newPopulation[child, ] = parents[2,]
                            newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
                            newEvalInt[child] = sortedInt[parentIDs[2]]
                        } else if (crossOverPoint == vars) {
                            newPopulation[child, ] = parents[1,]
                            newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
                            newEvalInt[child] = sortedInt[parentIDs[1]]
                        } else {
                            newPopulation[child, ] =
                                c(parents[1,][1:crossOverPoint],
                                  parents[2,][(crossOverPoint+1):vars])
                        }
                    }
    } else { # otherwise nothing to crossover
                    if (verbose) cat("  cannot crossover (#vars=1), using new randoms...\n");
                    # random fill the rest
                    newPopulation[(elitism+1):popSize,] = 
                        sortedPopulation[sample(1:popSize, popSize-elitism),];
                }
                
                population = newPopulation;
                evalVals   = newEvalVals;
                evalInt    = newEvalInt;
                
                # do mutation
                if (mutationChance > 0) {
                    if (verbose) cat("  applying mutations... ");
                    mutationCount = 0;
                    for (object in (elitism+1):popSize) { # don't mutate the best
                        for (var in 1:vars) {
                            if (runif(1) < mutationChance) { # ok, do mutation
                                # OPTION 1
                                # mutate to something random
                                # mutation = stringMin[var] +
                                #    runif(1)*(stringMax[var]-stringMin[var]);
                                # 
                                # OPTION 2
                                # mutate around solution
                              #-------- Modification :
                              if (mutationFactor == 0) {
                                dempeningFactor = (iters-iter)/iters
                                direction       = sample(c(-1,1),1)
                                mutationVal     = (stringMax[var]-stringMin[var])*0.01
                                mutation = population[object,var] + direction*mutationVal
                                # * dempeningFactor
                              }
                              if (mutationFactor > 0) {
                                direction = sample(c(-1,1),1)
                                mutation = population[object,var] * (1+(mutationFactor*direction))
                                # mutation = population[object,var] * (1-mutationFactor)

                              }
                              # but in domain. if not, then take random
                              if (mutation < stringMin[var]) 
                                mutation = stringMin[var] +
                                  runif(1)*(stringMax[var]-stringMin[var]);
                              if (mutation > stringMax[var]) 
                                mutation = stringMin[var] +
                                  runif(1)*(stringMax[var]-stringMin[var]);                                
                                # apply mutation, and delete known evalutation value
                                population[object,var] = mutation;
                                evalVals[object] = NA;
                                evalInt[object] = NA;
                                mutationCount = mutationCount + 1;
                            }
                        }
                    }
                    if (verbose) cat(paste(mutationCount, "mutations applied\n"));
                }
            }
        }
    }

    # report on GA settings
    result = list(type="floats chromosome", 
                  stringMin=stringMin, stringMax=stringMax,
                  popSize=popSize, iters=iters, suggestions=suggestions,
                  population=population, elitism=elitism, mutationChance=mutationChance,
                  evaluations=evalVals, best=bestEvals, mean=meanEvals);
    class(result) = "rbga";

    return(result);
}

