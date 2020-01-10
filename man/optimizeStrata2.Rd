\name{optimizeStrata2}
\Rdversion{1.3-1}
\alias{optimizeStrata2}
\title{
Best stratification of a sampling frame for multipurpose surveys (only with continuous stratification variables)
}
\description{
This function runs a set of other functions to optimise the stratification of a sampling frame, only when stratification variables are of the continuous type. It differentiates from 'optimizeStrata' that accepts both continuous and categorical variables
}
\usage{
	optimizeStrata2(
            errors, 
            framesamp,
            framecens = NULL, 
            strcens = FALSE, 
            model = NULL, 
            alldomains = TRUE, 
            dom = NULL, 
            nStrata = c(5), 
            minnumstr = 2, 
            iter = 50, 
            pops = 20, 
            mut_chance = NA, 
            elitism_rate = 0.2, 
            highvalue = 1e+08, 
            suggestions = NULL, 
            realAllocation = TRUE, 
            writeFiles = FALSE, 
            showPlot = TRUE, 
            parallel = TRUE, 
            cores = NA
	)
}
\arguments{
  \item{errors}{
  This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
  related to target variables of the survey.
}
  \item{framesamp}{
  This is the (mandatory) dataframe containing the information related to the sampling frame.
}
  \item{framecens}{
  This the (optional) dataframe containing the units to be selected in any case. It has same structure than "frame" dataframe.
}
  \item{strcens}{
  Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
}
\item{model}{
  In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
  to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
  target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters 
  beta, sig2, gamma (> 1 in case of heteroscedasticity). 
  Default is NULL.
}
  \item{alldomains}{
  Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
  then a value must be given to parameter 'dom'.
}
  \item{dom}{
  Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval 
  (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
}
  \item{nStrata}{
  Indicates the number of strata for each variable.
}
  \item{minnumstr}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
}
  \item{iter}{
  Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
}
  \item{pops}{
  The dimension of each generations in terms of individuals. Default is 20.
}
  \item{mut_chance}{
  Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
  of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
  but a slower convergence, while low values permit a faster convergence, but the final solution can be
  distant from the optimal one. Default is NA, in correspondence of which it is computed as 1/(vars+1) 
  where vars is the length of elements in the solution.
}
  \item{elitism_rate}{
  This  parameter indicates the rate of better solutions that must be preserved from one generation
	to another. Default is 0.2 (20%).
}
  \item{highvalue}{
  Parameter for genetic algorithm. In should not be changed
}
  \item{suggestions}{
  Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL. 
}
  \item{realAllocation}{
  If FALSE, the allocation is based on INTEGER values; if TRUE, the allocation is based on REAL values.
  Default is TRUE.
}  
  \item{writeFiles}{
  Indicates if the various dataframes and plots produced during the execution have to be written in the working
  directory.
  Default is FALSE.
}
  \item{showPlot}{
  Indicates if the plot showing the trend in the value of the objective function has to be shown or not. In parallel = TRUE, this defaults to FALSE
  Default is TRUE.
}
  \item{parallel}{
  Should the analysis be run in parallel.
  Default is TRUE.
}
  \item{cores}{
  If the analysis is run in parallel, how many cores should be used.
  If not specified n-1 of total available cores are used OR if number of domains < (n-1) cores, then number of cores equal to number of domains are used.
}
}

\value{
A list containing 
(1) the vector of the solution, 
(2) the optimal aggregated strata,
(3) the total sampling frame with the label of aggregated strata
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{
library(SamplingStrata)
data(swissmunicipalities)
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
frame <- buildFrameDF(df = swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("POPTOT","HApoly"),
                      Y = c("Surfacesbois", "Airind"))
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
                         CV1 = rep(0.1,ndom),
                         CV2 = rep(0.1,ndom),
                         domainvalue = c(1:ndom)))
################################
# Example of "continuous" method
################################
kmean <- KmeansSolution2(frame = frame, 
                         errors = cv, 
                         maxclusters = 10)
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
sugg <- prepareSuggestion(kmean = kmean, 
                          frame = frame, 
                          nstrat = nstrat)
solution <- optimizeStrata2(framesamp = frame,
                        errors = cv,
                        nStrata = nstrat,
                        iter = 50,
                        pops = 10,
                        suggestions = sugg)
framenew <- solution$framenew
outstrata <- solution$aggr_strata
s <- selectSample(framenew,outstrata)
}
}
\keyword{ survey }

