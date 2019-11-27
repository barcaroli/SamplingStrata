\name{optimStrata}
\Rdversion{1.5}
\alias{optimStrata}
\title{
Best stratification of a sampling frame for multipurpose surveys 
}
\description{
Wrapper function to call the different optimization functions: (i) optimizeStrata (method = "atomic"); (ii) optimizeStrata2 (method = "continuous"); (iii) optimizeStrata (method = "spatial").
}
\usage{
optimStrata <- function(method=c("atomic","continuous","spatial"),
                        # common parameters
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
                        cores,
                        # parameters only for optimizeStrata
                        strata=NULL,
                        initialStrata=NA,
                        addStrataFactor=NA,
						            cens=NULL,
                        # parameters only for optimizeStrata2 and optimizeStrataSpatial
                        framesamp=NULL,
                        framecens=NULL,
                        model=NULL,
                        nStrata=5,
                        # parameters only for optimizeStrataSpatial
                        fitting=NA,
                        range=NA,
                        kappa=NA
                        ) 
}
\arguments{
  \item{method}{
  This parameter allows to choose the method to be applied in the optimization step: (i) optimizeStrata (method = "atomic"); (ii) optimizeStrata (method = "continuous"); (iii) optimizeStrata (method = "spatial")
}
  \item{errors}{
  This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
  related to target variables of the survey.
}
  \item{alldomains}{
  Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
  then a value must be given to parameter 'dom'.
}
  \item{dom}{
  Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval 
  (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
}
  \item{strcens}{
  Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
}
\item{minnumstr}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
}
  \item{iter}{
  Indicates the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
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
  \item{suggestions}{
  Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL. 
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
  \item{strata}{
  This is the (mandatory) dataframe containing the information related to "atomic" strata, i.e. the strata obtained by 
  the Cartesian product of all auxiliary variables X's. Information concerns the identifiability of strata
  (values of X's) and variability of Y's (for each Y, mean and standard error in strata).
}
  \item{initialStrata}{
  This is the initial limit on the number of strata in the different domains for each solution. Default is NA, and in this case it is set equal to the
  number of atomic strata in each domain.
}
  \item{addStrataFactor}{
  This parameter indicates the probability that at each mutation the number of strata may increase with 
  respect to the current value. Default is 0.0.
}
  \item{cens}{
  This the dataframe containing the takeall strata, those strata whose units must be selected in whatever sample. It has same structure than "strata" dataframe.
}
  \item{framesamp}{
  This is the dataframe containing the information related to the sampling frame.
}
  \item{framecens}{
  This the dataframe containing the units to be selected in any case. It has same structure than "framesamp" dataframe.
}
  \item{nStrata}{
  Indicates the number of strata for each variable.
}
  \item{model}{
  In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
  to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
  target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters 
  beta, sig2, gamma (> 1 in case of heteroscedasticity). 
  Default is NULL.
}
\item{fitting}{
	Fitting of the model(s) (in terms of R squared)
}
\item{range}{
	Maximum range for spatial autocorrelation. It is a vector with as many elements as the number of target variables Y.
}	
\item{kappa}{
	Factor used in evaluating spatial autocorrelation. 
}	

}

\value{
In case of method = "atomic" a list containing 
(1) the vector of the solution, 
(2) the optimal aggregated strata,

In case of method = "continuous" or "spatial" a list containing 
(1) the vector of the solution, 
(2) the optimal aggregated strata,
(3) the total sampling frame with the label of aggregated strata
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{

}
}
\keyword{ survey }

