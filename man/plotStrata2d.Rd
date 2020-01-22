\name{plotStrata2d}
\alias{plotStrata2d}
\title{Plot bivariate distibutions in strata}
\description{
  Plots a 2d graph showing obtained strata 
}
\usage{
plotStrata2d (x,outstrata,domain,vars,labels)
}
\arguments{
  \item{x}{the sampling frame}
  \item{outstrata}{the optimized strata}
  \item{domain}{a domain in the frame}
  \item{vars}{vars to appear in x and y axis}
  \item{labels}{labels to appear in x and y axis}
 }
\value{
A formatted output containing information on the strata in the given domain
}
\examples{
\dontrun{
library(SamplingStrata)
data("swissmunicipalities")
swissmunicipalities = swissmunicipalities[swissmunicipalities$REG==1,]
data("errors")
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
swissmunicipalities$domain = 1
frame <- buildFrameDF(swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("Surfacesbois","Surfacescult"),
                      Y = c("Pop020", "Pop2040")
)
solution <- optimStrata (method = "continuous",
                        errors = errors,, 
                        framesamp = frame,
                        nStrata = 8, 
                        iter = 25, 
                        pops = 10) 
p <- plotStrata2d(solution$framenew, 
                  solution$aggr_strata,
                  domain = 1, 
                  vars = c("X1","X2"),
                  labels = c("Surfacesbois","Surfacescult"))
p
}
}
\keyword{survey}
