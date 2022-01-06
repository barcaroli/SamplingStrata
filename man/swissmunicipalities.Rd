\name{swissmunicipalities}
\alias{swissmunicipalities}
\docType{data}
\title{The Swiss municipalities population}
\description{This population provides information about the Swiss municipalities 
in 2003.
}
\usage{data(swissmunicipalities)}
\format{
  A data frame with 2896 observations on the following 22 variables:
  \describe{
	\item{id}{Municipality unique identifier.}
    \item{CT}{Swiss canton.}
    \item{REG}{Swiss region.}
    \item{COM}{municipality number.}
    \item{Nom}{municipality name.}
    \item{HApoly}{municipality area.}
    \item{Surfacesbois}{wood area.}
    \item{Surfacescult}{area under cultivation.}
    \item{Alp}{mountain pasture area.}
    \item{Airbat}{area with buildings.}
    \item{Airind}{industrial area.}
    \item{P00BMTOT}{number of men.}
    \item{P00BWTOT}{number of women.}
    \item{Pop020}{number of men and women aged between 0 and 19.}
    \item{Pop2040}{number of men and women aged between 20 and 39.}
    \item{Pop4065}{number of men and women aged between 40 and 64.}
    \item{Pop65P}{number of men and women aged between 65 and over.}
    \item{H00PTOT}{number of households.}
    \item{H00P01}{number of households with 1 person.}
    \item{H00P02}{number of households with 2 persons.}
    \item{H00P03}{number of households with 3 persons.}
    \item{H00P04}{number of households with 4 persons.}
    \item{POPTOT}{total population.}
    \item{lat}{latitude.}
    \item{long}{longitude.}
  }
}
\source{Swiss Federal Statistical Office.
}
\examples{
# data(swissmunicipalities)
# hist(swissmunicipalities$POPTOT)
}
\keyword{datasets}



