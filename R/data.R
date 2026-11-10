# ==============================================================================
# File:    data.R
# Purpose: roxygen documentation of the datasets distributed with the package
#          (folder "data"). This file contains no executable code.
# ==============================================================================

#' @title Precision constraints (maximum CVs) as input for Bethel allocation
#'
#' @description
#' Dataframe containing precision levels (expressed in terms of acceptable CV's)
#'
#' @format
#' The constraint data frame (errors) contains a row per each domain value with the following variables:
#' \describe{
#'   \item{DOM}{
#'   Type of domain code (factor)
#'   }
#'   \item{CV1}{
#'   Planned coefficient of variation for first variable Y1 (numeric)
#'   }
#'   \item{CVj}{
#'   Planned coefficient of variation for j-th variable Yj (numeric)
#'   }
#'   \item{CVn}{
#'   Planned coefficient of variation for last variable Yn (numeric)
#'   }
#'   \item{domainvalue}{
#'   Value of the domain to which the constraints refer (numeric)
#'   }
#' }
#'
#' @details
#' Note: the names of the variables must be the ones indicated above
#'
#' @examples
#' ## data(errors)
#' ## errors
#'
#' @keywords datasets
#' @docType data
#' @name errors
NULL

#' @title Dataset 'nations'
#'
#' @description
#' Dataset containing data on 207 countries (from 'nations.txt' in Rcmdr, with missing values imputed)
#'
#' @format
#' A data frame with 207 observations on the following variables:
#' \describe{
#'   \item{Country}{Name of the country}
#'   \item{TFR}{Total Fertility Rate}
#'   \item{contraception}{Rate of women making use of contraceptive means}
#'   \item{infant.mortality}{Infant mortality rate}
#'   \item{GDP}{Gross Domestic Product ($ per capita)}
#'   \item{region}{Continent (description)}
#'   \item{Continent}{Continent (code)}
#' }
#'
#' @source
#' Rcmdr package
#'
#' @keywords datasets
#' @docType data
#' @name nations
NULL

#' @title Dataframe containing information on strata in the frame
#'
#' @description
#' Dataframe containing information on strata in the frame
#'
#' @format
#' The strata data frame (strata) contains a row per stratum with the following variables:
#' \describe{
#'   \item{stratum}{
#' Identifier of the stratum (numeric)
#' }
#'   \item{N}{
#' Number of population units in the stratum (numeric)
#' }
#'   \item{X1}{
#' Value of first auxiliary variable X1 in the stratum (factor)
#' }
#'   \item{Xi}{
#' Value of i-th auxiliary variable Xi in the stratum (factor)
#' }
#'   \item{Xk}{
#' Value of last auxiliary variable Xk in the stratum (factor)
#' }
#'   \item{M1}{
#' Mean in the stratum of the first variable Y1 (numeric)
#' }
#'   \item{Mj}{
#' Mean in the stratum of the j-th variable Yt (numeric)
#' }
#'   \item{Mn}{
#' Mean in the stratum of the last variable Y (numeric)
#' }
#'   \item{S1}{
#' Standard deviation in the stratum of the first variable Y (numeric)
#' }
#'   \item{Sj}{
#' Standard deviation in the stratum of the j-th variable Yt (numeric)
#' }
#'   \item{Sn}{
#' Standard deviation in the stratum of the last variable Y (numeric)
#' }
#'   \item{cens}{
#' Flag (1 indicates a take all straum, 0 a sampling stratum) (numeric)
#' Default = 0
#' }
#'   \item{cost}{
#' Cost per interview in each stratum.  Default = 1 (numeric)
#' }
#'   \item{DOM1}{
#' Value of domain to which the stratum belongs (factor or numeric)
#' }
#' }
#'
#' @details
#' Note: the names of the variables must be the ones indicated above
#'
#' @examples
#' # data(strata)
#' # head(strata)
#'
#' @keywords datasets
#' @docType data
#' @name strata
NULL

#' @title Precision constraints (maximum CVs) as input for Bethel allocation
#'
#' @description
#' Dataframe containing precision levels (expressed in terms of acceptable CV's)
#'
#' @format
#' The constraint data frame (swisserrors) contains a row per each domain value with the following variables:
#' \describe{
#'   \item{DOM}{
#'   Type of domain code (factor)
#'   }
#'   \item{CV1}{
#'   Planned coefficient of variation for first variable Y1 (number of men and women aged between 0 and 19) (numeric)
#'   }
#'   \item{CV2}{
#'   Planned coefficient of variation for second variable Y2 (number of men and women aged between 20 and 39) (numeric)
#'   }
#'   \item{CV3}{
#'   Planned coefficient of variation for third variable Y3 (number of men and women aged between 40 and 64) (numeric)
#'   }
#'   \item{CV4}{
#'   Planned coefficient of variation for forth variable Y4 (number of men and women aged between 65 and over) (numeric)
#'   }
#'   \item{domainvalue}{
#'   Value of the domain to which the constraints refer (numeric)
#'   }
#' }
#'
#' @examples
#' ## data(swisserrors)
#' ## swisserrors
#'
#' @keywords datasets
#' @docType data
#' @name swisserrors
NULL

#' @title Dataframe containing information on all units in the population of reference that can be considered as the final sampling unit (this example is related to Swiss municipalities)
#'
#' @description
#' Dataframe containing information on all municipalities in Swiss
#' (it is a derivation of dataframe "swissmunicipalities" in "sampling" package)
#'
#' @format
#' The "swissframe" dataframe contains a row per each Swiss municipality with the following variables:
#' \describe{
#'   \item{progr}{
#' Progressive associated to the frame unit (numeric)
#' }
#' \item{id}{
#' Name of the frame unit (character)
#' }
#'   \item{X1}{
#' Classes of total population in the municipality (factor with 18 values)
#' }
#'   \item{X2}{
#' Classes of wood area in the municipality (factor with 3 values)
#' }
#'   \item{X3}{
#' Classes of area under cultivation in the municipality (factor with 3 values)
#' }
#' \item{X4}{
#' Classes of mountain pasture area in the municipality (factor with 3 values)
#' }
#' \item{X5}{
#' Classes of area with buildings in the municipality (factor with 3 values)
#' }
#' \item{X6}{
#' Classes of industrial area in the municipality (factor with 3 values)
#' }
#'   \item{Y1}{
#' Number of men and women aged between 0 and 19 (numeric)
#' }
#'   \item{Y2}{
#' Number of men and women aged between 20 and 39 (numeric)
#' }
#'   \item{Y3}{
#' Number of men and women aged between 40 and 64 (numeric)
#' }
#' \item{Y4}{
#' Number of men and women aged between 65 and over (numeric)
#' }
#'   \item{domainvalue}{
#' Value of domain to which the municipality belongs (factor or numeric)
#' }
#' }
#'
#' @examples
#' #data(swissframe)
#' #head(strata)
#'
#' @keywords datasets
#' @docType data
#' @name swissframe
NULL

#' @title The Swiss municipalities population
#'
#' @description
#' This population provides information about the Swiss municipalities
#' in 2003.
#'
#' @format
#' A data frame with 2896 observations on the following 22 variables:
#' \describe{
#' \item{id}{Municipality unique identifier.}
#'   \item{CT}{Swiss canton.}
#'   \item{REG}{Swiss region.}
#'   \item{COM}{municipality number.}
#'   \item{Nom}{municipality name.}
#'   \item{HApoly}{municipality area.}
#'   \item{Surfacesbois}{wood area.}
#'   \item{Surfacescult}{area under cultivation.}
#'   \item{Alp}{mountain pasture area.}
#'   \item{Airbat}{area with buildings.}
#'   \item{Airind}{industrial area.}
#'   \item{P00BMTOT}{number of men.}
#'   \item{P00BWTOT}{number of women.}
#'   \item{Pop020}{number of men and women aged between 0 and 19.}
#'   \item{Pop2040}{number of men and women aged between 20 and 39.}
#'   \item{Pop4065}{number of men and women aged between 40 and 64.}
#'   \item{Pop65P}{number of men and women aged between 65 and over.}
#'   \item{H00PTOT}{number of households.}
#'   \item{H00P01}{number of households with 1 person.}
#'   \item{H00P02}{number of households with 2 persons.}
#'   \item{H00P03}{number of households with 3 persons.}
#'   \item{H00P04}{number of households with 4 persons.}
#'   \item{POPTOT}{total population.}
#'   \item{lat}{latitude.}
#'   \item{long}{longitude.}
#' }
#'
#' @source
#' Swiss Federal Statistical Office.
#'
#' @examples
#' # data(swissmunicipalities)
#' # hist(swissmunicipalities$POPTOT)
#'
#' @keywords datasets
#' @docType data
#' @name swissmunicipalities
NULL

#' @title Dataframe containing information on strata in the swiss municipalities frame
#'
#' @description
#' Dataframe containing information on strata in the swiss municipalities frame
#'
#' @format
#' The "swissstrata" dataframe contains a row per stratum with the following variables:
#' \describe{
#'   \item{STRATO}{
#' Identifier of the stratum (character)
#' }
#'   \item{N}{
#' Number of population units in the stratum (numeric)
#' }
#'   \item{X1}{
#' Classes of total population in the municipality (factor with 18 values)
#' }
#'   \item{X2}{
#' Classes of wood area in the municipality (factor with 3 values)
#' }
#'   \item{X3}{
#' Classes of area under cultivation in the municipality (factor with 3 values)
#' }
#' \item{X4}{
#' Classes of mountain pasture area in the municipality (factor with 3 values)
#' }
#' \item{X5}{
#' Classes of area with buildings in the municipality (factor with 3 values)
#' }
#' \item{X6}{
#' Classes of industrial area in the municipality (factor with 3 values)
#' }
#'   \item{M1}{
#' Mean in the stratum of Y1 (number of men and women aged between 0 and 19)(numeric)
#' }
#'   \item{M2}{
#' Mean in the stratum of Y2 (number of men and women aged between 20 and 39) (numeric)
#' }
#'   \item{M3}{
#' Mean in the stratum of Y3 (number of men and women aged between 40 and 64) (numeric)
#' }
#' \item{M4}{
#' Mean in the stratum of Y4 (number of men and women aged between 64 and over) (numeric)
#' }
#'   \item{S1}{
#' Standard deviation in the stratum of Y1 (number of men and women aged between 0 and 19)(numeric)
#' }
#'   \item{S2}{
#' Standard deviation in the stratum of Y2 (number of men and women aged between 20 and 39) (numeric)
#' }
#'   \item{S3}{
#' Standard deviation in the stratum of Y3 (number of men and women aged between 40 and 64) (numeric)
#' }
#' \item{S4}{
#' Standard deviation in the stratum of Y4 (number of men and women aged between 64 and over) (numeric)
#' }
#'   \item{cens}{
#' Flag (1 indicates a take all straum, 0 a sampling stratum) (numeric)
#' Default = 0
#' }
#'   \item{cost}{
#' Cost per interview in each stratum.  Default = 1 (numeric)
#' }
#'   \item{DOM1}{
#' Value of domain to which the stratum belongs Default = 1 (factor or numeric)
#' }
#' }
#'
#' @examples
#' # data(swissstrata)
#' # head(swissstrata)
#'
#' @keywords datasets
#' @docType data
#' @name swissstrata
NULL
