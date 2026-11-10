# ==============================================================================
# File:      computeGamma.R
# Functions: computeGamma() [exported]
#
# Purpose:
#   Estimates the heteroscedasticity of the residuals e of a model with
#   respect to the predicted values x: x is binned into nbins classes
#   (var.bin()), and log(sd of e) is regressed on log(mean of x) across the
#   classes. Returns gamma (slope), sigma (exp of the intercept) and the R^2 of
#   this regression, to be used in the "model" dataframe for the anticipated
#   variance in the strata. Optionally plots residuals and fitted model.
# ==============================================================================

# Estimation of the coefficient of heteroscedasticity and residuals variability 

# parameters
# e : residuals
# x : explanatory variable
# nbins : number of bins in which categorize errors and explanatory variable

#' @title Function that allows to calculate a heteroscedasticity index, together with associate prediction variance, to be used by the optimization step to correctly evaluate the standard deviation in the strata due to prediction errors.
#'
#' @description
#' When the anticipated variance has to be calculated during the execution
#' of the optimization step, his function allows to calculate a heteroscedasticity index,
#' together with associate prediction variance,
#' to be used to correctly evaluate the variance in the strata.
#' The function returns a list where the first object is the heteroscedasticity index and the
#' second is the associated standard deviation in the strata
#' due to prediction errors.
#' The two parameters are calculated in this way:
#' (i) residuals 'e' are grouped in clusters defined by values of the explanatory variable 'x';
#' (ii) a model is fitted by considering log(e) and log(mean(x)) values;
#' (iii) the intercept is the value of standard deviation of residuals;
#' (iv) the slope is the value of the heteroscedasticity index.
#' These two values can be passed as parameters of the model, or used to calculate
#' prediction errors for ach unit in the frame.
#'
#' @param e This is the variable that contains prediction errors (residuals) of the model.
#' @param x This is the variable that contains explanatory variable in the model.
#' @param nbins Number of bins to be passed to the 'var.bin' function. Default is 6.
#' @param showPlot Visualization of plots. Default is TRUE.
#'
#' @return
#' A list containing: (i) the value of the heteroscedasticity index, (ii) associated standard deviation, (iii) R^2 of the interpolating model.
#'
#' @author
#' Marco Ballin, Giulio Barcaroli
#'
#' @examples
#' \dontrun{
#' data("swissmunicipalities")
#' swiss_sample <- swissmunicipalities[sample(c(1:nrow(swissmunicipalities)),500),]
#' mod_Airind_POPTOT <- lm(swiss_sample$Airind ~ swiss_sample$POPTOT)
#' computeGamma(mod_Airind_POPTOT$residuals,
#'              swiss_sample$POPTOT,
#'              nbins = 8)
#' #     gamma     sigma  r.square
#' # 0.8029292 0.0150446 0.9598539
#' }
#'
#' @keywords survey
#' @export
computeGamma <- function (e, x, nbins = 6, showPlot = TRUE) 
{
  dataset <- as.data.frame(list(e = e, x = x), stringsAsFactors = TRUE)
  colnames(dataset)
  require(SamplingStrata)
  ris_gamma <- NULL
  dataset$p_bins <- var.bin(dataset$x, nbins)
  std_eps <- sqrt(tapply(dataset$e, dataset$p_bins, var))
  x_eps <- tapply(dataset$x, dataset$p_bins, mean)
  lm_gamma <- lm(log(std_eps) ~ log(x_eps))
  summary(lm_gamma)$r.square
  if (showPlot == TRUE) {
    par(mfrow = c(1, 2))
    boxplot(e ~ p_bins, data = dataset, ylab = "residuals", 
            xlab = "clusters", col = "orange")
    title("Distribution of residuals")
    plot(log(x_eps), log(std_eps), xlab = "log(mean(x))", 
         ylab = "log(stdev(residuals))")
    abline(lm_gamma)
    title("Model")
    par(mfrow = c(1, 1))
  }
  ris_gamma[1] <- lm_gamma$coefficients[2]
  ris_gamma[2] <- exp(lm_gamma$coefficients[1])
  ris_gamma[3] <- summary(lm_gamma)$r.square
  names(ris_gamma) <- c("gamma","sigma","r.square")
  return(ris_gamma)
}



