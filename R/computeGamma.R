# Estimation of the coefficient of heteroscedasticity and residuals variability 

# parameters
# e : residuals
# x : explanatory variable
# nbins : number of bins in which categorize errors and explanatory variable

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



