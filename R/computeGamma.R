# Estimation of the coefficient of heteroscedasticity and variability 

# parameters
# e : errors
# x : explanatory variable
# nbins : number of bins in which categorize errors and explanatory variable

computeGamma <- function (e,
                          x,
                          nbins=6) {
  dataset <- as.data.frame(list(e=e,x=x))
  colnames(dataset)
  require(SamplingStrata)
  ris_gamma <- NULL
  dataset$p_bins <- var.bin(dataset$x,nbins)
  std_eps <- sqrt(tapply(dataset$e,dataset$p_bins,var))
  boxplot(e ~ p_bins, data=dataset) 
  x_eps <- tapply(dataset$x,dataset$p_bins,mean)
  plot(x_eps,std_eps)
  plot(log(x_eps),log(std_eps))
  lm_gamma <- lm(log(std_eps)~log(x_eps))
  lm_gamma$coefficients
  gamma_est <- lm_gamma$coefficients[2]
  sigma <- exp(lm_gamma$coefficients[1])
  gamma_est
  ris_gamma[1] <- gamma_est
  ris_gamma[2] <- sigma
  return(ris_gamma)
}



