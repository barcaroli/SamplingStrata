# SamplingStrata 1.3

Changes in Version 1.3

  o A new function 'KmeansSolution' produces an initial solution using the kmeans algorithm by clustering atomic strata considering the values of the means of target variables in them. Also, if the parameter 'nstrata' is not indicated, the optimal number of clusters is determined inside each domain, and the overall solution is obtained by concatenating optimal clusters obtained in domains. By indicating this solution as a suggestion to the optimization step, this may greatly speed the convergence to the optimal solution.

  o A new function 'selectSampleSystematic' has been added. It allows to select a stratified sample with the systematic method, that is a selection that begins   selecting the first unit by an initial randomly chosen starting point, and proceeding in selecting other units by adding an interval that is the inverse of the sampling rate in the stratum. This selection method can be useful if associated to a particular ordering of the selection frame, where the ordering variable(s) can be considered as additional stratum variable(s). 
  
  o It is now possible to handle "anticipated variance" by introducing a model linking a proxy variable whose values are available for all units in the sampling frame, with the target variable whose values are not available. In this implementation only linear model can be chosen. When calling the 'buildStrataDF' function, a dataframe is given, containing three parameters (beta, sigma2 and gamma)
for each couple target / proxy. On the basis of these parameters, means and   standard deviations in sampling strata are calculated accordingly to given formulas.





