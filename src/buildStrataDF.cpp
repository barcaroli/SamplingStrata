// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <set>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Funzioni di supporto per il calcolo della deviazione standard

double stdev1(const std::vector<double>& x, const std::vector<double>& w) {
  double sumw = 0.0, sumxw = 0.0;
  for (size_t i = 0; i < x.size(); i++){
    sumw += w[i];
    sumxw += x[i] * w[i];
  }
  double mean = sumxw / sumw;
  double sumsq = 0.0;
  for (size_t i = 0; i < x.size(); i++){
    sumsq += w[i] * std::pow(x[i] - mean, 2);
  }
  return std::sqrt(sumsq / (sumw - 1));
}

double stdev2(const std::vector<double>& x, const std::vector<double>& w) {
  double sumw = 0.0, sumxw = 0.0;
  for (size_t i = 0; i < x.size(); i++){
    sumw += w[i];
    sumxw += x[i] * w[i];
  }
  double mean = sumxw / sumw;
  double sumsq = 0.0;
  for (size_t i = 0; i < x.size(); i++){
    sumsq += w[i] * std::pow(x[i] - mean, 2);
  }
  return std::sqrt(sumsq / sumw);
}

// (stdev3, the unweighted standard deviation with divisor n - 1 of
// beta1 * Y + beta2 * W, was removed: the "spatial" model now computes the
// standard deviation of the predictions with stdev1 / stdev2, as the "linear"
// model does for Y)

// stdev4: spatial component of the anticipated variance for the "spatial"
// model type (pairs i < j only, parallelised if OpenMP is available).
// Error variance of unit i: s2_i = var_eps * Y_i^(2 gamma).
// Covariance of the errors of two distinct units:
//   Cov(e_i, e_j) = (1 - nugget) * s_i * s_j * exp(-kappa * d_ij / range)
// (no correlation with range <= 0 or nugget >= 1).
// Returns sqrt( sum_{i != j} [s2_i + s2_j - 2 Cov(e_i, e_j)] / (2 n^2) ).
double stdev4(const std::vector<double>& Y,
              const std::vector<double>& LON,
              const std::vector<double>& LAT,
              double var_eps,
              double range,
              double gamma,
              double kappa,
              double nugget) {
  size_t n = Y.size();
  if(n == 0) return 0.0;
  std::vector<double> var_ntimes(n);
  for (size_t i = 0; i < n; i++){
    var_ntimes[i] = var_eps * std::pow(Y[i], 2 * gamma);
  }
  bool spatial = (range > 0.0) && (nugget < 1.0);
  double inv_r = spatial ? kappa / range : 0.0;
  double w = 1.0 - nugget;
  double sumD2 = 0.0;
  // Itera solo su coppie con j > i
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sumD2) schedule(static)
#endif
  for (size_t i = 0; i < n; i++){
    for (size_t j = i+1; j < n; j++){
      double dlon = LON[i] - LON[j];
      double dlat = LAT[i] - LAT[j];
      double d = std::sqrt(dlon * dlon + dlat * dlat);
      double sum_couples_var = var_ntimes[i] + var_ntimes[j];
      double prod_std = std::sqrt(var_ntimes[i] * var_ntimes[j]);
      double spatial_autocov = spatial ? w * prod_std * std::exp(- d * inv_r) : 0.0;
      double D2 = sum_couples_var - 2 * spatial_autocov;
      sumD2 += D2;
    }
  }
  // Poiché abbiamo sommato solo le coppie (i,j) con i < j, moltiplichiamo per 2
  double var2 = (2 * sumD2) / (2 * n * n); // equivalente a sumD2/(n*n)
  return std::sqrt(var2);
}

// Struttura per raccogliere i risultati per ciascun dominio
struct DomainResult {
  std::vector<std::string> STRATO;
  std::vector<double> N;
  std::vector< std::vector<double> > M; // dimensione: nvarY
  std::vector< std::vector<double> > S; // dimensione: nvarY
  std::vector<double> COST;
  std::vector<double> CENS;
  std::vector<int> DOM1;            // MODIFICA: ora è un vettore di int
  std::vector< std::vector<double> > X; // dimensione: nvarX
};

//' @title Builds the "strata" dataframe containing information on target variables Y's distributions in the different strata, starting from sample data or from a frame
//'
//' @description
//' This function allows to build the information regarding strata in the population required as
//' an input by the algorithm of Bethel for the optimal allocation.
//' In order to estimate means and standard deviations for target variables Y's, we need data coming from:
//' (1) a previous round of the survey whose sample we want to plan;
//' (2) sample data from a survey with variables that are proxy to the ones we are interested to;
//' (3) a frame containing values of Y's variables (or proxy variables) for all the population.
//' In all cases, each unit in the dataset must contain auxiliary information (X's variables)
//' and also target variables Y's (or proxy variables) values: under these conditions it is possible
//' to build the dataframe "strata", containing information on the distribution of Y's in the different strata
//' (namely, means and standard deviations), together with information on strata (total population,
//' if it is to be censused or not, the cost per single interview).
//' If the information is contained in a sample dataset, a variable named WEIGHT is expected to be
//' present. In case of a frame, no such variable is given, and the function will define a WEIGHT variable
//' for each unit, whose value is always '1'.
//' Missing values for each Y variable will not be taken into account in the computation of means and standard
//' deviations (in any case, NA's can be present in the dataset).
//' The dataframe "strata" is written to an external file (tab delimited, extension "txt"), and will be
//' used as an input by the function "optimizeStrata".
//'
//' @param dataset This is the name of the dataframe containing the sample data, or the frame data.
//'   It is strictly required that auxiliary information is organised in variables named
//'   as X1, X2, ... , Xm (there should be at least one of them) and the target variables
//'   are denoted by Y1, Y2, ... , Yn.
//'   In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe,
//'   containing   the weigths associated to each sampling unit
//' @param model_ In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
//'   to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
//'   target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters
//'   beta, sig2, gamma (> 1 in case of heteroscedasticity).
//'   For type "spatial" (target modelled as beta * Y + beta2 * W + e, with error variance sig2 * Y^(2 gamma))
//'   the variance in a stratum is s1^2 / fitting + s2^2, where s1^2 is the variance of the
//'   predictions (divisor N for frame data, sum of weights - 1 for sample data) and s2^2 the mean, over the pairs of distinct units, of
//'   s2_i + s2_j - 2 (1 - nugget) s_i s_j exp(-kappa d_ij / range). Required columns: beta2, range,
//'   gamma. Optional columns: fitting (default 1), nugget (relative nugget c0/(c0+c1) of the
//'   variogram of the errors, default 0), kappa (default 3, i.e. range is the practical range of an
//'   exponential variogram; 1 if it is the range parameter, as returned by gstat).
//'   Since the error variance already compensates the levelling of the predictions, fitting must be 1
//'   when sig2 > 0 (additive correction, de Gruijter et al. 2016): fitting = R2 < 1 together with
//'   sig2 > 0 would count the error variance twice, and is rejected.
//'   Default is NULL.
//' @param progress If set to TRUE, a progress bar is visualised during the execution. Default is TRUE.
//' @param verbose If set to TRUE, information is given about the number of strata generated. Default is TRUE.
//'
//' @return
//' A dataframe containing strata
//'
//' @author
//' Giulio Barcaroli
//'
//' @examples
//' \dontrun{
//' # Plain example without model
//' data(swissframe)
//' strata <- buildStrataDF(dataset=swissframe,model=NULL)
//' head(strata)
//'
//' # More complex example with models
//' library(SamplingStrata)
//' data(swissmunicipalities)
//' swiss <- swissmunicipalities[,c("HApoly","Surfacesbois","Airind","POPTOT")]
//' Y1 = swiss$Surfacesbois
//' X1 = swiss$HApoly
//' mod1 <- lm( Y1 ~ X1 )
//' summary(mod1)
//' mod1$coefficients[2]
//' summary(mod1)$sigma
//'
//' Y2 = swiss$Airind
//' X2 = swiss$POPTOT
//' plot(log(X2[X2>0]),log(Y2[X2>0]))
//' mod2 <- lm( log(Y2[X2 > 0 & Y2>0]) ~ log(X2[X2 > 0 & Y2>0]) )
//' summary(mod2)
//' mod2$coefficients[2]
//' summary(mod2)$sigma
//'
//' swiss$id <- c(1:nrow(swiss))
//' swiss$dom <- 1
//' frame <- buildFrameDF(swiss,id="id",X="id",Y=c("HApoly","POPTOT"),domainvalue="dom")
//'
//' model <- NULL
//' model$type[1] <- "linear"
//' model$beta[1] <- mod1$coefficients[2]
//' model$sig2[1] <- summary(mod1)$sigma
//' model$gamma[1] = 2
//' model$type[2] <- "loglinear"
//' model$beta[2] <- mod2$coefficients[2]
//' model$sig2[2] <- summary(mod2)$sigma
//' model$gamma[2] = NA
//' model <- as.data.frame(model)
//'
//' strata <- buildStrataDF(dataset=frame, model=model)
//' }
//'
//' @keywords survey
//' @export
// [[Rcpp::export]]
DataFrame buildStrataDF(DataFrame dataset, Nullable<DataFrame> model_ = R_NilValue, 
                             bool progress = true, bool verbose = true) {
  // Upper-case column names on a COPY of the names and of the column list:
  // modifying dataset.names() in place changed the names of the caller's data
  // frame (and of every R object sharing it). The columns themselves are not
  // copied (shallow list of pointers).
  int nrows = dataset.nrows();
  CharacterVector colNames = clone(as<CharacterVector>(dataset.names()));
  for (int i = 0; i < colNames.size(); i++){
    std::string name = as<std::string>(colNames[i]);
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    colNames[i] = name;
  }
  List ds(dataset.size());
  for (int i = 0; i < dataset.size(); i++) ds[i] = dataset[i];
  ds.attr("names") = colNames;
  
  // Conta le colonne che iniziano per X e Y
  int nvarX = 0, nvarY = 0;
  for (int i = 0; i < colNames.size(); i++){
    std::string name = as<std::string>(colNames[i]);
    if(!name.empty()){
      if(name[0]=='X') nvarX++;
      if(name[0]=='Y') nvarY++;
    }
  }
  
  // Estrazione della colonna WEIGHT (o default a 1)
  bool weightProvided = false;
  NumericVector weight;
  if (ds.containsElementNamed("WEIGHT"))
    weight = ds["WEIGHT"], weightProvided = true;
  else {
    weight = NumericVector(nrows, 1.0);
    weightProvided = false;
  }
  
  // Gestione del modello
  bool useModel = false;
  DataFrame model;
  if(model_.isNotNull()){
    model = model_.get();
    useModel = true;
    // fitting < 1 together with an error variance counts it twice
    if(model.containsElementNamed("fitting") && model.containsElementNamed("type") &&
       model.containsElementNamed("sig2")) {
      CharacterVector tp = model["type"];
      NumericVector ft = model["fitting"];
      NumericVector sg = model["sig2"];
      for (int k = 0; k < model.nrows(); k++) {
        if (as<std::string>(tp[k]) == "spatial" && !NumericVector::is_na(ft[k]) &&
            ft[k] < 1.0 && !NumericVector::is_na(sg[k]) && sg[k] > 0.0) {
          stop("model of type 'spatial': fitting < 1 cannot be used together with sig2 > 0 "
               "(the error variance would be counted twice). Omit fitting (default 1).");
        }
      }
    }
    if(model.nrows() != nvarY)
      stop("A model for each Y variable must be specified");
    // Validation of the model parameters. It must be done here, before the
    // OpenMP parallel region: an Rcpp::stop() thrown inside the parallel loop
    // cannot be caught by R and terminates the session.
    if(model.containsElementNamed("type")) {
      CharacterVector tp = model["type"];
      for (int k = 0; k < model.nrows(); k++) {
        std::string tk = as<std::string>(tp[k]);
        if (tk != "linear" && tk != "loglinear" && tk != "spatial")
          stop("Type of model misspecified");
        if (tk != "spatial") continue;
        if (!model.containsElementNamed("range") || !model.containsElementNamed("gamma"))
          stop("Missing parameters in model for spatial type (range, gamma)");
        if (!model.containsElementNamed("beta2"))
          stop("beta2 missing in model for spatial type");
        if (model.containsElementNamed("fitting")) {
          NumericVector v = model["fitting"];
          if (!NumericVector::is_na(v[k]) && v[k] <= 0.0)
            stop("fitting must be > 0 in model for spatial type");
        }
        if (model.containsElementNamed("nugget")) {
          NumericVector v = model["nugget"];
          if (!NumericVector::is_na(v[k]) && (v[k] < 0.0 || v[k] > 1.0))
            stop("nugget must be in [0, 1] in model for spatial type");
        }
        if (model.containsElementNamed("kappa")) {
          NumericVector v = model["kappa"];
          if (!NumericVector::is_na(v[k]) && (!R_finite(v[k]) || v[k] <= 0.0))
            stop("kappa must be a finite value > 0 in model for spatial type");
        }
      }
    }
  }
  
  if(!ds.containsElementNamed("DOMAINVALUE"))
    stop("DOMAINVALUE column missing in dataset");
  CharacterVector domainValues = ds["DOMAINVALUE"];
  
  // Raccolta dei domini unici
  std::set<std::string> domainsSet;
  for (int i = 0; i < domainValues.size(); i++){
    std::string dom = as<std::string>(domainValues[i]);
    domainsSet.insert(dom);
  }
  std::vector<std::string> domains(domainsSet.begin(), domainsSet.end());
  // numeric order of the domains (std::set orders strings alphabetically:
  // "1","10","11","2",...), as expected by the R functions that assign the
  // solution of each domain to the strata (e.g. updateStrata())
  std::sort(domains.begin(), domains.end(),
            [](const std::string &a, const std::string &b) { return std::stod(a) < std::stod(b); });
  
  // ---------------------------------------------------------------------------
  // Everything the domain loop needs is extracted here, in the main thread, into
  // plain C++ containers: the loop below may run in an OpenMP parallel region,
  // where R/Rcpp API calls (allocation, as<>(), stop(), Rcout) are not safe.
  // All the checks on the data and on the model are therefore done here.
  // ---------------------------------------------------------------------------
  std::vector<std::string> domStr(nrows);
  for (int i = 0; i < nrows; i++) domStr[i] = as<std::string>(domainValues[i]);

  std::vector< std::vector<double> > Xcols;
  for (int j = 1; j <= nvarX; j++){
    std::string colName = "X" + std::to_string(j);
    if(!ds.containsElementNamed(colName.c_str()))
      stop(("Column " + colName + " missing in dataset").c_str());
    Xcols.push_back(as< std::vector<double> >(ds[colName.c_str()]));
  }
  std::vector< std::vector<double> > Ycols;
  for (int j = 1; j <= nvarY; j++){
    std::string colName = "Y" + std::to_string(j);
    if(!ds.containsElementNamed(colName.c_str()))
      stop(("Column " + colName + " missing in dataset").c_str());
    Ycols.push_back(as< std::vector<double> >(ds[colName.c_str()]));
  }
  std::vector<double> w_all = as< std::vector<double> >(weight);
  bool hasCost = ds.containsElementNamed("COST");
  std::vector<double> cost_all;
  if(hasCost) cost_all = as< std::vector<double> >(ds["COST"]);

  // model parameters, one value per Y variable
  std::vector<std::string> m_type(nvarY);
  std::vector<double> m_beta(nvarY, 0.0), m_sig2(nvarY, 0.0), m_gamma(nvarY, 0.0),
                      m_beta2(nvarY, 0.0), m_range(nvarY, 0.0),
                      m_fitting(nvarY, 1.0), m_nugget(nvarY, 0.0), m_kappa(nvarY, 3.0);
  std::vector< std::vector<double> > Wcols(nvarY);
  std::vector<double> lon_all, lat_all;
  if(useModel){
    if(!model.containsElementNamed("type") || !model.containsElementNamed("beta") ||
       !model.containsElementNamed("sig2"))
      stop("model must contain the columns type, beta and sig2");
    CharacterVector tp = model["type"];
    NumericVector bt = model["beta"];
    NumericVector sg = model["sig2"];
    bool anySpatial = false;
    for (int k = 0; k < nvarY; k++) {
      m_type[k] = as<std::string>(tp[k]);
      m_beta[k] = bt[k];
      m_sig2[k] = sg[k];
      if (m_type[k] == "linear") {
        if(!model.containsElementNamed("gamma"))
          stop("gamma missing in model for linear type");
        NumericVector v = model["gamma"]; m_gamma[k] = v[k];
      } else if (m_type[k] == "spatial") {
        // presence of range, gamma, beta2 and validity of fitting, nugget,
        // kappa have been checked above
        anySpatial = true;
        NumericVector g = model["gamma"];  m_gamma[k] = g[k];
        NumericVector r = model["range"];  m_range[k] = r[k];
        NumericVector b2 = model["beta2"]; m_beta2[k] = b2[k];
        if(model.containsElementNamed("fitting")) {
          NumericVector v = model["fitting"];
          if(!NumericVector::is_na(v[k])) m_fitting[k] = v[k];
        }
        if(model.containsElementNamed("nugget")) {
          NumericVector v = model["nugget"];
          if(!NumericVector::is_na(v[k])) m_nugget[k] = v[k];
        }
        if(model.containsElementNamed("kappa")) {
          NumericVector v = model["kappa"];
          if(!NumericVector::is_na(v[k])) m_kappa[k] = v[k];
        }
        std::string colW = "W" + std::to_string(k + 1);
        if(ds.containsElementNamed(colW.c_str()))
          Wcols[k] = as< std::vector<double> >(ds[colW.c_str()]);
        else
          Wcols[k] = std::vector<double>(nrows, 0.0);
      }
    }
    if (anySpatial) {
      if(!ds.containsElementNamed("LON") || !ds.containsElementNamed("LAT"))
        stop("Missing coordinates on sampling frame (LON, LAT) for model of spatial type");
      lon_all = as< std::vector<double> >(ds["LON"]);
      lat_all = as< std::vector<double> >(ds["LAT"]);
      for (int i = 0; i < nrows; i++)
        if (ISNAN(lon_all[i]) || ISNAN(lat_all[i]))
          stop("Missing values in the coordinates (LON, LAT) of the sampling frame");
    }
  }

  // Prepara un vettore di risultati per ciascun dominio
  std::vector<DomainResult> domain_results(domains.size());
  
  // Parallelizza sul loop sui domini (ogni dominio viene processato in modo
  // indipendente). No R/Rcpp API calls inside the loop.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int d = 0; d < (int)domains.size(); d++){
    DomainResult local;
    const std::string &dom = domains[d];
    // Indici delle righe per il dominio corrente
    std::vector<int> idx;
    for (int i = 0; i < nrows; i++){
      if(domStr[i] == dom)
        idx.push_back(i);
    }
    // Costruisce lo "strato" concatenando i valori delle colonne X
    std::vector<std::string> strato(idx.size());
    for (size_t k = 0; k < idx.size(); k++){
      int i = idx[k];
      std::ostringstream oss;
      // 17 significant digits: distinct X values give distinct strata and the
      // values parsed back from the key are exact (the default precision of 6
      // digits merged e.g. 1234567 and 1234568)
      oss << std::setprecision(17);
      for (int j = 0; j < nvarX; j++){
        oss << Xcols[j][i];
        if(j < nvarX - 1)
          oss << "*";
      }
      strato[k] = oss.str();
    }
    // Raggruppa per strato
    std::unordered_map<std::string, std::vector<int>> group;
    for (size_t k = 0; k < idx.size(); k++){
      group[strato[k]].push_back(idx[k]);
    }
    
    // Prepara i vettori locali per i risultati: M e S sono vettori di vettori con dimensione nvarY
    local.M.resize(nvarY);
    local.S.resize(nvarY);
    local.X.resize(nvarX);
    
    // Itera sui gruppi (strati) del dominio corrente
    for (auto const &g : group) {
      const std::string &strata_id = g.first;
      const std::vector<int> &group_indices = g.second;
      double sumW = 0.0;
      double sumCost = 0.0;
      for (int i : group_indices) {
        sumW += w_all[i];
        if(hasCost)
          sumCost += w_all[i] * cost_all[i];
      }
      double N_val = sumW;
      double cost_val = hasCost ? (sumCost / sumW) : 1.0;
      
      std::vector<double> M_vals(nvarY, 0.0);
      std::vector<double> S_vals(nvarY, 0.0);
      
      // Calcola M e S per ciascuna variabile Y
      for (int y = 0; y < nvarY; y++){
        const std::vector<double> &Ycol = Ycols[y];
        const std::string &mtype = m_type[y];
        bool spatialY = useModel && mtype == "spatial";
        // units with missing Y are skipped; for the spatial model W and the
        // coordinates are taken for the same units, so that they stay aligned
        std::vector<double> Y_vals, w_vals, W_vals, lon_vals, lat_vals;
        for (int i : group_indices) {
          double y_val = Ycol[i];
          if(!ISNAN(y_val)) {
            Y_vals.push_back(y_val);
            w_vals.push_back(w_all[i]);
            if(spatialY) {
              W_vals.push_back(Wcols[y][i]);
              lon_vals.push_back(lon_all[i]);
              lat_vals.push_back(lat_all[i]);
            }
          }
        }
        
        if(Y_vals.size() == 0){
          M_vals[y] = 0.0;
          S_vals[y] = 0.0;
          continue;
        }
        
        if(!useModel) {
          double sumYW = 0.0, sumw = 0.0;
          for (size_t j = 0; j < Y_vals.size(); j++){
            sumYW += Y_vals[j] * w_vals[j];
            sumw += w_vals[j];
          }
          M_vals[y] = sumYW / sumw;
          S_vals[y] = weightProvided ? stdev1(Y_vals, w_vals) : stdev2(Y_vals, w_vals);
        } else {
          double beta = m_beta[y];
          double sig2 = m_sig2[y];
          if(mtype == "linear") {
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += Y_vals[j] * beta * w_vals[j];
              sumw += w_vals[j];
            }
            M_vals[y] = sumYW / sumw;
            double s_raw = weightProvided ? stdev1(Y_vals, w_vals) : stdev2(Y_vals, w_vals);
            double sum_pow = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sum_pow += std::pow(Y_vals[j], 2 * m_gamma[y]);
            }
            double gammas = sum_pow / Y_vals.size();
            S_vals[y] = std::sqrt(s_raw * s_raw * beta * beta + sig2 * gammas);
          } else if(mtype == "loglinear") {
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += std::pow(Y_vals[j], beta) * w_vals[j];
              sumw += w_vals[j];
            }
            M_vals[y] = sumYW / sumw;
            double count_positive = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              if(Y_vals[j] > 0) count_positive += 1.0;
            }
            double ph = count_positive / Y_vals.size();
            double sumY2 = 0.0, sumYbeta = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumY2 += std::pow(Y_vals[j], 2 * beta) * w_vals[j];
              sumYbeta += std::pow(Y_vals[j], beta) * w_vals[j];
            }
            double avgYbeta = sumYbeta / Y_vals.size();
            S_vals[y] = std::sqrt(ph * ((std::exp(sig2) * (sumY2 / Y_vals.size())) - ph * (avgYbeta * avgYbeta)));
          } else { // "spatial" (the type has been checked before the loop)
            double beta2 = m_beta2[y];
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += (Y_vals[j] * beta + W_vals[j] * beta2) * w_vals[j];
              sumw += w_vals[j];
            }
            M_vals[y] = sumYW / sumw;
            // standard deviation of the predictions beta * Y + beta2 * W:
            // divisor sum(w) - 1 for sample data (WEIGHT given), sum(w) for
            // frame data, as for the other model types
            std::vector<double> pred_vals(Y_vals.size());
            for (size_t j = 0; j < Y_vals.size(); j++)
              pred_vals[j] = Y_vals[j] * beta + W_vals[j] * beta2;
            double s1 = weightProvided ? stdev1(pred_vals, w_vals) : stdev2(pred_vals, w_vals);
            double s2 = stdev4(Y_vals, lon_vals, lat_vals, sig2, m_range[y], m_gamma[y],
                               m_kappa[y], m_nugget[y]);
            S_vals[y] = std::sqrt((s1 * s1 / m_fitting[y]) + s2 * s2);
          }
        }
      }
      
      // Salva i risultati locali per lo strato corrente
      // stratum label with 15 significant digits (as R does), e.g. "1*2"
      {
        std::istringstream iss_lab(strata_id);
        std::string tok;
        std::ostringstream lab;
        lab << std::setprecision(15);
        bool first = true;
        while (std::getline(iss_lab, tok, '*')) {
          if (!first) lab << "*";
          lab << std::stod(tok);
          first = false;
        }
        local.STRATO.push_back(lab.str());
      }
      local.N.push_back(N_val);
      for (int j = 0; j < nvarY; j++){
        local.M[j].push_back(M_vals[j]);
        local.S[j].push_back(S_vals[j]);
      }
      local.COST.push_back(cost_val);
      local.CENS.push_back(0.0);
      // MODIFICA: converte il dominio (stringa) in intero
      local.DOM1.push_back(std::stoi(dom));
      // Estrae i valori X dal "strato" (split della stringa)
      std::istringstream iss(strata_id);
      std::string token;
      std::vector<double> Xvals;
      while(std::getline(iss, token, '*')) {
        Xvals.push_back(std::stod(token));
      }
      for (int j = 0; j < nvarX; j++){
        double val = (j < (int)Xvals.size() ? Xvals[j] : 0.0);
        local.X[j].push_back(val);
      }
    } // Fine ciclo per ciascun gruppo
    
    domain_results[d] = local;
  } // Fine ciclo sui domini
  
  // progress messages printed here, in the main thread
  if(progress) {
    for (size_t d = 0; d < domains.size(); d++)
      Rcout << "Processed domain: " << domains[d] << "\n";
  }
  
  // Unione dei risultati di tutti i domini
  std::vector<std::string> out_STRATO;
  std::vector<double> out_N;
  std::vector< std::vector<double> > out_M(nvarY);
  std::vector< std::vector<double> > out_S(nvarY);
  std::vector<double> out_COST;
  std::vector<double> out_CENS;
  std::vector<int> out_DOM1;           // MODIFICA: vettore di int
  std::vector< std::vector<double> > out_X(nvarX);
  
  for (size_t d = 0; d < domain_results.size(); d++){
    DomainResult &local = domain_results[d];
    for (size_t i = 0; i < local.STRATO.size(); i++){
      out_STRATO.push_back(local.STRATO[i]);
      out_N.push_back(local.N[i]);
      for (int j = 0; j < nvarY; j++){
        out_M[j].push_back(local.M[j][i]);
        out_S[j].push_back(local.S[j][i]);
      }
      out_COST.push_back(local.COST[i]);
      out_CENS.push_back(local.CENS[i]);
      out_DOM1.push_back(local.DOM1[i]);  // MODIFICA: copia dei valori interi
      for (int j = 0; j < nvarX; j++){
        out_X[j].push_back(local.X[j][i]);
      }
    }
  }
  
  List out;
  out["STRATO"] = out_STRATO;
  out["N"] = out_N;
  for (int j = 0; j < nvarY; j++){
    std::string name = "M" + std::to_string(j+1);
    out[name] = out_M[j];
  }
  for (int j = 0; j < nvarY; j++){
    std::string name = "S" + std::to_string(j+1);
    out[name] = out_S[j];
  }
  out["COST"] = out_COST;
  out["CENS"] = out_CENS;
  out["DOM1"] = out_DOM1;   // MODIFICA: ora sarà un vettore di interi
  for (int j = 0; j < nvarX; j++){
    std::string name = "X" + std::to_string(j+1);
    out[name] = out_X[j];
  }
  
  if(verbose) {
    Rcout << "\nNumber of strata: " << out_STRATO.size() << "\n";
    int one_unit = 0;
    for (double val : out_N) {
      if(val == 1.0) one_unit++;
    }
    Rcout << "... of which with only one unit: " << one_unit << "\n";
  }
  
  return as<DataFrame>(out);
}
