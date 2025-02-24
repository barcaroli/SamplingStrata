// [[Rcpp::plugins(openmp)]]


#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
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

double stdev3(const std::vector<double>& Y, const std::vector<double>& W, double beta1, double beta2) {
  size_t n = Y.size();
  if(n < 2) return 0.0;
  double meanY = 0.0, meanW = 0.0;
  for (size_t i = 0; i < n; i++){
    meanY += Y[i];
    meanW += W[i];
  }
  meanY /= n;
  meanW /= n;
  double varY = 0.0, varW = 0.0, covYW = 0.0;
  for (size_t i = 0; i < n; i++){
    varY += (Y[i] - meanY) * (Y[i] - meanY);
    varW += (W[i] - meanW) * (W[i] - meanW);
    covYW += (Y[i] - meanY) * (W[i] - meanW);
  }
  varY /= (n - 1);
  varW /= (n - 1);
  covYW /= (n - 1);
  double value = beta1 * beta1 * varY + 2 * beta1 * beta2 * covYW + beta2 * beta2 * varW;
  return std::sqrt(value);
}

// Ottimizzazione di stdev4: iterazione solo per i < j e parallelizzazione (se il gruppo è grande)
double stdev4(const std::vector<double>& Y,
              const std::vector<double>& LON,
              const std::vector<double>& LAT,
              double var_eps,
              double range,
              double gamma) {
  size_t n = Y.size();
  if(n == 0) return 0.0;
  std::vector<double> var_ntimes(n);
  for (size_t i = 0; i < n; i++){
    var_ntimes[i] = var_eps * std::pow(Y[i], 2 * gamma);
  }
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
      double spatial_autocov = prod_std * std::exp(- d / (range + 1e-7));
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
  std::vector<std::string> DOM1;
  std::vector< std::vector<double> > X; // dimensione: nvarX
};

// [[Rcpp::export]]
DataFrame buildStrataDF(DataFrame dataset, Nullable<DataFrame> model_ = R_NilValue, 
                             bool progress = true, bool verbose = true) {
  CharacterVector colNames = dataset.names();
  int nrows = dataset.nrows();
  
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
  if (dataset.containsElementNamed("WEIGHT"))
    weight = dataset["WEIGHT"], weightProvided = true;
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
    if(model.nrows() != nvarY)
      stop("A model for each Y variable must be specified");
  }
  
  if(!dataset.containsElementNamed("DOMAINVALUE"))
    stop("DOMAINVALUE column missing in dataset");
  CharacterVector domainValues = dataset["DOMAINVALUE"];
  
  // Raccolta dei domini unici
  std::set<std::string> domainsSet;
  for (int i = 0; i < domainValues.size(); i++){
    std::string dom = as<std::string>(domainValues[i]);
    domainsSet.insert(dom);
  }
  std::vector<std::string> domains(domainsSet.begin(), domainsSet.end());
  
  // Estrae le colonne X
  std::vector< NumericVector > Xcols;
  for (int j = 1; j <= nvarX; j++){
    std::string colName = "X" + std::to_string(j);
    if(!dataset.containsElementNamed(colName.c_str()))
      stop(("Column " + colName + " missing in dataset").c_str());
    NumericVector col = dataset[colName.c_str()];
    Xcols.push_back(col);
  }
  
  // Prepara un vettore di risultati per ciascun dominio
  std::vector<DomainResult> domain_results(domains.size());
  
  // Parallelizza sul loop sui domini (ogni dominio viene processato in modo indipendente)
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int d = 0; d < (int)domains.size(); d++){
    DomainResult local;
    std::string dom = domains[d];
    // Indici delle righe per il dominio corrente
    std::vector<int> idx;
    for (int i = 0; i < nrows; i++){
      std::string dval = as<std::string>(domainValues[i]);
      if(dval == dom)
        idx.push_back(i);
    }
    // Costruisce lo "strato" concatenando i valori delle colonne X
    std::vector<std::string> strato(idx.size());
    for (size_t k = 0; k < idx.size(); k++){
      int i = idx[k];
      std::ostringstream oss;
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
      std::string strata_id = g.first;
      std::vector<int> group_indices = g.second;
      double sumW = 0.0;
      double sumCost = 0.0;
      bool hasCost = dataset.containsElementNamed("COST");
      NumericVector COSTcol;
      if(hasCost) COSTcol = dataset["COST"];
      
      for (int i : group_indices) {
        sumW += weight[i];
        if(hasCost)
          sumCost += weight[i] * COSTcol[i];
      }
      double N_val = sumW;
      double cost_val = hasCost ? (sumCost / sumW) : 1.0;
      
      std::vector<double> M_vals(nvarY, 0.0);
      std::vector<double> S_vals(nvarY, 0.0);
      
      // Calcola M e S per ciascuna variabile Y
      for (int y_idx = 1; y_idx <= nvarY; y_idx++){
        std::string colY = "Y" + std::to_string(y_idx);
        if(!dataset.containsElementNamed(colY.c_str()))
          stop(("Column " + colY + " missing").c_str());
        NumericVector Ycol = dataset[colY.c_str()];
        
        NumericVector Wcol;
        if(useModel){
          std::string colW = "W" + std::to_string(y_idx);
          if(dataset.containsElementNamed(colW.c_str()))
            Wcol = dataset[colW.c_str()];
          else
            Wcol = NumericVector(nrows, 0.0);
        }
        
        std::vector<double> Y_vals, w_vals, W_vals;
        if(useModel){
          CharacterVector model_type = model["type"];
          std::string mtype = as<std::string>(model_type[y_idx-1]);
          for (int i : group_indices) {
            double y_val = Ycol[i];
            if(!R_IsNA(y_val)) {
              Y_vals.push_back(y_val);
              w_vals.push_back(weight[i]);
              if(mtype == "spatial")
                W_vals.push_back(Wcol[i]);
            }
          }
        } else {
          for (int i : group_indices) {
            double y_val = Ycol[i];
            if(!R_IsNA(y_val)) {
              Y_vals.push_back(y_val);
              w_vals.push_back(weight[i]);
            }
          }
        }
        
        if(Y_vals.size() == 0){
          M_vals[y_idx-1] = 0.0;
          S_vals[y_idx-1] = 0.0;
          continue;
        }
        
        if(!useModel) {
          double sumYW = 0.0, sumw = 0.0;
          for (size_t j = 0; j < Y_vals.size(); j++){
            sumYW += Y_vals[j] * w_vals[j];
            sumw += w_vals[j];
          }
          double M_val = sumYW / sumw;
          M_vals[y_idx-1] = M_val;
          double s_val = weightProvided ? stdev1(Y_vals, w_vals) : stdev2(Y_vals, w_vals);
          S_vals[y_idx-1] = s_val;
        } else {
          CharacterVector model_type = model["type"];
          NumericVector model_beta = model["beta"];
          NumericVector model_sig2 = model["sig2"];
          std::string mtype = as<std::string>(model_type[y_idx-1]);
          double beta = model_beta[y_idx-1];
          double sig2 = model_sig2[y_idx-1];
          
          if(mtype == "linear") {
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += Y_vals[j] * beta * w_vals[j];
              sumw += w_vals[j];
            }
            double M_val = sumYW / sumw;
            M_vals[y_idx-1] = M_val;
            double s_raw = weightProvided ? stdev1(Y_vals, w_vals) : stdev2(Y_vals, w_vals);
            if(!model.containsElementNamed("gamma"))
              stop("gamma missing in model for linear type");
            NumericVector model_gamma = model["gamma"];
            double gamma_val = model_gamma[y_idx-1];
            double sum_pow = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sum_pow += std::pow(Y_vals[j], 2 * gamma_val);
            }
            double gammas = sum_pow / Y_vals.size();
            S_vals[y_idx-1] = std::sqrt(s_raw * s_raw * beta * beta + sig2 * gammas);
          } else if(mtype == "loglinear") {
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += std::pow(Y_vals[j], beta) * w_vals[j];
              sumw += w_vals[j];
            }
            double M_val = sumYW / sumw;
            M_vals[y_idx-1] = M_val;
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
            S_vals[y_idx-1] = std::sqrt(ph * ((std::exp(sig2) * (sumY2 / Y_vals.size())) - ph * (avgYbeta * avgYbeta)));
          } else if(mtype == "spatial") {
            if(W_vals.size() != Y_vals.size())
              stop("Mismatch in Y and W lengths for spatial model");
            if(!model.containsElementNamed("beta2"))
              stop("beta2 missing in model for spatial type");
            NumericVector model_beta2 = model["beta2"];
            double beta2 = model_beta2[y_idx-1];
            double sumYW = 0.0, sumw = 0.0;
            for (size_t j = 0; j < Y_vals.size(); j++){
              sumYW += (Y_vals[j] * beta + W_vals[j] * beta2) * w_vals[j];
              sumw += w_vals[j];
            }
            double M_val = sumYW / sumw;
            M_vals[y_idx-1] = M_val;
            double s1 = stdev3(Y_vals, W_vals, beta, beta2);
            if(!dataset.containsElementNamed("LON") || !dataset.containsElementNamed("LAT"))
              stop("Missing coordinates on sampling frame");
            NumericVector LONcol = dataset["LON"];
            NumericVector LATcol = dataset["LAT"];
            std::vector<double> lon_vals, lat_vals;
            for (int i : group_indices) {
              lon_vals.push_back(LONcol[i]);
              lat_vals.push_back(LATcol[i]);
            }
            if(!model.containsElementNamed("range") || !model.containsElementNamed("gamma") ||
               !model.containsElementNamed("fitting"))
               stop("Missing parameters in model for spatial type");
            NumericVector model_range = model["range"];
            NumericVector model_gamma = model["gamma"];
            NumericVector model_fitting = model["fitting"];
            double range_val = model_range[y_idx-1];
            double gamma_val = model_gamma[y_idx-1];
            double fitting = model_fitting[y_idx-1];
            double s2 = stdev4(Y_vals, lon_vals, lat_vals, sig2, range_val, gamma_val);
            S_vals[y_idx-1] = std::sqrt((s1 * s1 / fitting) + s2 * s2);
          } else {
            stop("Type of model misspecified");
          }
        }
      }
      
      // Salva i risultati locali per lo strato corrente
      local.STRATO.push_back(strata_id);
      local.N.push_back(N_val);
      for (int j = 0; j < nvarY; j++){
        local.M[j].push_back(M_vals[j]);
        local.S[j].push_back(S_vals[j]);
      }
      local.COST.push_back(cost_val);
      local.CENS.push_back(0.0);
      local.DOM1.push_back(dom);
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
    if(progress) {
#ifdef _OPENMP
#pragma omp critical
#endif
      Rcout << "Processed domain: " << dom << "\n";
    }
  } // Fine ciclo sui domini
  
  // Unione dei risultati di tutti i domini
  std::vector<std::string> out_STRATO;
  std::vector<double> out_N;
  std::vector< std::vector<double> > out_M(nvarY);
  std::vector< std::vector<double> > out_S(nvarY);
  std::vector<double> out_COST;
  std::vector<double> out_CENS;
  std::vector<std::string> out_DOM1;
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
      out_DOM1.push_back(local.DOM1[i]);
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
  out["DOM1"] = out_DOM1;
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
