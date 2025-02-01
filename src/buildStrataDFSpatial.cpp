#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <sstream>

// Per comodità
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Funzione di utilità: trasforma un std::string in maiuscolo
// -----------------------------------------------------------------------------
static std::string strToUpper(const std::string &s) {
  std::string out = s;
  std::transform(out.begin(), out.end(), out.begin(), ::toupper);
  return out;
}

// -----------------------------------------------------------------------------
// Funzione di utilità ottimizzata: calcola la "stdev spaziale" su un subset
// di righe definito da 'indices', lavorando direttamente su LON, LAT, Y, VAR.
// Usa un loop ridotto (r1 < r2) e raddoppia le somme, per ridurre i calcoli.
// -----------------------------------------------------------------------------
static double stdevOptim(const std::vector<int> &indices,
                         const NumericVector &lon,
                         const NumericVector &lat,
                         const NumericVector &y,
                         const NumericVector &var,
                         double fittingVal,
                         double rangeVal,
                         double kappa)
{
  int n = indices.size();
  if (n <= 1) {
    // Se c'è 1 o 0 righe, stdev = 0
    return 0.0;
  }
  
  long double sum_zz = 0.0;
  long double sum_var_minus_cov = 0.0;
  
  // Pre-calcolo per l'esponenziale
  double inv_r = kappa / rangeVal;  // => exp(- dist * inv_r)
  
  for (int i = 0; i < n; i++) {
    int r1 = indices[i];
    for (int j = i + 1; j < n; j++) {
      int r2 = indices[j];
      
      double dlon = lon[r1] - lon[r2];
      double dlat = lat[r1] - lat[r2];
      double dist = std::sqrt(dlon * dlon + dlat * dlat);
      
      double diffy = (y[r1] - y[r2]);
      double z_z   = diffy * diffy;
      sum_zz += z_z;
      
      double scv = var[r1] + var[r2];
      double pcs = 0.0;
      if (var[r1] >= 0.0 && var[r2] >= 0.0) {
        pcs = std::sqrt(var[r1] * var[r2]);
      }
      double spatial_cov = pcs * std::exp(- dist * inv_r);
      
      sum_var_minus_cov += (scv - 2.0 * spatial_cov);
    }
  }
  
  // Abbiamo calcolato la somma solo per (r1 < r2), raddoppiamo
  sum_zz *= 2.0;
  sum_var_minus_cov *= 2.0;
  
  // denom = 2 * n^2 (come nel tuo codice originale)
  long double denom = 2.0L * (long double)n * (long double)n;
  double sd1 = std::sqrt((long double)sum_zz / denom);
  double sd2 = std::sqrt((long double)sum_var_minus_cov / denom);
  
  double var_strato = sd1*sd1 / fittingVal + sd2*sd2;
  if (var_strato < 0.0) {
    var_strato = 0.0;
  }
  return std::sqrt(var_strato);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame buildStrataDFSpatial(
    DataFrame dataset,
    NumericVector fitting = NumericVector::create(1.0),
    NumericVector range   = NumericVector::create(0.0),
    double kappa = 3.0,
    bool progress = false,
    bool verbose  = false)
{
  // 1) Converte i nomi delle colonne in maiuscolo
  CharacterVector cn = dataset.names();
  for (int i = 0; i < cn.size(); i++) {
    cn[i] = strToUpper(Rcpp::as<std::string>(cn[i]));
  }
  dataset.names() = cn;
  
  // 2) Verifica che esista DOMAINVALUE
  if (!dataset.containsElementNamed("DOMAINVALUE")) {
    stop("Manca la colonna DOMAINVALUE nel dataset!");
  }
  NumericVector domainValue = dataset["DOMAINVALUE"];
  
  // Controlla LON, LAT
  if (!dataset.containsElementNamed("LON") || !dataset.containsElementNamed("LAT")) {
    stop("Mancano le colonne LON o LAT nel dataset!");
  }
  NumericVector lon = dataset["LON"];
  NumericVector lat = dataset["LAT"];
  
  // 3) Identifica le posizioni di X, Y, VAR
  std::vector<int> XcolsIndex, YcolsIndex, VARcolsIndex;
  for (int i = 0; i < cn.size(); i++) {
    std::string name = Rcpp::as<std::string>(cn[i]);
    if (name.rfind("X", 0) == 0) {
      XcolsIndex.push_back(i);
    } else if (name.rfind("Y", 0) == 0) {
      YcolsIndex.push_back(i);
    } else if (name.rfind("VAR", 0) == 0) {
      VARcolsIndex.push_back(i);
    }
  }
  
  int nvarX = XcolsIndex.size();
  int nvarY = YcolsIndex.size();
  
  if (verbose) {
    Rcout << "[buildStrataDFSpatial2] nvarX=" << nvarX
          << ", nvarY=" << nvarY << "\n";
    Rcout << "Computations are being done on population data...\n";
  }
  
  // 4) Prepara i vettori Y e VAR (colonne) in modo da accedervi per indice
  std::vector<NumericVector> Ycols(nvarY), VARcols(nvarY);
  for (int iY = 0; iY < nvarY; iY++) {
    // Converte la entry di cn in un Rcpp::String, poi lo usa come indice di dataset
    Rcpp::String yName   = cn[ YcolsIndex[iY] ];
    Rcpp::String varName = cn[ VARcolsIndex[iY] ];
    
    Ycols[iY]   = dataset[yName];
    VARcols[iY] = dataset[varName];
  }
  
  // 5) Prepara i vettori X
  std::vector<NumericVector> XcolsVec(nvarX);
  for (int ix = 0; ix < nvarX; ix++) {
    Rcpp::String xName = cn[ XcolsIndex[ix] ];
    XcolsVec[ix] = dataset[xName];
  }
  
  // 6) Costruisce la mappa: domainValue -> elenco di indici riga
  std::map<double, std::vector<int>> domain2rows;
  int nrows = dataset.nrows();
  for (int r = 0; r < nrows; r++) {
    double dom = domainValue[r];
    domain2rows[dom].push_back(r);
  }
  int numdom = domain2rows.size();
  
  // 7) Prepara vettori di output
  std::vector<std::string> out_STRATO;
  std::vector<double> out_N;
  std::vector<double> out_DOM1;
  std::vector<double> out_COST;
  std::vector<double> out_CENS;
  
  // M e S per ogni Y
  std::vector< std::vector<double> > out_M(nvarY), out_S(nvarY);
  
  // X1..XnvarX
  std::vector< std::vector<double> > out_X(nvarX);
  
  // 8) Loop su ciascun dominio
  int domCount = 0;
  for (auto &kv : domain2rows) {
    domCount++;
    double domVal = kv.first;
    const std::vector<int> &rowsDom = kv.second;  // Indici di riga
    
    if (progress) {
      Rcout << " [Domain " << domVal << " : " << domCount << "/" << numdom << "]\n";
    }
    
    // Per questo dominio, raggruppa per strato (concatenando X1, X2, ecc.)
    std::unordered_map<std::string, std::vector<int>> strato2rows;
    for (int r : rowsDom) {
      std::ostringstream oss;
      for (int ix = 0; ix < nvarX; ix++) {
        if (ix > 0) oss << "*";
        oss << XcolsVec[ix][r];
      }
      std::string sKey = oss.str();
      strato2rows[sKey].push_back(r);
    }
    
    // 9) Per ogni strato, calcola i valori e aggiorna i vettori di output
    for (auto &st : strato2rows) {
      std::string sKey = st.first;
      const std::vector<int> &srows = st.second;
      int sN = srows.size();
      
      // Calcolo M_i, S_i
      std::vector<double> M_i(nvarY, 0.0), S_i(nvarY, 0.0);
      
      for (int iY = 0; iY < nvarY; iY++) {
        // mean
        double sumY = 0.0;
        int countY = 0;
        for (int idx : srows) {
          double val = Ycols[iY][idx];
          if (!NumericVector::is_na(val)) {
            sumY += val;
            countY++;
          }
        }
        double meanY = (countY > 0) ? (sumY / (double)countY) : 0.0;
        M_i[iY] = meanY;
        
        // stdev spaziale
        double fittingVal = fitting[iY];  // Presupponendo lunghezza sufficiente
        double rangeVal   = range[iY];
        double sdev = stdevOptim(srows, lon, lat, Ycols[iY], VARcols[iY],
                                 fittingVal, rangeVal, kappa);
        S_i[iY] = sdev;
      }
      
      // Decompone sKey in X1..XnvarX
      std::vector<double> splittedX(nvarX, 0.0);
      {
        std::stringstream ss(sKey);
        std::string token;
        int ix = 0;
        while (std::getline(ss, token, '*')) {
          splittedX[ix++] = std::stod(token);
        }
      }
      
      // Salva i risultati per lo strato
      out_STRATO.push_back(sKey);
      out_N.push_back((double)sN);
      out_DOM1.push_back(domVal);
      out_COST.push_back(1.0); // come nell'originale
      out_CENS.push_back(0.0); // come nell'originale
      
      // M_i, S_i
      for (int iY = 0; iY < nvarY; iY++) {
        out_M[iY].push_back(M_i[iY]);
        out_S[iY].push_back(S_i[iY]);
      }
      // X1..XnvarX
      for (int ix = 0; ix < nvarX; ix++) {
        out_X[ix].push_back(splittedX[ix]);
      }
    }
  }
  
  // 10) Assemble finale in DataFrame
  int totalSize = out_STRATO.size();
  CharacterVector STRATO_(totalSize);
  NumericVector   N_(totalSize), DOM1_(totalSize), COST_(totalSize), CENS_(totalSize);
  
  for (int i = 0; i < totalSize; i++) {
    STRATO_[i] = out_STRATO[i];
    N_[i]      = out_N[i];
    DOM1_[i]   = out_DOM1[i];
    COST_[i]   = out_COST[i];
    CENS_[i]   = out_CENS[i];
  }
  
  // M_i, S_i
  std::vector<NumericVector> Mcols(nvarY), Scols(nvarY);
  for (int iY = 0; iY < nvarY; iY++) {
    Mcols[iY] = NumericVector(totalSize);
    Scols[iY] = NumericVector(totalSize);
    for (int j = 0; j < totalSize; j++) {
      Mcols[iY][j] = out_M[iY][j];
      Scols[iY][j] = out_S[iY][j];
    }
  }
  
  // X1..XnvarX
  std::vector<NumericVector> XcolsOut(nvarX);
  for (int ix = 0; ix < nvarX; ix++) {
    XcolsOut[ix] = NumericVector(totalSize);
    for (int j = 0; j < totalSize; j++) {
      XcolsOut[ix][j] = out_X[ix][j];
    }
  }
  
  // Costruisci la lista di colonne
  List outList;
  outList["STRATO"] = STRATO_;
  outList["N"]      = N_;
  
  for (int iY = 0; iY < nvarY; iY++) {
    std::string mName = "M" + std::to_string(iY+1);
    std::string sName = "S" + std::to_string(iY+1);
    outList[mName] = Mcols[iY];
    outList[sName] = Scols[iY];
  }
  
  outList["COST"] = COST_;
  outList["CENS"] = CENS_;
  outList["DOM1"] = DOM1_;
  
  for (int ix = 0; ix < nvarX; ix++) {
    std::string xName = "X" + std::to_string(ix+1);
    outList[xName] = XcolsOut[ix];
  }
  
  DataFrame stratatot = DataFrame(outList);
  
  if (verbose) {
    Rcout << "\nNumber of strata: " << totalSize << std::endl;
    int countN1 = 0;
    for (double nn : out_N) {
      if (nn == 1.0) countN1++;
    }
    Rcout << "...of which with only one unit: " << countN1 << std::endl;
  }
  
  return stratatot;
}
