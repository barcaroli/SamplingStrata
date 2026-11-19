#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

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
// Optimised utility: "spatial" standard deviation of the stratum made of the
// rows in 'indices', computed directly on LON, LAT, Y, VAR.
// Loop on pairs r1 < r2 only, sums doubled.
//
// Model-expected squared difference between two distinct units i, j
// (de Gruijter et al. 2015, Eq. 10; de Gruijter et al. 2016, Eq. 15):
//
//   D2_ij = (y_i - y_j)^2 / fitting + s2_i + s2_j - 2 Cov(e_i, e_j)
//   Cov(e_i, e_j) = (1 - nugget) * s_i * s_j * rho(d_ij)
//
// where rho is the correlation function of the variogram model fitted on the
// prediction errors (r = range):
//   exponential  rho(d) = exp(-kappa * d / r)
//   spherical    rho(d) = 1 - 1.5 (d / r) + 0.5 (d / r)^3 for d < r, 0 for d >= r
//   gaussian     rho(d) = exp(-kappa * (d / r)^2)
// kappa is not used by the spherical model (r is the distance at which the
// variogram reaches the sill).
//
// where s2 = VAR is the prediction-error variance of each unit and nugget is
// the RELATIVE nugget c0 / (c0 + c1) of the variogram of the prediction
// errors (0 = no nugget). The stratum variance is
//   S2_h = sum_{i != j} D2_ij / (2 n^2),
// i.e. the population variance (divisor n) of the stratum, consistently with
// buildStrataDF() for frame data.
//
// fitting: when the VAR columns contain the full prediction-error variances
// fitting MUST be 1 (additive correction): the levelling of the predictions
// is already compensated by s2_i + s2_j, and dividing also by R^2 counts the
// error variance twice. fitting = R^2 (< 1) is meant only for the rescaling
// correction, i.e. with all VAR = 0.
//
// range <= 0 or nugget >= 1: no spatial correlation (Cov = 0 for i != j).
// -----------------------------------------------------------------------------

// Variogram models of the prediction errors
enum { VGM_EXP = 0, VGM_SPH = 1, VGM_GAU = 2 };

// Code of a variogram model: full name or gstat abbreviation, case-insensitive
static int vgmModelCode(const std::string &s) {
  std::string m = s;
  std::transform(m.begin(), m.end(), m.begin(), ::tolower);
  if (m == "exponential" || m == "exp") return VGM_EXP;
  if (m == "spherical"   || m == "sph") return VGM_SPH;
  if (m == "gaussian"    || m == "gau") return VGM_GAU;
  stop("buildStrataDFSpatial: 'vgm_model' must be \"exponential\" (\"Exp\"), "
       "\"spherical\" (\"Sph\") or \"gaussian\" (\"Gau\"); found \"" + s + "\"");
  return VGM_EXP;  // not reached
}

// Correlation of the prediction errors of two units at distance d (r > 0);
// kr = kappa / r is precomputed for the exponential model
static inline double spatialCorr(double d, double r, double kappa, double kr, int model) {
  switch (model) {
  case VGM_SPH: {
    if (d >= r) return 0.0;
    double h = d / r;
    return 1.0 - 1.5 * h + 0.5 * h * h * h;
  }
  case VGM_GAU: {
    double h = d / r;
    return std::exp(- kappa * h * h);
  }
  default:
    return std::exp(- d * kr);
  }
}

static double stdevOptim(const std::vector<int> &indices,
                         const NumericVector &lon,
                         const NumericVector &lat,
                         const NumericVector &y,
                         const NumericVector &var,
                         double fittingVal,
                         double rangeVal,
                         double kappa,
                         double nuggetVal,
                         int vgmModel)
{
  int n = indices.size();
  if (n <= 1) {
    // Se c'è 1 o 0 righe, stdev = 0
    return 0.0;
  }
  
  long double sum_zz = 0.0;
  long double sum_var_minus_cov = 0.0;
  
  // Spatial correlation switched off with range <= 0 or nugget >= 1
  bool spatial = (rangeVal > 0.0) && (nuggetVal < 1.0);
  double inv_r = spatial ? kappa / rangeVal : 0.0;  // exponential: exp(- dist * inv_r)
  double w = 1.0 - nuggetVal;                       // structured share of the error variance
  
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
      // covariance of the prediction errors of two DISTINCT units: the
      // nugget part is uncorrelated even at distance 0 (pairs i = j are not
      // in the loop and contribute 0 to the stratum variance)
      double spatial_cov = spatial ?
        w * pcs * spatialCorr(dist, rangeVal, kappa, inv_r, vgmModel) : 0.0;
      
      sum_var_minus_cov += (scv - 2.0 * spatial_cov);
    }
  }
  
  // Abbiamo calcolato la somma solo per (r1 < r2), raddoppiamo
  sum_zz *= 2.0;
  sum_var_minus_cov *= 2.0;
  
  // denom = 2 * n^2: population variance of the stratum
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
//' @param fitting Default is 1, which must be used when the frame contains the
//'   prediction-error variances (VAR columns): the additive correction
//'   (de Gruijter et al. 2016) already accounts for the levelling of the predictions.
//'   A value < 1 (the coefficient of determination R2 of the prediction model)
//'   divides the squared differences of the predictions, and is accepted only with all
//'   VAR = 0 (rescaling correction): with positive variances the error variance would be
//'   counted twice, and the function stops. If only R2 is known, the additive correction
//'   can be used with fitting = 1 and a constant variance var(predictions) * (1 - R2) / R2.
//' @param range Range of the correlation function of the prediction errors
//'   (in the units of the coordinates), as defined by the variogram model given in
//'   vgm_model. 0 means no spatial correlation.
//' @param kappa Factor used in evaluating spatial autocorrelation with the
//'   exponential and gaussian models: the correlation between the prediction errors
//'   of two units at distance d is (1 - nugget) * exp(-kappa * d / range) (exponential)
//'   or (1 - nugget) * exp(-kappa * (d / range)^2) (gaussian). Use kappa = 1 when range
//'   is the range parameter of the variogram (as returned by gstat), kappa = 3 when
//'   range is the practical range. It is not used by the spherical model.
//'   Any finite value > 0 is accepted; 0 or negative values give an error. Default is 3.
//' @param progress If set to TRUE, a progress bar is visualised during the execution. Default is FALSE.
//' @param verbose If set to TRUE, information is given about the number of strata generated. Default is FALSE.
//' @param nugget Relative nugget c0/(c0+c1) of the variogram of the prediction errors,
//'   in [0, 1]: the share of the error variance that is not spatially correlated.
//'   Default is 0 (no nugget).
//' @param vgm_model Variogram model fitted on the prediction errors, which defines
//'   the correlation function rho(d) used in the covariance
//'   (1 - nugget) * s_i * s_j * rho(d_ij):
//'   "exponential" (default): rho(d) = exp(-kappa * d / range);
//'   "spherical": rho(d) = 1 - 1.5 (d / range) + 0.5 (d / range)^3 for d < range,
//'   0 for d >= range (range is the distance at which the variogram reaches the sill);
//'   "gaussian": rho(d) = exp(-kappa * (d / range)^2).
//'   The gstat abbreviations "Exp", "Sph" and "Gau" are also accepted (case-insensitive).
//'   A single value, or one value for each target variable Y.
//'
//' @return
//' A dataframe containing strata
//'
//' @author
//' Giulio Barcaroli
//'
//' @examples
//' # synthetic spatial frame: one domain, one target variable with its
//' # prediction variance and the coordinates of each unit
//' set.seed(1234)
//' frame <- data.frame(id = 1:100, X1 = runif(100), Y1 = runif(100, 10, 20),
//'                     var1 = 0.5, lon = runif(100), lat = runif(100),
//'                     domainvalue = 1)
//' frame$X1 <- var.bin(frame$X1, bins = 3)   # stratification variable
//' strata <- buildStrataDFSpatial(dataset = frame, fitting = 1, range = 0.5,
//'                                kappa = 3, nugget = 0.3)
//' strata
//' # same frame, spherical variogram of the prediction errors
//' strata_sph <- buildStrataDFSpatial(dataset = frame, fitting = 1, range = 0.5,
//'                                    nugget = 0.3, vgm_model = "spherical")
//' strata_sph
//'
//' @keywords survey
//' @export
// [[Rcpp::export]]
DataFrame buildStrataDFSpatial(
    DataFrame dataset,
    NumericVector fitting = NumericVector::create(1.0),
    NumericVector range   = NumericVector::create(0.0),
    double kappa = 3.0,
    bool progress = false,
    bool verbose  = false,
    NumericVector nugget = NumericVector::create(0.0),
    CharacterVector vgm_model = CharacterVector::create("exponential"))
{
  // 'fitting' appears below (in stdevOptim) as sd1*sd1/fittingVal in the
  // within-stratum variance (see the note on the additive and rescaling
  // corrections above stdevOptim):
  // fitting = 0 is a division by zero, and any fitting close to 0 correctly
  // forces a very large required sample size (a model that explains nothing
  // needs near-census sampling to meet any precision target). Reject it here
  // instead of letting Inf/NaN propagate silently into bethel_cpp and every
  // genetic-algorithm evaluation, which used to surface much later as an
  // opaque "'ylim' values must be finite" crash inside plot.window().
  for (int i = 0; i < fitting.size(); i++) {
    if (fitting[i] == 0.0) {
      stop("buildStrataDFSpatial: 'fitting' cannot be 0. If your intent is to "
           "disable the spatial-correlation adjustment, use fitting = 1 "
           "together with range = 0 instead of fitting = 0.");
    }
  }

  // 1) Converte i nomi delle colonne in maiuscolo
  // Upper-case names on a COPY of the names and of the column list: modifying
  // dataset.names() in place changed the names of the caller's data frame.
  // The columns are not copied (shallow list of pointers).
  int nrows = dataset.nrows();
  CharacterVector cn = clone(as<CharacterVector>(dataset.names()));
  for (int i = 0; i < cn.size(); i++) {
    cn[i] = strToUpper(Rcpp::as<std::string>(cn[i]));
  }
  List ds(dataset.size());
  for (int i = 0; i < dataset.size(); i++) ds[i] = dataset[i];
  ds.attr("names") = cn;
  
  // 2) Verifica che esista DOMAINVALUE
  if (!ds.containsElementNamed("DOMAINVALUE")) {
    stop("Manca la colonna DOMAINVALUE nel dataset!");
  }
  NumericVector domainValue = ds["DOMAINVALUE"];
  
  // Controlla LON, LAT
  if (!ds.containsElementNamed("LON") || !ds.containsElementNamed("LAT")) {
    stop("Mancano le colonne LON o LAT nel dataset!");
  }
  NumericVector lon = ds["LON"];
  NumericVector lat = ds["LAT"];
  
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
  // one VAR column per Y (otherwise VARcolsIndex is read out of bounds)
  if ((int)VARcolsIndex.size() != nvarY)
    stop("buildStrataDFSpatial: the number of VAR columns must equal the number of Y columns");
  // one value of fitting and range per Y: a single value is recycled; any other
  // length would be read out of bounds (garbage values, NaN/Inf)
  if (fitting.size() != 1 && fitting.size() != nvarY)
    stop("buildStrataDFSpatial: 'fitting' must have length 1 or equal to the number of Y variables");
  if (range.size() != 1 && range.size() != nvarY)
    stop("buildStrataDFSpatial: 'range' must have length 1 or equal to the number of Y variables");
  if (nugget.size() != 1 && nugget.size() != nvarY)
    stop("buildStrataDFSpatial: 'nugget' must have length 1 or equal to the number of Y variables");
  if (vgm_model.size() != 1 && vgm_model.size() != nvarY)
    stop("buildStrataDFSpatial: 'vgm_model' must have length 1 or equal to the number of Y variables");
  NumericVector fit_y(nvarY), range_y(nvarY), nugget_y(nvarY);
  std::vector<int> vgm_y(nvarY);
  for (int iY = 0; iY < nvarY; iY++) {
    fit_y[iY]    = (fitting.size() == 1) ? fitting[0] : fitting[iY];
    range_y[iY]  = (range.size() == 1)   ? range[0]   : range[iY];
    nugget_y[iY] = (nugget.size() == 1)  ? nugget[0]  : nugget[iY];
    if (NumericVector::is_na(nugget_y[iY]) || nugget_y[iY] < 0.0 || nugget_y[iY] > 1.0)
      stop("buildStrataDFSpatial: 'nugget' must be in [0, 1] (relative nugget c0/(c0+c1))");
    String vm = (vgm_model.size() == 1) ? vgm_model[0] : vgm_model[iY];
    if (vm == NA_STRING)
      stop("buildStrataDFSpatial: 'vgm_model' cannot be NA");
    vgm_y[iY] = vgmModelCode(std::string(vm.get_cstring()));
  }
  
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
    
    Ycols[iY]   = ds[yName];
    VARcols[iY] = ds[varName];
  }
  
  // kappa must be finite and > 0 (no warning here: the function is called at
  // every evaluation of the genetic algorithm; the usual values 1 and 3 are
  // checked by optimizeStrataSpatial and KmeansSolutionSpatial)
  if (!R_finite(kappa) || kappa <= 0.0)
    stop("buildStrataDFSpatial: 'kappa' must be a finite value > 0");

  // fitting < 1 together with positive prediction-error variances would count
  // the error variance twice: rejected
  for (int iY = 0; iY < nvarY; iY++) {
    if (fit_y[iY] < 1.0) {
      for (int r = 0; r < nrows; r++) {
        double v = VARcols[iY][r];
        if (!NumericVector::is_na(v) && v > 0.0)
          stop("buildStrataDFSpatial: fitting < 1 cannot be used together with positive "
               "prediction-error variances (the error variance would be counted twice). "
               "Either omit fitting (default 1) and keep the variances, or, if only the R2 "
               "of the prediction model is known, use fitting = 1 with a constant variance "
               "equal to var(predictions) * (1 - R2) / R2.");
      }
    }
  }

  // 5) Prepara i vettori X
  std::vector<NumericVector> XcolsVec(nvarX);
  for (int ix = 0; ix < nvarX; ix++) {
    Rcpp::String xName = cn[ XcolsIndex[ix] ];
    XcolsVec[ix] = ds[xName];
  }
  
  // 6) Costruisce la mappa: domainValue -> elenco di indici riga
  std::map<double, std::vector<int>> domain2rows;
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
      // 17 significant digits: exact key (default precision 6 merged distinct values)
      oss << std::setprecision(17);
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
        double fittingVal = fit_y[iY];
        double rangeVal   = range_y[iY];
        double nuggetVal  = nugget_y[iY];
        double sdev = stdevOptim(srows, lon, lat, Ycols[iY], VARcols[iY],
                                 fittingVal, rangeVal, kappa, nuggetVal, vgm_y[iY]);
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
      // stratum label with 15 significant digits (as R does)
      {
        std::ostringstream lab;
        lab << std::setprecision(15);
        for (int ix = 0; ix < nvarX; ix++) {
          if (ix > 0) lab << "*";
          lab << splittedX[ix];
        }
        out_STRATO.push_back(lab.str());
      }
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
