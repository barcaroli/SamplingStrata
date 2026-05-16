//-------------------------------------------------------
// C++ script implementing Bethel algorithm (optimized)
// Author: Giulio Barcaroli (original), improvements by ...
// Comments in English
//-------------------------------------------------------

#include <Rcpp.h>
#include <iostream>  // for possible debugging (use Rcpp::Rcout in production)
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//-------------------------------------------------------
// 1) select_variables
//    Extracts columns with names prefix + (1..nvar) from a DataFrame
//-------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix select_variables(DataFrame dati, 
                               std::string prefix, 
                               int nvar) {
  StringVector cols(nvar);
  for (int i = 0; i < nvar; i++) {
    // Build the column name: prefix + (i+1)
    cols[i] = prefix + std::to_string(i + 1);
  }
  DataFrame subset = dati[cols];
  // Convert the subset to a NumericMatrix
  NumericMatrix mat = internal::convert_using_rfunction(subset, "as.matrix");  
  return mat;
}

//-------------------------------------------------------
// 2) disjoint
//    Creates a disjoint (dummy) matrix based on 'dom' values:
//    If dom[j] = k, then disj(j, k-1) = 1
//-------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix disjoint(const NumericVector& dom) { 
  // We assume 'dom' contains integer-like values >= 1
  int nstrat = dom.size();
  // maxdom is the maximum domain value
  int maxdom = (int) Rcpp::max(dom); 
  IntegerMatrix disj(nstrat, maxdom);
  
  // Fill the matrix: if dom[j] = k, we set disj(j, k-1) = 1
  for (int j = 0; j < nstrat; j++) {
    int val = (int)dom[j] - 1;
    if (val >= 0 && val < maxdom) {
      disj(j, val) = 1;
    }
  }
  return disj;
}

//-------------------------------------------------------
// 3) m_s
//    Multiplies each row of 'mat' by the corresponding
//    indicator in 'disj' and combines them horizontally.
//    If disj(i,k) = 1, columns for that domain are copied; otherwise 0
//-------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix m_s(const IntegerMatrix& disj,
                  const NumericMatrix& mat) {
  int nc = disj.ncol();      // number of domains
  int nstrat = disj.nrow();  // number of strata
  int nvar = mat.ncol();     // number of variables
  
  // The output matrix will have size (nstrat) x (nc * nvar)
  NumericMatrix out(nstrat, nc * nvar);
  
  // For each domain (k), fill the block of columns [k*nvar : (k+1)*nvar - 1]
  for (int k = 0; k < nc; k++) {
    int offset = k * nvar;
    for (int i = 0; i < nstrat; i++) {
      double multiplier = (disj(i, k) == 1) ? 1.0 : 0.0;
      for (int j = 0; j < nvar; j++) {
        out(i, offset + j) = mat(i, j) * multiplier;
      }
    }
  }
  return out;
}

//-------------------------------------------------------
// 4) rowSums_Rcpp
//    Calculates row sums of a NumericMatrix
//-------------------------------------------------------
// [[Rcpp::export]]
NumericVector rowSums_Rcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int i = 0; i < nr; i++) {
    double sum = 0.0;
    for (int j = 0; j < nc; j++) {
      sum += x(i, j);
    }
    ans[i] = sum;
  }
  return ans;
}

//-------------------------------------------------------
// 5) colSums_Rcpp
//    Calculates column sums of a NumericMatrix
//-------------------------------------------------------
// [[Rcpp::export]]
NumericVector colSums_Rcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += x(i, j);
    }
    ans[j] = sum;
  }
  return ans;
}

//-------------------------------------------------------
// 6) cv_Rcpp
//    Example of building a 1 x (nvar*ndom) matrix of CVs.
//    If you really want the same CV repeated for each domain,
//    this approach replicates the same row n times horizontally.
//-------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix cv_Rcpp(DataFrame errors, 
                      int ndom, 
                      int nvar) {
  // Select columns "CV1", "CV2", ..., "CVnvar" from 'errors'
  NumericMatrix cvx = select_variables(errors, "CV", nvar);
  // cvx presumably has dimension 1 x nvar (depending on 'errors')
  
  // We create a 1 x (nvar * ndom) output
  NumericMatrix out(1, nvar * ndom);
  
  // Replicate cvx horizontally ndom times
  for (int d = 0; d < ndom; d++) {
    int offset = d * nvar;
    for (int j = 0; j < nvar; j++) {
      out(0, offset + j) = cvx(0, j);
    }
  }
  return out;
}

//-------------------------------------------------------
// 7) crea_a
//    Constructs the matrix 'a' used by the Chromy algorithm
//    based on m, s, nocens, N, cv, etc.
//-------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix crea_a(NumericMatrix& m,
                     NumericMatrix& s,
                     NumericVector& nocens,
                     NumericVector& N,
                     NumericVector& cv,
                     double& epsilon) {
  int nr = m.nrow(), nc = m.ncol();
  NumericMatrix numA(nr, nc);
  
  // numA(i,j) = N(i)^2 * s(i,j)^2 * nocens(i)
  for (int i = 0; i < nr; i++) {
    double Ni = N[i];
    double Ni2 = Ni * Ni;
    double noCi = nocens[i];
    for (int j = 0; j < nc; j++) {
      double sij = s(i, j);
      numA(i, j) = Ni2 * sij * sij * noCi;
    }
  }
  
  // Build y = N(i)*m(i,j)
  NumericMatrix y(nr, nc);
  for (int i = 0; i < nr; i++) {
    double Ni = N[i];
    for (int j = 0; j < nc; j++) {
      y(i, j) = Ni * m(i, j);
    }
  }
  
  // We'll transpose and multiply by cv
  NumericMatrix yT(nc, nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      yT(j, i) = y(i, j);
    }
  }
  // Multiply each row j by cv[j]
  for (int j = 0; j < nc; j++) {
    double cvj = cv[j];
    for (int i = 0; i < nr; i++) {
      yT(j, i) *= cvj;
    }
  }
  // Transpose back to y
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      y(i, j) = yT(j, i);
    }
  }
  
  // denA1 = (colSums of y)^2
  NumericVector denA1 = colSums_Rcpp(y);
  for (int j = 0; j < nc; j++) {
    denA1[j] = denA1[j] * denA1[j];
  }
  
  // denA2 = colSums( N(i)*s(i,j)^2 * nocens(i) )
  NumericMatrix w(nr, nc);
  for (int i = 0; i < nr; i++) {
    double Ni = N[i];
    double noCi = nocens[i];
    for (int j = 0; j < nc; j++) {
      double sij = s(i, j);
      w(i, j) = Ni * sij * sij * noCi;
    }
  }
  NumericVector denA2 = colSums_Rcpp(w);
  
  // denA = denA1 + denA2 + epsilon
  NumericVector denA(nc);
  for (int j = 0; j < nc; j++) {
    denA[j] = denA1[j] + denA2[j] + epsilon;
  }
  
  // a = transpose(numA) / denA => then transpose again
  NumericMatrix z(nc, nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      z(j, i) = numA(i, j);
    }
  }
  for (int j = 0; j < nc; j++) {
    double dA = denA[j];
    for (int i = 0; i < nr; i++) {
      z(j, i) /= dA;
    }
  }
  // Transpose z back into 'a'
  NumericMatrix a(nr, nc);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      a(i, j) = z(j, i);
    }
  }
  return a;
}

//-------------------------------------------------------
// 8) chromy_Rcpp
//    Implements the Chromy algorithm to find allocations
//-------------------------------------------------------
// [[Rcpp::export]]
NumericVector chromy_Rcpp(NumericMatrix a,
                          double alfatot, 
                          double diff, 
                          int iter, 
                          NumericVector alfa, 
                          NumericVector alfanext, 
                          NumericVector x,
                          NumericVector cost,
                          int nvar,
                          bool realAllocation
) {
  int maxiter = 200; 
  double epsilon = 1e-11;
  int nr = a.nrow();
  int nc = nvar;
  
  NumericVector n(nr);  // final output
  
  while (diff > epsilon && iter < maxiter) {
    
    // 1) den1 = sqrt( rowSums( transpose( transpose(a)*alfa ) ) )
    // We'll create a transposed version of a in 'b'
    NumericMatrix b(nc, nr);
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
        b(j, i) = a(i, j) * alfa[j];
      }
    }
    // c = transpose(b)
    NumericMatrix c(nr, nc);
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
        c(i, j) = b(j, i);
      }
    }
    NumericVector den1 = rowSums_Rcpp(c);
    for (int i = 0; i < nr; i++) {
      den1[i] = std::sqrt(den1[i]);
    }
    
    // 2) den2 = sum( sqrt( rowSums( t(t(a * cost) * alfa) ) ) )
    for (int i = 0; i < nr; i++) {
      double ci = cost[i];
      for (int j = 0; j < nc; j++) {
        c(i, j) = a(i, j) * ci;
      }
    }
    // b = transpose(c), multiply by alfa
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
        b(j, i) = c(i, j) * alfa[j];
      }
    }
    // c = transpose(b) again
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
        c(i, j) = b(j, i);
      }
    }
    NumericVector d = rowSums_Rcpp(c);
    for (int i = 0; i < nr; i++) {
      d[i] = std::sqrt(d[i]);
    }
    double den2 = 0.0;
    for (int i = 0; i < nr; i++) {
      den2 += d[i];
    }
    
    // x(i) = sqrt(cost(i)) / (den1(i)*den2 + epsilon)
    for (int i = 0; i < nr; i++) {
      x[i] = std::sqrt(cost[i]) / (den1[i]*den2 + epsilon);
    }
    
    // alfatot = sum( alfa[j] * (t(a) %*% x)[j]^2 )
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
        b(j, i) = a(i, j);
      }
    }
    NumericVector e(nc, 0.0);
    for (int j = 0; j < nc; j++) {
      double sum_ = 0.0;
      for (int i = 0; i < nr; i++) {
        sum_ += b(j, i)* x[i];
      }
      e[j] = sum_;
    }
    double new_alfatot = 0.0;
    for (int j = 0; j < nc; j++) {
      new_alfatot += alfa[j] * (e[j]*e[j]);
    }
    if (new_alfatot == 0) new_alfatot = epsilon;
    
    // alfanext = (alfa[j]*e[j]^2) / alfatot
    for (int j = 0; j < nc; j++) {
      alfanext[j] = (alfa[j]* e[j]*e[j]) / new_alfatot;
    }
    
    // diff = max( abs( alfanext - alfa ) )
    double maxdiff = 0.0;
    for (int j = 0; j < nc; j++) {
      double tmp = std::fabs(alfanext[j] - alfa[j]);
      if (tmp > maxdiff) maxdiff = tmp;
    }
    diff = maxdiff;
    alfatot = new_alfatot;
    
    // update alfa
    for (int j = 0; j < nc; j++) {
      alfa[j] = alfanext[j];
    }
    iter++;
  }
  
  // n = 1/x if realAllocation = true, or ceil(1/x) otherwise
  if (realAllocation) {
    for (int i = 0; i < nr; i++) {
      n[i] = 1.0 / x[i];
    }
  } else {
    for (int i = 0; i < nr; i++) {
      n[i] = std::ceil(1.0 / x[i]);
    }
  }
  return n;
}

//-------------------------------------------------------
// 9) check_n
//    Ensures that n(i) respects [minnumstrat, N(i)] boundaries
//-------------------------------------------------------
// [[Rcpp::export]]
NumericVector check_n(NumericVector n,
                      NumericVector N,
                      int minnumstrat) {
  int nstrat = n.size();
  NumericVector n1(nstrat);
  
  for (int i = 0; i < nstrat; i++) {
    double ni = n[i];
    double Ni = N[i];
    if (ni > Ni) {
      n1[i] = Ni;
    } else if (ni < minnumstrat) {
      n1[i] = (Ni >= minnumstrat) ? (double)minnumstrat : Ni;
    } else {
      n1[i] = ni;
    }
  }
  return n1;
}

//-------------------------------------------------------
// 10) bethel
//     Main function to compute Bethel allocation
//-------------------------------------------------------
// [[Rcpp::export]]
NumericVector bethel_cpp(DataFrame strata,
                          DataFrame errors,
                          int minnumstrat,
                          bool realAllocation = true) {
  
  double epsilon = 1e-11;
  int nstrat = strata.nrows();
  
  // Determine number of CV variables (subtract domain columns)
  StringVector errors_names = errors.names();
  int nvar = errors_names.size() - 2;
  bool isPresent = (std::find(errors_names.begin(), errors_names.end(), "domainvalue") != errors_names.end());
  if (!isPresent) {
    nvar = errors_names.size() - 1;
  }
  
  // Extract mean (M...) and standard deviation (S...) from 'strata'
  NumericMatrix med = select_variables(strata, "M", nvar);
  NumericMatrix esse = select_variables(strata, "S", nvar);
  
  // Other vectors from 'strata'
  NumericVector N = strata["N"];
  NumericVector cens = strata["CENS"];
  NumericVector nocens(nstrat);
  for (int i = 0; i < nstrat; i++) {
    nocens[i] = 1.0 - cens[i];
  }
  NumericVector cost = strata["COST"];
  
  // 'DOM1' and creation of disjoint matrix
  NumericVector dom = strata["DOM1"];
  int ndom = (int) Rcpp::max(dom);
  
  IntegerMatrix disj = disjoint(dom);
  NumericMatrix m = m_s(disj, med);
  NumericMatrix s = m_s(disj, esse);
  
  // Retrieve CV matrix
  NumericMatrix cv = cv_Rcpp(errors, ndom, nvar);
  nvar = cv.ncol();  // should remain the same
  
  // Build matrix a
  NumericMatrix a = crea_a(m, s, nocens, N, cv, epsilon);
  
  // Chromy parameters
  double alfatot = 0.0;
  double diff = 999.0;
  int iter = 0;
  NumericVector alfa(nvar);
  NumericVector alfanext(nvar);
  NumericVector x(nstrat, 0.1);
  
  // Initialize alfa
  for (int j = 0; j < nvar; j++) {
    alfa[j] = 1.0 / nvar;
  }
  
  // First Chromy call
  NumericVector n = chromy_Rcpp(a, alfatot, diff, iter, alfa, alfanext, x, cost, nvar, realAllocation);
  
  // Check if n(i) > N(i)
  int contx = 0;
  for (int i = 0; i < nstrat; i++) {
    if (n[i] > N[i]) {
      contx++;
      cens[i] = 1.0;
      nocens[i] = 0.0;
    }
  }
  n = check_n(n, N, minnumstrat);
  
  // Iterate if overshoot occurs
  int iter1 = 0;
  int maxiter1 = 25;
  while (contx > 0 && iter1 < maxiter1) {
    iter1++;
    a = crea_a(m, s, nocens, N, cv, epsilon);
    n = chromy_Rcpp(a, alfatot, diff, iter, alfa, alfanext, x, cost, nvar, realAllocation);
    
    contx = 0;
    for (int i = 0; i < nstrat; i++) {
      if (n[i] > N[i]) {
        contx++;
        cens[i] = 1.0;
        nocens[i] = 0.0;
        n[i] = N[i];
      }
    }
    n = check_n(n, N, minnumstrat);
  }
  
  // Final step: n <- n * nocens + N * cens
  for (int i = 0; i < nstrat; i++) {
    n[i] = n[i]*nocens[i] + N[i]*cens[i];
  }
  
  return n;
}
