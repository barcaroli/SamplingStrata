#' Single-domain minimal integer allocation (internal)
#'
#' Worker of [min_int_allocation()]: computes the allocation for one domain,
#' using the first row of `errors`.
#'
#' @inheritParams min_int_allocation
#' @return A list with `n_h`, `n`, `cv`, `n_ceiling`, `n_real` and `method`.
#' @noRd

.min_int_allocation_1dom <- function(strata, errors, minnumstrat = 2,
                                     divisor = c("N_h-1", "N_h"),
                                     compare_ceiling = TRUE) {
  divisor <- match.arg(divisor)
  names(strata) <- toupper(names(strata))
  S_cols <- grep("^S[0-9]+$", names(strata), value = TRUE)
  S_cols <- S_cols[order(as.integer(sub("^S", "", S_cols)))]
  k <- length(S_cols)
  M_cols <- paste0("M", seq_len(k))
  if (!all(M_cols %in% names(strata))) stop("columns M1..M", k, " not found")
  cv_t <- if (is.data.frame(errors)) {
    as.numeric(errors[1, paste0("CV", seq_len(k))])
  } else rep_len(as.numeric(errors), k)

  Nh <- strata$N
  W  <- Nh / sum(Nh)
  S2 <- sapply(S_cols, function(s) strata[[s]]^2)
  if (!is.matrix(S2)) S2 <- matrix(S2, nrow = length(Nh))
  if (divisor == "N_h-1") S2 <- S2 * ifelse(Nh > 1, Nh / (Nh - 1), 1)
  Ybar <- sapply(M_cols, function(m) sum(W * strata[[m]]))
  Vt   <- (cv_t * Ybar)^2                      # target variances
  A    <- W^2 * S2                             # H x k: V_j = sum A_hj (1/n_h - 1/N_h)
  Vfun <- function(n) colSums(A * (1 / n - 1 / Nh))
  tol  <- 1 + 1e-10

  ## 1-2. marginal allocation from the lower bounds
  lo <- pmin(minnumstrat, Nh)
  n  <- lo
  V  <- Vfun(n)
  while (any(V > Vt * tol)) {
    viol <- V > Vt * tol
    gain <- as.vector((A[, viol, drop = FALSE] %*% (1 / Vt[viol])) *
                        (1 / n - 1 / (n + 1)))
    gain[n >= Nh] <- -Inf
    if (all(!is.finite(gain))) stop("constraints cannot be met (all strata take-all)")
    h <- which.max(gain)
    n[h] <- n[h] + 1
    V <- V - A[h, ] * (1 / (n[h] - 1) - 1 / n[h])
  }

  ## 3. removal of redundant units (useful only with several targets)
  if (k > 1) repeat {
    cand <- which(n > lo)
    if (length(cand) == 0) break
    ok <- sapply(cand, function(h) {
      Vn <- V + A[h, ] * (1 / (n[h] - 1) - 1 / n[h])
      if (any(Vn > Vt * tol)) NA else max(Vn / Vt)
    })
    if (all(is.na(ok))) break
    h <- cand[which.min(ok)]
    n[h] <- n[h] - 1
    V <- Vfun(n)
  }

  out <- list(n_h = n, n = sum(n), cv = sqrt(Vfun(n)) / abs(Ybar),
              n_ceiling = NA, n_real = NA, method = "greedy")

  ## comparison with the Bethel allocation rounded up stratum by stratum
  if (compare_ceiling && requireNamespace("SamplingStrata", quietly = TRUE)) {
    s <- strata
    s[S_cols] <- lapply(S_cols, function(c) sqrt(S2[, match(c, S_cols)]))
    if (!"COST" %in% names(s)) s$COST <- 1
    if (!"CENS" %in% names(s)) s$CENS <- 0
    if (!"DOM1" %in% names(s)) s$DOM1 <- 1
    err <- data.frame(DOM = "DOM1", t(cv_t), domainvalue = 1)
    names(err)[1 + seq_len(k)] <- paste0("CV", seq_len(k))
    nb <- as.numeric(SamplingStrata::bethel(s, err, minnumstrat = minnumstrat,
                                            realAllocation = TRUE))
    nc <- pmin(Nh, pmax(lo, ceiling(nb - 1e-9)))
    out$n_real <- sum(nb)
    out$n_ceiling <- sum(nc)
    if (sum(nc) < out$n && all(Vfun(nc) <= Vt * tol)) {
      out$n_h <- nc; out$n <- sum(nc); out$cv <- sqrt(Vfun(nc)) / abs(Ybar)
      out$method <- "ceiling"
    }
  }
  out
}

#' Minimal integer allocation for a given stratification
#'
#' Finds the smallest integer sample size \eqn{n = \sum_h n_h}, with
#' \eqn{\min(m, N_h) \le n_h \le N_h} (\eqn{m} = `minnumstrat`), such that the
#' expected CV of the estimated mean (total) of every target variable does not
#' exceed its target. Single- and multi-domain stratifications are supported.
#'
#' @details
#' For each target variable \eqn{j} the constraint is
#' \deqn{V_j(n) = \sum_h W_h^2 \left(\frac{1}{n_h} - \frac{1}{N_h}\right) S_{jh}^2
#'   \le (CV_j \bar{Y}_j)^2,}
#' with \eqn{W_h = N_h / N} and \eqn{\bar{Y}_j = \sum_h W_h M_{jh}}
#' (Equation 1 of the paper).
#'
#' The algorithm is a marginal allocation:
#' \enumerate{
#'   \item start from the lower bounds \eqn{n_h = \min(m, N_h)};
#'   \item while some constraint is violated, add one unit to the stratum with
#'     the largest reduction of the violated variances, each normalised by its
#'     target: \eqn{\sum_{j\ violated} W_h^2 S_{jh}^2 (1/n_h - 1/(n_h+1)) / V_j^*};
#'   \item (several targets) remove, one at a time, any unit whose removal keeps
#'     all the constraints satisfied.
#' }
#' With a single target the result is the optimal integer allocation (marginal
#' gains decrease with \eqn{n_h}, so the greedy allocation is optimal for every
#' total). With several targets the result is a heuristic: it always satisfies
#' the constraints and, if `compare_ceiling = TRUE`, is never larger than the
#' Bethel allocation rounded up stratum by stratum.
#'
#' **Domains.** If `strata` has a column `DOM1` (case insensitive) and `errors`
#' is a data frame with a column `domainvalue`, the allocation is computed
#' separately in every domain with its own CV targets (the optimisation is
#' independent across domains, as in Bethel's algorithm). Otherwise a single
#' domain is assumed and only the first row of `errors` is used.
#'
#' @param strata Data frame with one row per stratum and columns `N`,
#'   `M1`..`Mk`, `S1`..`Sk` (and `DOM1` for several domains), as returned by
#'   [SamplingStrata::buildStrataDF()] or by the `aggr_strata` element of
#'   [SamplingStrata::optimStrata()]. Column names are not case sensitive.
#' @param errors Data frame with columns `CV1`..`CVk` (and `domainvalue` for
#'   several domains), or a numeric vector of `k` target CVs (single domain).
#' @param minnumstrat Minimum number of units per stratum (default 2).
#' @param divisor `"N_h-1"`: the `S_jh` are rescaled by
#'   \eqn{\sqrt{N_h/(N_h-1)}}, i.e. variances have divisor \eqn{N_h - 1} as in
#'   the common evaluation of the paper; `"N_h"`: the `S_jh` are used as given
#'   (SamplingStrata convention).
#' @param compare_ceiling Logical. If `TRUE` and SamplingStrata is available,
#'   the Bethel allocation rounded up is computed too and returned if smaller.
#'
#' @return A list. For a single domain:
#' \describe{
#'   \item{`n_h`}{integer allocation by stratum}
#'   \item{`n`}{total sample size}
#'   \item{`cv`}{achieved CVs of the target variables}
#'   \item{`n_ceiling`}{Bethel allocation rounded up (`NA` if not computed)}
#'   \item{`n_real`}{continuous Bethel sample size (`NA` if not computed)}
#'   \item{`method`}{`"greedy"` or `"ceiling"`}
#' }
#' For several domains:
#' \describe{
#'   \item{`n_h`}{allocation, aligned with the rows of `strata`}
#'   \item{`n`}{total sample size}
#'   \item{`by_domain`}{named list of single-domain results}
#'   \item{`n_dom`}{sample size by domain}
#'   \item{`cv`}{matrix (domains x targets) of achieved CVs}
#'   \item{`n_ceiling`, `n_real`}{totals over domains (`NA` if not computed)}
#'   \item{`method`}{method chosen in each domain}
#' }
#'
#' @seealso [SamplingStrata::optimStrata()], [SamplingStrata::bethel()]
#'
#' @examples
#' \dontrun{
#' # Evaluate solution
#' outstrata <- solution$aggr_strata
#' expected_CV(outstrata)
#' # Minimal integer solution
#' source("min_int_allocation.R")
#' # Domains are handled internally (DOM1 in outstrata, domainvalue in errors)
#' newalloc <- min_int_allocation(strata = outstrata, errors = errors, minnumstrat = 2,
#'                                divisor = "N_h", compare_ceiling = TRUE)
#' # Verify again
#' outstrata$SOLUZ <- newalloc$n_h
#' expected_CV(outstrata)
#' }
#' @export
min_int_allocation <- function(strata, errors, minnumstrat = 2,
                               divisor = c("N_h-1", "N_h"),
                               compare_ceiling = TRUE) {
  divisor <- match.arg(divisor)
  dom_s <- strata[[which(toupper(names(strata)) == "DOM1")[1]]]
  multi <- length(dom_s) > 0 && is.data.frame(errors) &&
    "domainvalue" %in% names(errors)
  if (!multi)
    return(.min_int_allocation_1dom(strata, errors, minnumstrat, divisor,
                                    compare_ceiling))

  doms <- sort(unique(dom_s))
  n_h <- rep(NA_real_, nrow(strata))
  res <- setNames(vector("list", length(doms)), as.character(doms))
  for (d in doms) {
    idx <- which(dom_s == d)
    e <- errors[as.character(errors$domainvalue) == as.character(d), ,
                drop = FALSE]
    if (nrow(e) == 0) stop("no precision constraints for domain ", d)
    res[[as.character(d)]] <- .min_int_allocation_1dom(
      strata[idx, , drop = FALSE], e, minnumstrat, divisor, compare_ceiling)
    n_h[idx] <- res[[as.character(d)]]$n_h
  }
  tot <- function(f) { v <- sapply(res, `[[`, f); if (anyNA(v)) NA else sum(v) }
  list(n_h = n_h, n = sum(n_h),
       by_domain = res,
       n_dom = sapply(res, `[[`, "n"),
       cv = do.call(rbind, lapply(res, `[[`, "cv")),
       n_ceiling = tot("n_ceiling"), n_real = tot("n_real"),
       method = sapply(res, `[[`, "method"))
}

