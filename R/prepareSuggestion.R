# ==============================================================================
# File:      prepareSuggestion.R
# Functions: prepareSuggestion() [exported]; evalSuggestion() [exported];
#            .ps2_cuts() [internal]
#
# Purpose:
#   prepareSuggestion(): builds the two initial chromosomes (cut points of the X
#   variables, normalised by their maximum within each domain) to be passed as
#   "suggestions" to optimStrata(method = "continuous"). Seeds: 1-D k-means,
#   cumulative root frequency, quantiles, a single "univariate" variable, or a
#   translation of a KmeansSolution2() solution; X variables can be made
#   inactive.
#   evalSuggestion(): partition implied by one chromosome (same rule as
#   strataGenalg2()) and its Bethel sample size in each domain.
#   .ps2_cuts(): cut points of one variable according to a given rule.
# ==============================================================================

#' Initial solutions (suggestions) for the "continuous" method
#'
#' @description
#' Builds two initial solutions (chromosomes) to be passed, through the
#' argument \code{suggestions}, to \code{\link{optimStrata}} with
#' \code{method = "continuous"} (i.e. to \code{\link{optimizeStrata2}}).
#' Each chromosome contains, for every domain, the cut points of the
#' stratification variables X1, X2, ..., expressed as fractions of the maximum
#' of each variable within the domain. This function replaces the previous
#' version of \code{prepareSuggestion}.
#'
#' @details
#' \strong{How the "continuous" method forms the strata.}
#' In \code{strataGenalg2()} (function \code{evaluate()}) the cut points of each
#' variable X_j split its range into H closed intervals; a unit equal to a cut
#' point falls in the upper interval. The stratum of a unit is the
#' \emph{maximum}, over the variables, of the index of the interval in which it
#' falls. With several X variables the strata are therefore nested
#' "staircase" regions, not the cross-classification of the intervals.
#'
#' \strong{What this function does.}
#' \itemize{
#'   \item the seeds are built directly in the space of the X variables, so
#'     that each chromosome is an exact, valid partition under the rule above;
#'   \item cut points are normalised by the maximum of X_j \emph{within each
#'     domain}, as \code{evaluate()} does (the previous version used the
#'     overall maximum);
#'   \item cut points are placed half-way between the last value of a stratum
#'     and the next observed value, so that units on a boundary are not
#'     shifted to the upper stratum;
#'   \item the two chromosomes can come from different rules (default: 1-D
#'     k-means on each X, and the cumulative root frequency rule of Dalenius
#'     and Hodges on a single X).
#' }
#'
#' \strong{Active and inactive variables.}
#' Cutting a variable that is irrelevant for the targets is harmful under the
#' "maximum index" rule: its cuts push units into high strata almost at
#' random. Only the variables listed in \code{active} (default: all) receive
#' cut points; the cuts of the other variables are set to 1 (i.e. at their
#' maximum), which makes them practically inactive (only the unit(s) equal to
#' the maximum of an inactive variable are moved to the top stratum). The
#' genetic algorithm remains free to activate them during the search.
#'
#' \strong{Available seeds} (argument \code{seeds}):
#' \describe{
#'   \item{\code{"univariate"}}{cut points on \emph{one} active variable only
#'     (\code{univariate_var}, by default the active X with the highest mean
#'     absolute correlation with the target variables Y* in \code{frame}),
#'     computed with the rule \code{univariate_rule}; all the other variables
#'     are inactive. Since the genetic algorithm keeps its best solution
#'     (elitism), it cannot end worse than this univariate stratification.}
#'   \item{\code{"kmeans1d"}}{1-D k-means with H clusters on each X_j (the
#'     clusters of a single variable are intervals, so the translation into
#'     cut points is exact).}
#'   \item{\code{"cumrootf"}}{cumulative root frequency rule (Dalenius and
#'     Hodges, 1959) on each X_j.}
#'   \item{\code{"quantile"}}{equal-frequency intervals on each X_j.}
#'   \item{\code{"kmean"}}{approximate translation of a
#'     \code{\link{KmeansSolution2}} solution (argument \code{kmean}), close to
#'     the previous \code{prepareSuggestion}: for each X_j the clusters are
#'     ordered by the median of X_j and their maxima are used as upper
#'     bounds.}
#' }
#'
#' \strong{Limitations.}
#' \itemize{
#'   \item the X variables must be non-negative (the genetic algorithm encodes
#'     the cut points in [0, 1] times the maximum of X_j);
#'   \item \code{optimizeStrata2()} counts as X variables \emph{all} the
#'     columns whose (upper case) name contains the letter "X": the frame
#'     should not contain other such columns (a warning is issued).
#' }
#'
#' @param frame The sampling frame, as produced by \code{\link{buildFrameDF}}:
#'   a data frame with the stratification variables X1, X2, ... (numeric,
#'   non-negative, no missing values), the target variables Y1, Y2, ... and
#'   the column \code{domainvalue} (consecutive integers 1, ..., D). The
#'   column \code{id} is required only by the seed \code{"kmean"}. Column
#'   names are not case sensitive.
#' @param nStrata Number of strata H: either a single integer, used for all
#'   domains, or a vector with one value (at least 2) per domain.
#' @param seeds Rules used to build the two chromosomes: one or two of
#'   \code{"kmeans1d"}, \code{"cumrootf"}, \code{"quantile"}, \code{"kmean"},
#'   \code{"univariate"}. If only one rule is given, it is used for both
#'   chromosomes. Default \code{c("kmeans1d", "univariate")}.
#' @param active Character vector with the names of the X variables that
#'   receive cut points (e.g. \code{c("X1", "X4")}). Default \code{NULL}: all
#'   the X variables are active.
#' @param univariate_var Name of the X variable used by the seed
#'   \code{"univariate"}; it must be one of the active variables. Default
#'   \code{NULL}: the active X with the highest mean absolute correlation with
#'   the Y variables (the first active X if there are no Y variables).
#' @param univariate_rule Rule used by the seed \code{"univariate"}: one of
#'   \code{"cumrootf"} (default), \code{"kmeans1d"}, \code{"quantile"}.
#' @param kmean Data frame with columns \code{id} and \code{suggestions}, as
#'   returned by \code{\link{KmeansSolution2}}; required only by the seed
#'   \code{"kmean"}.
#' @param nclass Number of equal-width classes used to compute the cumulative
#'   root frequencies in the rule \code{"cumrootf"} (default 1000; reduced to
#'   the number of distinct values when this is smaller).
#'
#' @return A data frame with one row per cut point (H - 1 for each X
#'   variable, in the order X1, X2, ..., within each domain) and columns
#'   \code{domainvalue}, \code{suggestions1} (first chromosome) and
#'   \code{suggestions2} (second chromosome), with values in [0, 1]. The
#'   attribute \code{"seeds"} contains a short description of the rule that
#'   generated each chromosome.
#'
#' @references
#' Dalenius, T. and Hodges, J. L. (1959). Minimum variance stratification.
#' \emph{Journal of the American Statistical Association}, 54(285), 88-101.
#'
#' @seealso \code{\link{evalSuggestion}}, \code{\link{optimStrata}},
#'   \code{\link{KmeansSolution2}}, \code{\link{buildFrameDF}}
#'
#' @examples
#' # ---------------------------------------------------------------------
#' # Preparation of the frame from the "iris" dataset
#' # ---------------------------------------------------------------------
#' # Domains: the three species (50 units each)
#' # Stratification variables (X): Petal.Length, Petal.Width, Sepal.Width
#' # Target variables (Y):         Sepal.Length, Petal.Length
#' data(iris)
#' iris$id <- seq_len(nrow(iris))
#' frame <- buildFrameDF(df = iris,
#'                       id = "id",
#'                       X = c("Petal.Length", "Petal.Width", "Sepal.Width"),
#'                       Y = c("Sepal.Length", "Petal.Length"),
#'                       domainvalue = "Species")
#' head(frame)
#' # Precision constraints: CV of 2% on Y1 and 3% on Y2 in each species
#' errors <- data.frame(DOM = "DOM1",
#'                      CV1 = 0.02,
#'                      CV2 = 0.03,
#'                      domainvalue = 1:3)
#' errors
#' # ---------------------------------------------------------------------
#' # Suggestions: only X1 (Petal.Length) and X2 (Petal.Width) are active
#' # ---------------------------------------------------------------------
#' set.seed(1234)    # the seed "kmeans1d" uses stats::kmeans()
#' sugg <- prepareSuggestion(frame, nStrata = 3, active = c("X1", "X2"))
#' attr(sugg, "seeds")                          # description of the seeds
#' evalSuggestion(frame, errors, sugg, which = 1)  # 1-D k-means seed
#' evalSuggestion(frame, errors, sugg, which = 2)  # univariate seed
#' \donttest{
#' # ---------------------------------------------------------------------
#' # Optimisation with the "continuous" method starting from the seeds
#' # ---------------------------------------------------------------------
#' set.seed(1234)
#' solution <- optimStrata(method = "continuous",
#'                         framesamp = frame,
#'                         errors = errors,
#'                         nStrata = c(3, 3, 3),
#'                         suggestions = sugg,
#'                         iter = 50,
#'                         pops = 10,
#'                         writeFiles = FALSE,
#'                         showPlot = FALSE,
#'                         parallel = FALSE)
#' sum(solution$aggr_strata$SOLUZ)              # total sample size
#' }
#'
#' @export
prepareSuggestion <- function(frame,
                               nStrata,
                               seeds = c("kmeans1d", "univariate"),
                               active = NULL,
                               univariate_var = NULL,
                               univariate_rule = c("cumrootf", "kmeans1d", "quantile"),
                               kmean = NULL,
                               nclass = 1000) {
  frame <- as.data.frame(frame)
  colnames(frame) <- toupper(colnames(frame))
  if (!"DOMAINVALUE" %in% colnames(frame))
    stop("'frame' must contain the column 'domainvalue'")

  xvars <- grep("^X[0-9]+$", colnames(frame), value = TRUE)
  xvars <- xvars[order(as.integer(sub("^X", "", xvars)))]
  nX <- length(xvars)
  if (nX == 0) stop("No stratification variables X1, X2, ... in 'frame'")
  if (length(grep("X", colnames(frame))) != nX)
    warning("optimizeStrata2() counts every column whose name contains 'X' ",
            "as a stratification variable: remove or rename the columns ",
            paste(setdiff(grep("X", colnames(frame), value = TRUE), xvars),
                  collapse = ", "))
  if (anyNA(frame[xvars])) stop("Missing values in the X variables")
  if (any(frame[xvars] < 0)) stop("The X variables must be non-negative")

  seeds <- match.arg(seeds, c("kmeans1d", "cumrootf", "quantile", "kmean",
                              "univariate"), several.ok = TRUE)
  if (length(seeds) == 1) seeds <- rep(seeds, 2)
  if (length(seeds) != 2) stop("'seeds' must contain one or two rules")
  univariate_rule <- match.arg(univariate_rule)

  # active variables
  if (is.null(active)) active <- xvars
  active <- toupper(active)
  if (!all(active %in% xvars))
    stop("'active' must be a subset of ", paste(xvars, collapse = ", "))
  # variable of the univariate seed: the active X most correlated with the Y*
  if ("univariate" %in% seeds) {
    if (is.null(univariate_var)) {
      yvars <- grep("^Y[0-9]+$", colnames(frame), value = TRUE)
      if (length(yvars) == 0) {
        univariate_var <- active[1]
      } else {
        sc <- vapply(active, function(x)
          mean(abs(suppressWarnings(stats::cor(frame[[x]], frame[yvars])))),
          numeric(1))
        sc[is.na(sc)] <- 0
        univariate_var <- active[which.max(sc)]
      }
    }
    univariate_var <- toupper(univariate_var)
    if (!univariate_var %in% active)
      stop("'univariate_var' must be one of the active variables")
  }

  domains <- sort(unique(as.integer(as.character(frame$DOMAINVALUE))))
  if (!identical(domains, seq_along(domains)))
    stop("'domainvalue' must contain consecutive integers 1, ..., D")
  nStrata <- as.integer(nStrata)
  if (length(nStrata) == 1 && length(domains) > 1)
    nStrata <- rep(nStrata, length(domains))
  if (length(nStrata) != length(domains))
    stop("'nStrata' must have one element per domain")
  if (any(nStrata < 2)) stop("'nStrata' must be at least 2 in every domain")

  km_lab <- NULL
  if ("kmean" %in% seeds) {
    if (is.null(kmean)) stop("Seed 'kmean' requires the argument 'kmean'")
    if (!"ID" %in% colnames(frame)) stop("Seed 'kmean' requires the column 'id' in 'frame'")
    colnames(kmean) <- tolower(colnames(kmean))
    km_lab <- as.character(kmean$suggestions)[match(as.character(frame$ID),
                                                    as.character(kmean$id))]
    if (anyNA(km_lab)) stop("Some units of 'frame' are not in 'kmean'")
  }

  out <- vector("list", length(domains))
  for (d in domains) {
    sel <- as.integer(as.character(frame$DOMAINVALUE)) == d
    H <- nStrata[d]
    if (!is.null(km_lab)) {
      nk <- length(unique(km_lab[sel]))
      if (nk != H)
        warning("Domain ", d, ": 'kmean' has ", nk, " clusters but nStrata = ", H)
    }
    s1 <- s2 <- NULL
    for (x in xvars) {
      v <- frame[[x]][sel]
      mx <- max(v)
      lab <- if (is.null(km_lab)) NULL else km_lab[sel]
      cc <- lapply(seeds, function(sd) {
        if (!x %in% active) return(rep(mx, H - 1))           # inactive
        if (sd == "univariate") {
          if (x != univariate_var) return(rep(mx, H - 1))    # inactive
          return(.ps2_cuts(v, H, univariate_rule, lab, nclass))
        }
        .ps2_cuts(v, H, sd, lab, nclass)
      })
      c1 <- cc[[1]]
      c2 <- cc[[2]]
      if (mx > 0) {
        c1 <- c1 / mx
        c2 <- c2 / mx
      } else {                       # constant zero variable: no information
        c1 <- c2 <- rep(1, H - 1)
      }
      s1 <- c(s1, pmin(pmax(c1, 0), 1))
      s2 <- c(s2, pmin(pmax(c2, 0), 1))
    }
    out[[d]] <- data.frame(domainvalue = d, suggestions1 = s1, suggestions2 = s2)
  }
  res <- do.call(rbind, out)
  seed_label <- function(sd) if (sd == "univariate")
    paste0("univariate ", univariate_rule, " on ", univariate_var) else
    paste0(sd, " on ", paste(active, collapse = "+"))
  attr(res, "seeds") <- vapply(seeds, seed_label, character(1), USE.NAMES = FALSE)
  res
}


# Upper bounds of the H - 1 lower strata of variable v according to a rule,
# converted into cut points compatible with strataGenalg2()
#' @noRd
.ps2_cuts <- function(v, H, rule, lab = NULL, nclass = 1000) {
  u <- sort(unique(v))
  if (length(u) <= H) {
    upper <- u[seq_len(min(H - 1, length(u)))]
  } else if (rule == "quantile") {
    upper <- stats::quantile(v, probs = seq_len(H - 1) / H, type = 1, names = FALSE)
  } else if (rule == "cumrootf") {
    nc <- min(nclass, length(u))
    br <- seq(u[1], u[length(u)], length.out = nc + 1)
    f <- tabulate(findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE),
                  nbins = nc)
    csf <- cumsum(sqrt(f))
    k <- vapply(csf[nc] * seq_len(H - 1) / H,
                function(t) which(csf >= t)[1], integer(1))
    upper <- vapply(br[k + 1], function(e) max(v[v <= e]), numeric(1))
  } else if (rule == "kmeans1d") {
    km <- stats::kmeans(v, centers = H, nstart = 10, iter.max = 100)
    ord <- order(km$centers)
    upper <- vapply(ord[-H], function(k) max(v[km$cluster == k]), numeric(1))
  } else if (rule == "kmean") {
    med <- tapply(v, lab, stats::median)
    ord <- names(sort(med))
    mxs <- vapply(ord, function(k) max(v[lab == k]), numeric(1))
    upper <- cummax(mxs)[seq_len(min(H - 1, length(mxs)))]
  }
  upper <- sort(upper)
  if (length(upper) < H - 1)                       # pad: empty upper strata
    upper <- c(upper, rep(u[length(u)], H - 1 - length(upper)))
  # cut half-way between the upper bound and the next observed value, so that
  # the unit equal to the upper bound stays in the lower stratum
  vapply(upper, function(b) {
    nxt <- u[u > b]
    if (length(nxt) == 0) b else (b + nxt[1]) / 2
  }, numeric(1))
}


#' Evaluation of a suggestion produced by prepareSuggestion
#'
#' @description
#' Computes the partition of the frame implied by one of the two chromosomes
#' returned by \code{\link{prepareSuggestion}}, using the same rule as
#' \code{strataGenalg2()}, and the corresponding Bethel sample size in each
#' domain. Useful to check the seeds before running the genetic algorithm,
#' and to compare them with its final solution.
#'
#' @param frame The sampling frame used to build \code{sugg} (see
#'   \code{\link{prepareSuggestion}}).
#' @param errors Data frame of precision constraints (columns \code{DOM},
#'   \code{CV1}, \code{CV2}, ..., \code{domainvalue}), as required by
#'   \code{\link{optimStrata}}.
#' @param sugg Data frame returned by \code{\link{prepareSuggestion}}.
#' @param which Chromosome to evaluate: 1 (column \code{suggestions1}) or 2
#'   (column \code{suggestions2}).
#' @param minnumstr Minimum number of units allocated to each stratum
#'   (passed to \code{\link{bethel}} as \code{minnumstrat}). Default 2.
#'
#' @return A data frame with one row per domain and columns
#'   \code{domainvalue}, \code{strata_requested} (H),
#'   \code{strata_nonempty} (number of non-empty strata produced by the
#'   chromosome) and \code{n} (Bethel sample size, rounded to one decimal).
#'
#' @seealso \code{\link{prepareSuggestion}}, \code{\link{bethel}}
#'
#' @examples
#' data(iris)
#' iris$id <- seq_len(nrow(iris))
#' frame <- buildFrameDF(df = iris, id = "id",
#'                       X = c("Petal.Length", "Petal.Width", "Sepal.Width"),
#'                       Y = c("Sepal.Length", "Petal.Length"),
#'                       domainvalue = "Species")
#' errors <- data.frame(DOM = "DOM1", CV1 = 0.02, CV2 = 0.03, domainvalue = 1:3)
#' set.seed(1234)    # the seed "kmeans1d" uses stats::kmeans()
#' sugg <- prepareSuggestion(frame, nStrata = 3, active = c("X1", "X2"))
#' evalSuggestion(frame, errors, sugg, which = 1)
#' evalSuggestion(frame, errors, sugg, which = 2)
#'
#' @export
evalSuggestion <- function(frame, errors, sugg, which = 1, minnumstr = 2) {
  frame <- as.data.frame(frame)
  colnames(frame) <- toupper(colnames(frame))
  colnames(errors) <- toupper(colnames(errors))
  xvars <- grep("^X[0-9]+$", colnames(frame), value = TRUE)
  xvars <- xvars[order(as.integer(sub("^X", "", xvars)))]
  nX <- length(xvars)
  col <- paste0("suggestions", which)
  res <- NULL
  for (d in sort(unique(as.integer(as.character(frame$DOMAINVALUE))))) {
    fr <- frame[as.integer(as.character(frame$DOMAINVALUE)) == d, ]
    sg <- sugg[sugg$domainvalue == d, col]
    H <- length(sg) / nX + 1
    lab <- rep(1L, nrow(fr))
    for (j in seq_len(nX)) {
      x <- fr[[xvars[j]]]
      cuts <- sort(sg[((j - 1) * (H - 1) + 1):(j * (H - 1))] * max(x))
      lab <- pmax(lab, findInterval(x, cuts) + 1L)   # ties go to the upper stratum
    }
    fr2 <- fr[, setdiff(colnames(fr), xvars), drop = FALSE]
    fr2$X1 <- as.integer(factor(lab))
    st <- SamplingStrata::buildStrataDF(fr2, progress = FALSE, verbose = FALSE)
    err <- errors[as.integer(as.character(errors$DOMAINVALUE)) == d, , drop = FALSE]
    nh <- SamplingStrata::bethel(st, err, minnumstrat = minnumstr,
                                 realAllocation = TRUE)
    res <- rbind(res, data.frame(domainvalue = d, strata_requested = H,
                                 strata_nonempty = length(unique(lab)),
                                 n = round(sum(nh), 1)))
  }
  res
}
