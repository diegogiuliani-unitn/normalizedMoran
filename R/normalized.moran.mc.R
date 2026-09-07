#' Permutation test for normalized Moran's \eqn{I} statistic
#'
#' A permutation test for the "\eqn{I_6}" normalized Moran's \eqn{I} statistic by Tillé et al. (2026) calculated by using `nsim` random permutations of `x` for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the `nsim` simulated values.
#'
#' @param x a numeric vector the same length as the neighbours list in `listw`
#' @param listw a `listw` object created for example by nb2listw
#' @param nsim number of permutations
#' @param zero.policy default `attr(listw, "zero.policy")` as set when `listw` was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param alternative a character string specifying the alternative hypothesis, must be one of `greater` (default), `less` or `two.sided`
#' @param na.action a function (default `na.fail`), can also be `na.omit` or `na.exclude` - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set `zero.policy` to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the `glist` argument to `nb2listw` may be subsetted. If `na.pass` is used, zero is substituted for NA values in calculating the spatial lag
#' @param spChk should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use `get.spChkOption()`
#' @param return_boot return an object of class `boot` from the equivalent permutation bootstrap rather than an object of class `htest`
#' @param adjust.n default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted
#'
#' @returns A list with class `htest` and `mc.sim` containing the following components: \describe{
#' \item{`statistic`}{the value of the observed normalized Moran's I.}
#' \item{`parameter`}{the rank of the observed normalized Moran's I.}
#' \item{`p.value`}{the pseudo p-value of the test.}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string giving the method used.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' \item{`res`}{`nsim` simulated values of statistic, final value is observed statistic.}
#' }
#' @seealso [normalized.moran], [normalized.moran.test]
#' @section Author: Diego Giuliani \email{diego.giuliani@@unitn.it}
#' @section References: Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the Normalization of the Moran's \eqn{I} index: A Correlation-based Approach with Inference, \emph{International Statistical Review}
#' @export
#'
#' @examples
#' library(spdep)
#' data(oldcol)
#' colw <- nb2listw(COL.nb, style="W")
#' nsim <- 99
#' set.seed(1234)
#' sim1 <- normalized.moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
#' sim1
#' mean(sim1$res[1:nsim])
#' var(sim1$res[1:nsim])
#' summary(sim1$res[1:nsim])
#' colold.lags <- nblag(COL.nb, 3)
#' set.seed(1234)
#' sim2 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
#'                                                     style="W"), nsim=nsim)
#' summary(sim2$res[1:nsim])
#' sim3 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
#'                                                     style="W"), nsim=nsim)
#' summary(sim3$res[1:nsim])
#' crime <- COL.OLD$CRIME
#' is.na(crime) <- sample(1:length(crime), 10)
#' try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
#'                         na.action=na.fail))
#' normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
#'                     na.action=na.omit)
#' normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
#'                     return_boot=TRUE, na.action=na.omit)
#' normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
#'                     na.action=na.exclude)
#' normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
#'                     return_boot=TRUE, na.action=na.exclude)
#' try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
normalized.moran.mc <- function(x, listw, nsim, zero.policy = attr(listw, "zero.policy"),
                                alternative = "greater", na.action = na.fail, spChk = NULL,
                                return_boot = FALSE, adjust.n = TRUE)
{
  alternative <- match.arg(alternative, c("greater", "less",
                                          "two.sided"))
  wname <- deparse(substitute(listw))
  if (!inherits(listw, "listw"))
    stop(wname, "is not a listw object")
  xname <- deparse(substitute(x))
  if (!is.numeric(x))
    stop(xname, "is not a numeric vector")
  if (is.null(attr(listw$weights, "W")))
    message(wname, " is not a listw object of style W, then it has been row-standardized")
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  stopifnot(length(zero.policy) == 1L)
  if (missing(nsim))
    stop("nsim must be given")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  stopifnot(is.logical(spChk))
  stopifnot(length(spChk) == 1L)
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  cards <- card(listw$neighbours)
  if (!zero.policy && any(cards == 0))
    stop("regions with no neighbours found")
  stopifnot(length(na.action) == 1L)
  if (deparse(substitute(na.action)) == "na.pass")
    stop("na.pass not permitted")
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    if (return_boot)
      message("NA observations omitted: ", paste(na.act,
                                                 collapse = ", "))
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  gamres <- suppressWarnings(nsim > gamma(n + 1))
  if (gamres)
    stop("nsim too large for this number of observations")
  if (nsim < 1)
    stop("nsim too small")
  if (adjust.n)
    n <- n - sum(cards == 0L)
  if (is.null(attr(listw$weights, "W"))) {
    W <- spdep::listw2mat(listw)
    listw <- spdep::mat2listw(W, style="W", zero.policy = zero.policy)
  }
  if (return_boot) {
    moran_boot <- function(var, i, ...) {
      var <- var[i]
      return(normalized.moran(x = var, ...))
    }
    p_setup <- spdep:::parallel_setup(NULL)
    parallel <- p_setup$parallel
    ncpus <- p_setup$ncpus
    cl <- p_setup$cl
    res <- boot::boot(x, statistic = moran_boot, R = nsim, sim = "permutation",
                      listw = listw, zero.policy = zero.policy,
                      parallel = parallel, ncpus = ncpus, cl = cl)
    return(res)
  }
  res <- numeric(length = nsim + 1)
  for (i in 1:nsim) res[i] <- normalized.moran(sample(x), listw,
                                               zero.policy)
  res[nsim + 1] <- normalized.moran(x, listw, zero.policy)
  rankres <- rank(res)
  xrank <- rankres[length(res)]
  diff <- nsim - xrank
  diff <- ifelse(diff > 0, diff, 0)
  if (alternative == "less")
    pval <- punif((diff + 1)/(nsim + 1), lower.tail = FALSE)
  else if (alternative == "greater")
    pval <- punif((diff + 1)/(nsim + 1))
  else pval <- punif(abs(xrank - (nsim + 1)/2)/(nsim + 1),
                     0, 0.5, lower.tail = FALSE)
  if (!is.finite(pval) || pval < 0 || pval > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  statistic <- res[nsim + 1]
  names(statistic) <- "statistic"
  parameter <- xrank
  names(parameter) <- "observed rank"
  method <- "Monte-Carlo simulation of Normalized Moran I"
  data.name <- paste(xname, "\nweights:", wname, ifelse(is.null(na.act),
                                                        "", paste("\nomitted:", paste(na.act, collapse = ", "))),
                     "\nnumber of simulations + 1:", nsim + 1, "\n")
  lres <- list(statistic = statistic, parameter = parameter,
               p.value = pval, alternative = alternative, method = method,
               data.name = data.name, res = res)
  if (!is.null(na.act))
    attr(lres, "na.action") <- na.act
  class(lres) <- c("htest", "mc.sim")
  lres
}
