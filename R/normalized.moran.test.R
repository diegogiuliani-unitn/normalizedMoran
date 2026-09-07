#' Normalized Moran's \eqn{I} test for spatial autocorrelation
#'
#' Normalized "\eqn{I_6}" Moran's test for spatial autocorrelation proposed by Tillé et al. (2026) using a spatial weights matrix in weights list form. The test assumes asymptotic normality, and results may be checked against those of `normalized.moran.mc` permutations.
#'
#' @param x a numeric vector the same length as the neighbours list in listw
#' @param listw a `listw` object created for example by `nb2listw`
#' @param zero.policy default `attr(listw, "zero.policy")` as set when `listw` was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param alternative a character string specifying the alternative hypothesis, must be one of `greater` (default), `less` or `two.sided`
#' @param na.action a function (default `na.fail`), can also be `na.omit` or `na.exclude` - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set `zero.policy` to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the `glist` argument to `nb2listw` may be subsetted. If `na.pass` is used, zero is substituted for NA values in calculating the spatial lag
#' @param spChk should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use `get.spChkOption()`
#' @param adjust.n default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted
#'
#' @returns A list with class `htest` containing the following components: \describe{
#' \item{`statistic`}{the value of the standard deviate of normalized Moran's I.}
#' \item{`p.value`}{the p-value of the test.}
#' \item{`estimate`}{the value of the observed normalized Moran's I, its expectation and variance under normality.}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string giving the assumption used for calculating the standard deviate.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' }
#' @seealso [normalized.moran], [normalized.moran.mc]
#' @section Author: Diego Giuliani \email{diego.giuliani@@unitn.it}
#' @section References: Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the Normalization of the Moran's \eqn{I} index: A Correlation-based Approach with Inference, \emph{International Statistical Review}
#' @export
#'
#' @examples
#' library(spdep)
#' data(oldcol)
#' coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
#' normalized.moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
#' crime <- COL.OLD$CRIME
#' is.na(crime) <- sample(1:length(crime), 10)
#' res <- try(normalized.moran.test(crime, nb2listw(COL.nb, style="W"),
#'                                  na.action=na.fail))
#' normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
#'                       na.action=na.omit)
#' normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
#'                       na.action=na.exclude)
#' normalized.moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)
normalized.moran.test <- function (x, listw, zero.policy = attr(listw, "zero.policy"), alternative = "greater",
                                   na.action = na.fail,
                                   spChk = NULL, adjust.n = TRUE)
{
  alternative <- match.arg(alternative, c("greater", "less",
                                          "two.sided"))
  wname <- deparse(substitute(listw))
  if (!inherits(listw, "listw"))
    stop(wname, "is not a listw object")
  xname <- deparse(substitute(x))
  if (!is.numeric(x))
    stop(xname, " is not a numeric vector")
  if (is.null(attr(listw$weights, "W")))
    message(wname, " is not a listw object of style W, then it has been row-standardized")
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  stopifnot(length(zero.policy) == 1L)
  if (is.null(spChk))
    spChk <- get.spChkOption()
  stopifnot(is.logical(spChk))
  stopifnot(length(spChk) == 1L)
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  stopifnot(length(na.action) == 1L)
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  W <- listw2mat(listw)
  if (is.null(attr(listw$weights, "W")))
    listw <- spdep::mat2listw(W, style="W", zero.policy = zero.policy)
  I <- normalized.moran(x, listw, zero.policy = zero.policy,
                        NAOK = NAOK)
  one <- rep(1,n)
  In <- diag(rep(1,n))
  P <- In - one %*% t(one)/n
  D <- diag(rowSums(W))
  invD <- solve(D)
  B <- P %*% invD %*% W
  C <- t(W) %*% invD %*% P %*% invD %*% W
  B_tilde <- (B + t(B)) / 2

  tr <- function(M) sum(diag(M))
  trP <- tr(P)
  trC <- tr(C)
  trB <- tr(B)
  trB2_tilde <- tr(B_tilde %*% B_tilde)
  trP2 <- tr(P %*% P)
  trC2 <- tr(C %*% C)
  trBP <- tr(B_tilde %*% P)
  trBC <- tr(B_tilde %*% C)
  trPC <- tr(P %*% C)

  EI <- -1/sqrt((n-1)*trC)
  VI <- (2/(trP*trC)) * (trB2_tilde + trB^2/(4*trP^2)*trP2 + trB^2/(4*trC^2)*trC2 -
                           trB/trP*trBP - trB/trC*trBC +
                           trB^2/(2*trP*trC)*trPC)
  if (VI < 0)
    warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Normalized Moran I statistic standard deviate"
  if (alternative == "two.sided")
    PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater")
    PrI <- pnorm(ZI, lower.tail = FALSE)
  else PrI <- pnorm(ZI)
  if (!is.finite(PrI) || PrI < 0 || PrI > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  vec <- c(I, EI, VI)
  names(vec) <- c("Normalized Moran I statistic", "Expectation", "Variance")
  method <- "Normalized Moran I test under normality"
  data.name <- paste(xname, "\nweights:", wname, ifelse(is.null(na.act), "",
                                                        paste("\nomitted:", paste(na.act, collapse = ", "))),
                     ifelse(adjust.n && isTRUE(any(sum(card(listw$neighbours) ==
                                                         0L))), "\nn reduced by no-neighbour observations",
                            ""), "", "\n")
  res <- list(statistic = statistic, p.value = PrI, estimate = vec,
              alternative = alternative, method = method, data.name = data.name)
  if (!is.null(na.act))
    attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}
