#' Compute Normalized Moran's \eqn{I}
#'
#'A simple function to compute the normalized "\eqn{I_6}" Moran's \eqn{I} statistic proposed by Tillé et al. (2026), called by `normalized.moran.test` and `normalized.moran.mc`.
#'
#' @param x a numeric vector the same length as the neighbours list in listw
#' @param listw a `listw` object of style 'W' created for example by `nb2listw`
#' @param zero.policy default `attr(listw, "zero.policy")` as set when `listw` was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param NAOK if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.
#'
#' @returns The normalized "\eqn{I_6}" Moran's \eqn{I} of the values in `x` is computed, as a numeric vector of length one.
#' @seealso [normalized.moran.test], [normalized.moran.mc]
#' @section Author: Diego Giuliani \email{diego.giuliani@@unitn.it}
#' @section References: Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the Normalization of the Moran's \eqn{I} index: A Correlation-based Approach with Inference, \emph{International Statistical Review}
#' @export
#'
#' @examples
#' library(spdep)
#' data(oldcol)
#' col.W <- nb2listw(COL.nb, style="W")
#' crime <- COL.OLD$CRIME
#' normalized.moran(crime, col.W)
#' is.na(crime) <- sample(1:length(crime), 10)
#' normalized.moran(crime, col.W, NAOK=TRUE)
normalized.moran <- function (x, listw, zero.policy = attr(listw, "zero.policy"),
                              NAOK = FALSE)
{
  if (is.null(zero.policy))
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  x <- c(x)
  if (n1 != length(x))
    stop("objects of different length")
  xx <- mean(x, na.rm = NAOK)
  z <- x - xx
  lz <- spdep::lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  if(NAOK) {
    z_lz <- cbind(z,lz)
    z_lz <- na.omit(z_lz)
    z <- z_lz[,1]
    lz <- z_lz[,2]
  }
  I <- cor(z, lz)
  I
}
