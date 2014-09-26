#'
#' @title AKm02q0 derive m0 from q0
#' 
#' @description Derive m0 from q0 according to the relevant segment of the Andreeve-Kingkade formula. This is elegant because it's an analytic solution, but ugly because, man, look at it. Carl Boe got this from Maple I think. This formula is only necessary because AK start with q0 whereas the HMD starts with m0, so we needed to adapt.
#' 
#' @param m0 a value of m0, per the period lifetable derivation
#' @param constant the intercept of the relevant Andreev-Kingkade segment
#' @param slope the slope of the relevant Andreev-Kingkade segment
#' 
#' @return q0 the estimate of q0 according to the identity between a0, m0, q0
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 

AKm02q0 <- function(m0,constant,slope){
  -1 / slope / m0 * (-m0 +  (m0 * constant) - 0.1e1 + sqrt(((m0 ^ 2) - 2 * constant * (m0 ^ 2) + 2 * m0 + (constant ^ 2) * (m0 ^ 2) - 2 *  (m0 * constant) + 1 - 4 * slope * (m0 ^ 2)))) / 2
}

#m0 <- .02
#constant <- 0.1493
#slope <- -2.0367
#t1  <- (m0 * constant)
#t2  <- (m0 ^ 2)
#t6  <- (constant ^ 2)
#t12 <- sqrt((t2 - 2 * constant * t2 + 2 * m0 + t6 * t2 - 2 * t1 + 1 - 4 * slope * t2))
#t19 <- -0.1e1 / slope / m0 * (-m0 + t1 - 0.1e1 + t12) / 0.2e1
#t19
# AKm02q0(m0,constant,slope)
#-1 / slope / m0 * (-m0 + m0 * constant - 1 + 
#    sqrt((m0 ^ 2 - 2 * constant * m0 ^ 2 + 2 * m0 + constant ^ 2 * m0 ^ 2 - 2 * m0 * constant + 1 - 4 * slope * m0 ^ 2))) / 2
#
#-1 / slope / m0 * (-m0 + m0 * constant - 1 + 
#    sqrt((m0 ^ 2 - 2 * constant * m0 ^ 2 + 2 * m0 + constant ^ 2 * m0 ^ 2 - 2 * m0 * constant + 1 - 4 * slope * m0 ^ 2))) / 2
# ----------------------------------------------------------
# this is the version in the MP:
#(1 / (2*m0*slope)) *
#   (-m0 + m0 * constant - 1 + sqrt(m0 ^ 2 - 2 * constant * m0 ^ 2 + 2 * m0 + constant ^ 2 * m0 ^ 2 - 2 * m0 * constant + 1 - 4 * slope * m0 ^ 2))
#




