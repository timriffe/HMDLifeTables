#' @title \code{CDa0} estimate a0 using the Coale-Demeny rule of thumb.
#'
#' @description
#' \code{CDa0} is an auxiliary function used by the four HMD lifetable functions, \code{ltper_AxN()}, \code{ltcoh_AxN()}, \code{ltperBoth_AxN()}, \code{ltcohBoth_AxN()}, described in formulas 61 and 62 of MP version 5.  
#'
#' @param m0 a value or vector of values of m0, the death rate for age 0 infants.
#' @param sex either "m" or "f"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
# Author: triffe
###############################################################################
CDa0 <- function(m0, sex ){
  sex <- rep(sex, length(m0))
  ifelse(sex == "m", 
    ifelse(m0 >= 0.107, 0.330, {0.045 + 2.684 * m0}), # males
    ifelse(m0 >= 0.107, 0.350, {0.053 + 2.800 * m0})  # females
  ) 
}

