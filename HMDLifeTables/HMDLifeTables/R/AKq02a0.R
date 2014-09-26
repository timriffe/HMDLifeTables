#' @title \code{AKq02a0} estimates a0 using the Andreev-Kingkade rule of thumb.
#'
#' @description \code{AKq02a0} is an auxiliary function used by version 6 of the four HMD lifetable functions, \code{ltper_AxN()}, \code{ltcoh_AxN()}, \code{ltperBoth_AxN()}, \code{ltcohBoth_AxN()}.
#'
#' @param q0 a value or vector of values of q0, the death probability for age 0 infants.
#' @param sex either "m" or "f"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export

AKq02a0 <- function(q0, sex = "m"){
  sex <- rep(sex, length(q0))
  ifelse(sex == "m", 
    ifelse(q0 < .0226, {0.1493 - 2.0367 * q0},
      ifelse(q0 < 0.0785, {0.0244 + 3.4994 * q0},.2991)),
    ifelse(q0 < 0.0170, {0.1490 - 2.0867 * q0},
      ifelse(q0 < 0.0658, {0.0438 + 4.1075 * q0}, 0.3141))
  )
}
