#' @title \code{ExtractXXXfromWORKING} a function to infer the country abbreviation by parsing the working directory file path.
#' 
#' @description this function is written in the name of modularity, since this otherwise takes place at the beginning of every function. When \code{XXX = NULL} (default). Clearly, if the last path element is for some reason NOT the HMD abbreviation, then \code{XXX} should be specified in the function arguments.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. For instance: \code{"/data/commons/hmd/HMDWORK/ISL"}, though this function has no default.
#' 
#' @return character vector of length 1" HMD letter code abbrevation for the country (inferred)
#'
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
# Author: triffe
###############################################################################
ExtractXXXfromWORKING <- function(WORKING){
  parts         <- rev(unlist(strsplit(WORKING, split = .Platform$file.sep)))
  parts[!parts == ""][1] 
}
