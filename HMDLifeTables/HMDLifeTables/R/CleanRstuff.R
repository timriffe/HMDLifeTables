#' @title \code{CleanRstuff} DELETE the folders Rbin, RSTATS and RCHECKS permanently. Only used when doing testing. Function not to be included in shared package.
#'
#' @description Only used when doing testing. Function not to be included in shared package. Before running the full gamut of lifetable functions it is often desirable to have totally clean folders, since some functions use existing objects where possible. Likewise, the RCHECKS folder ought to be up to date, so this removes detritus. The lifetable functions overwrite output, so this function is not strictly necessary. The most important folder to have clear before a clean run (sometimes) would be Rbin. These folder names are fixed (Rbin, RSTATS, RCHECKS). Do this type of thing manually if you want anything else.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{NULL}, since this is a consequential function.
#' 
#' @return function used for its side effects. USE EXTREME CAUTION.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export
# Author: triffe
###############################################################################
CleanRstuff <- function(WORKING = NULL){
  # remove Rbin, RSTATS, RCHECKS folders 
  path1 <- file.path(WORKING, "Rbin")
  path2 <- file.path(WORKING, "RSTATS")
  path3 <- file.path(WORKING, "RCHECKS")
  if (file.exists(path1)){
    unlink(path1, recursive = TRUE)
  }
  if (file.exists(path2)){
    unlink(path2, recursive = TRUE)
  }
  if (file.exists(path3)){
    unlink(path3, recursive = TRUE)
  } 
}
