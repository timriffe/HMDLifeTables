#' @title \code{checkRuncoh} a function to do an automatic check for whether cohort lifetables are run
#'
#' @description this function checks to see whether any cohorts have been observed from birth until \code{OPENAGE} and returns a list with logical elements for males, females, etc. The object created by this function is used in the cohort lifetable functions, \code{ltcoh_AxN()} and \code{ltcohBoth_AxN()}, which use the argument \code{run.condition} to decide which list element (\code{"m"}, code{"f"}, \code{"both"}, \code{"either"}) to extract.
#'
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param males the cohort components \code{data.frame} object, as created by \code{getCohortComponents()} for males. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param females the cohort components \code{data.frame} object, as created by getCohortComponents() for females. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param OPENAGE the desired open age. Default value is 110
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as \code{Rbin/checkRuncoh.Rdata} as well?
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the past path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a list with
#' \itemize{
#'  \item{"timestamp"}{e.g. "18 Dec 2012" *character class), just so you know whether the check was recent or not.}
#'  \item{"m"}{logical, males appear to have at least 1 extinct cohort}
#'  \item{"f"}{logical, females appear to have at least 1 extinct cohort}
#'  \item{"both"}{logical, both males and females appear to have at least 1 extinct cohort}
#'  \item{"either"}{logical, either males or females appear to have at least 1 extinct cohort}
#'  \item{"note"}{"TRUE indicates that a cohort appears extinct"}
#' }
#'
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
checkRuncoh <- function(WORKING = getwd(),  
  males = NULL,
  females = NULL,
  OPENAGE = 110,
  save.bin = TRUE, 
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL){
  
  # 'XXX' (country abbrev) is used here and there
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH <- file.path(WORKING, "InputDB")
  }
  extinct <- list(
    timestamp = format(Sys.time(), "%d %b %Y"),
    m = NULL,
    f = NULL,
    both = NULL,
    either = NULL,
    note = "TRUE indicates that a cohort appears extinct"
  )
  # run each for males, females sex = "m"
  if (is.null(males) | is.null(females)){
    cohcomp.path.m  <- file.path(WORKING, "Rbin", "m_cohortComponents.Rdata")
    cohcomp.path.f  <- file.path(WORKING, "Rbin", "f_cohortComponents.Rdata")
    if (file.exists(cohcomp.path.m) & file.exists(cohcomp.path.f)){ 
      males         <- local(get(load(cohcomp.path.m)))
      females       <- local(get(load(cohcomp.path.f)))
    } else {
      # load cohort components function
      # generate data
      males <- getCohortComponents(
                  WORKING = WORKING, 
                  sex = "m",
                  save.bin = TRUE,
                  XXX = XXX,
                  LDBPATH = LDBPATH,
                  IDBPATH = IDBPATH)
      females <- getCohortComponents(
                  WORKING = WORKING, 
                  sex = "f",
                  save.bin = TRUE,
                  XXX = XXX,
                  LDBPATH = LDBPATH,
                  IDBPATH = IDBPATH)
    }
  }
    
# ---------------------------------------------------------------------------------
# acast() transforms the data into a matrix with cohorts in the columns and ages in the rows.
# NAs are imputed for incomplete cohorts. require(reshape2)
    pop.m       <- acast(males, Age ~ Year, value.var = "pop")
    pop.f       <- acast(females, Age ~ Year, value.var = "pop")
# ---------------------------------------------------------------------------------
# cohort range not the same: dimension to be equal
# ---------------------------------------------------------------------------------
    # 1) limit to cohorts observed from birth (leaving intermediate NAs, e.g. BEL)
    cohs.all    <- as.integer(colnames(pop.m))
    FYNZ        <- cohs.all[!is.na(pop.m[1, ])][1]
    coh.keep1   <- as.character(FYNZ:max(cohs.all))
    pop.m       <- pop.m[, coh.keep1]
    pop.f       <- pop.f[, coh.keep1]
       
# ---------------------------------------------------------------------------------
    extinct[["m"]]      <- any(apply(pop.m, 2, function(.x){
                                if (any(.x == 0, na.rm = TRUE)){
                                  .x <- .x[1:(which(.x == 0)[1])]
                                }
                                # now check to make sure COMPLETELY observed (accounts for BEL)
                                any(.x == 0) & !any(is.na(.x))
                              }
                            )
                          )
    extinct[["f"]]      <- any(apply(pop.f, 2, function(.x){
                                if (any(.x == 0, na.rm = TRUE)){
                                  .x <- .x[1:(which(.x == 0)[1])]
                                }
                                # now check to make sure COMPLETELY observed (accounts for BEL)
                                any(.x == 0) & !any(is.na(.x))
                              }
                            )
                          )
  # determine if both or either is extinct
    extinct[["both"]]   <- extinct[["m"]] & extinct[["f"]]
    extinct[["either"]] <- extinct[["m"]] | extinct[["f"]]
  
  # saving routine
  if (save.bin){
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
    }
    out.name0 <- "checkRuncoh.Rdata"
    out.path0 <- file.path(dir.path, out.name0)
    
    # saving happens here:
    save(extinct, file = out.path0)
    # now save mx
    #Sys.chmod( out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }
  invisible(extinct)
}
