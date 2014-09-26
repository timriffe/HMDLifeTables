#' @title \code{ltper_GetDx100} determines the age at which Kannisto smoothed mx values should be imputed over the raw Mx for each year of data.
#' 
#' @description Smoothed mx values are imputed starting at the lowest age where male or female deaths drop to at most 100 with a minimum age of 80 and a maximum age of 95. The same age is used for both males and females.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param OPENAGE the desired open age. Default value is 110
#' @param N year interval: 1, 5 or 10. Other intervals would also work in theory.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}.
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}.
#' 
#' @return a vector with the age index (age + 1) for each year where smoothed mx values should be imputed.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
ltper_GetDx100 <- function(WORKING = getwd(), 
  perComp,
  OPENAGE = 110, 
  N = 1, 
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL){
  
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH       <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH       <- file.path(WORKING, "InputDB")
  }
  GetDx100sex <- function(WORKING, sex, OPENAGE, N, XXX, LDBPATH, IDBPATH){
   
      percomp.path <- file.path(WORKING, "Rbin", paste0(sex, "_periodComponents.Rdata"))
      if (file.exists(percomp.path)){
        perComp    <- local(get(load(percomp.path)))
      } else {
        perComp <- getPeriodComponents(
                      WORKING = WORKING, 
                      sex = sex,
                      OPENAGE = OPENAGE,
                      save.bin = TRUE,
                      XXX = XXX, 
                      LDBPATH = LDBPATH,
                      IDBPATH = IDBPATH)
      }
    
    dl        <- acast(perComp, Age ~ Year, value.var = "dl")
    du        <- acast(perComp, Age ~ Year, value.var = "du")
    Dx        <- dl + du
# ---------------------------------------------------------------------------------
# matlab does x <= 100 head(Dx)
    extrap.ages.i <- apply(Dx, 2, function(x, ages = 0:OPENAGE){                      
                           if (all(is.na(x))){
                             return(NA)
                           } else {
                             age.extrap <- ages[x <= 100 & ages >= 80 & ages <= 95][1]
                             ifelse(is.na(age.extrap), 95, age.extrap) + 1
                           }
                         }
                       )
    extrap.ages.i
  }
  
  indm <- GetDx100sex(WORKING = WORKING, 
                      sex = "m", 
                      OPENAGE = 110, 
                      N = N, 
                      XXX = XXX, 
                      LDBPATH = LDBPATH, 
                      IDBPATH = IDBPATH)
  indf <- GetDx100sex(WORKING = WORKING, 
                      sex = "f", 
                      OPENAGE = 110, 
                      N = N, 
                      XXX = XXX, 
                      LDBPATH = LDBPATH, 
                      IDBPATH = IDBPATH)
  pmin(indm, indf)
}
