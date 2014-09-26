#' @title \code{Write_Births} a function to prepare, format and write \code{Births.txt}.
#'
#' @description This function pulls birth count data from the Lexis Database (age = 0, lexis = 1) and sets up the standard output with columns \code{Year}, \code{Female}, \code{Male}, and \code{Total}. This function makes no consequential calls to other functions, and is not called by any other function, nor does it require another function to have been run previously.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param PVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' 
#' @return function called for its side effect of creating the file \code{Births.txt}. No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export
# Author: triffe
###############################################################################
Write_Births <- function(
  WORKING = getwd(),
  STATSFOLDER = "RSTATS",
  LDBPATH = NULL,
  MPVERSION = 5,
  XXX = NULL){
  
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH       <- file.path(WORKING, "LexisDB")
  }
  # Lexis DB object
  ldb.path.f      <- file.path(LDBPATH, paste0("f", XXX, ".txt"))
  LDBobj.f        <- read.table(ldb.path.f, header = FALSE, sep = ",", 
                        col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  ldb.path.m      <- file.path(LDBPATH, paste0("m", XXX, ".txt"))
  LDBobj.m        <- read.table(ldb.path.m, header = FALSE, sep = ",", 
                        col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  
  Bf              <- with(LDBobj.f, Population[Age == 0 & Lexis == 1 & Population > -1])
  Bm              <- with(LDBobj.m, Population[Age == 0 & Lexis == 1 & Population > -1])
  Bt              <- Bf + Bm
  Year            <- with(LDBobj.m, Year[Age == 0 & Lexis == 1 & Population > -1])
 
  # for the metadata header: country long name
  # country.lookup should load; if throws error, try data(country.lookup)
  CountryLong     <- country.lookup[country.lookup[,1] == XXX,2]
  DateMod         <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ",")
  # Methods Protocol version
  MPvers          <- ifelse(MPVERSION == 5, " MPv5 (May07)", "MPv6 (in development)\n")
  DataType        <- ",  Births (1-year)"
 
  # save formatted .txt out to this folder, make sure exists
  STATS.path      <- file.path(WORKING, STATSFOLDER)
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #system(paste0("chgrp hmdcalc ", STATS.path))
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
  }
  write.out.file  <- file.path(STATS.path, "Births.txt")
  
  cat(
    # metadata header
    paste0(CountryLong, DataType, DateMod, MPvers,"\n"),
    # column headers, copied and pasted
    "Year          Female      Male     Total",
    # the data, rounded and formatted in place- no tabbing
    paste(
      sprintf(paste0("%-", 10, "s"), Year),
      sprintf(paste0("%", 10, "s"), Bf),
      sprintf(paste0("%", 10, "s"), Bm),
      sprintf(paste0("%", 10, "s"), Bt),
      sep = ""), 
    file = write.out.file, sep = "\n")
  #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
  #system(paste0("chgrp hmdcalc ", write.out.file))
}

