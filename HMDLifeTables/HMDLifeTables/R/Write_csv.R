#' @title \code{Write_csv} a function to write formatted .csv files containing
#'   1x1 lifetable data, plus occurances and exposures.

#' 
#' @description This function requires that all pertinent \code{ltper_AxN.Rdata}
#'   and \code{ltcoh_AxN()} data objects be present in the folder 
#'   \code{WORKING/Rbin/} for males, females and both-sex tables. Objects are 
#'   selected by grepping, which works for the present HMD file naming scheme. 
#'   \code{.Rdata} files are read in, rounded, formatted and written back out 
#'   into \code{WORKING/RSTATS/} by default, although this can be changed with 
#'   the \code{STATSFOLDER} argument.
#'   
#' @param WORKING path to working directory, which typically ends with the HMD 
#'   country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a 
#'   full path). Default \code{"RSTATS"}.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects file headers.
#'   
#'   
#' @param OldStyle.  Logical.  When true, output is .csv files analogous to
#'   those produced under Matlab V5 code (no headers), that are part of the
#'   input for SAS diagnostics.  When FALSE (default, not yet implemented) the output will be ....

#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is
#'   extracted from \code{WORKING} as the last path part.
#'   
#'   
#' @return function called for its side effect of creating the lifetable txt
#'   output files, e.g. \code{mltper_1x1.txt} and other time/sex configurations.
#'   No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export 
#  Author: cboe
###############################################################################
Write_csv <- function(
  WORKING = getwd(), 
  STATSFOLDER = "RSTATS", 
  MPVERSION = 5,
  OldStyle=TRUE,  # compatible with SAS diagnostics
  XXX = NULL){
  
  # MatlabRound() is for rounding output, should give same result as matlab, assuming that's important
  # by CB, updated by TR to take digits as arg.
  MatlabRoundFW <- function(x, digits = 0, pad = TRUE, Age = FALSE, totalL = 8){ 
    
    # this 1) rounds, and
    # 2) makes sure that the zeros stay on the end of the number string to the specified number of digits
    if (is.numeric(x)){
      fac         <- rep(10 ^ digits, length(x))
      x           <- sprintf(paste0("%.", digits, "f"), floor(x * fac + sign(x) * 0.5) / fac)
    }
    # strings will potentially vary in length due to integers: adds white space out to the
    # longest string (this is relevant e.g. for dx). ensures digit alignment.
    if (pad){
      maxw        <- max(nchar(x))  
      x           <- sprintf(paste0("%", ifelse(Age, {maxw - 1}, maxw), "s"), x)
      x           <- sprintf(paste0("%-", maxw, "s"), x) # this only had affect if Age == TRUE (dealing with '+')
    }
    # add optional left padding to specify total character width
    x             <- sprintf(paste0("%", totalL, "s"), x)
    x
  }
  # get country abbrev
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  # for the metadata header: country long name
  # ltper_AxN() should have been run with save.bin = TRUE, 
  # such that there's an Rbin folder waiting there with stuff in it
  Rbin.path       <- file.path(WORKING, "Rbin")
  # save formatted .txt out to this folder, make sure exists
  STATS.path      <- file.path(WORKING, STATSFOLDER)
  # get country long name
  CountryLong  <- country.lookup[country.lookup[,1] == XXX, 2]
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", STATS.path))
  }
  
  
  if(! OldStyle) 
    return("new CSV not implemented")
  
  # Deaths and Exposures have f,m,b together
  deaths.name <- file.path(Rbin.path, "Deaths_1x1.Rdata")
  exposures.name <- file.path(Rbin.path,"Exposures_1x1.Rdata")
     ## CAB: instead of default load, probably best to load into tmp.env and then get from that
  load(deaths.name) # Deaths
  load(exposures.name) # Exposures
  Deaths$Agei <- as.numeric( gsub("([0-9])\\+", "\\1", Deaths$Age) )
  Exposures$Agei <- as.numeric( gsub("([0-9])\\+", "\\1", Exposures$Age) )

  ## process f,m,b variants.  For old-style .csv
  for( sex in c("f", "m", "b")){
    sexlong      <- ifelse( sex == "f", "Female", ifelse (sex == "m", "Male", "Total"))
    lt.name <- file.path(Rbin.path, paste0(sex, "ltper_1x1", ".Rdata"))
    csv.name <- file.path(STATS.path, paste0(sex, "ltper_1x1", ".csv"))
    ## CAB: Rbin saves of ltper etc. save to object 'output' even though there is a single object.
    ## this needs fixing as it creates lots of unneeded workaround code and is dangerous
    
    tmp.env <- new.env() # create a temporary environment; will disappear when function exits?
    load(lt.name, envir=tmp.env) # load workspace into temporary environment
    lt <- get("output", pos=tmp.env) # get the object from temp environment
    lt$Agei <- as.numeric( gsub("([0-9])\\+", "\\1", lt$Age) )
    
    Deaths.to.merge <- Deaths[, c("Year", "Agei", sexlong)]
    names(Deaths.to.merge)[3] <- "Deaths"
    Exposures.to.merge <- Exposures[, c("Year", "Agei", sexlong)]
    names(Exposures.to.merge)[3] <- "Exposures"
    lt.merge <- merge(lt, Deaths.to.merge, sort=FALSE)  # add Deaths, keyed on Year, Agei
    lt.merge <- merge(lt.merge, Exposures.to.merge, sort=FALSE) #add Exposures, keyed on Year, Agei
    lt.merge <- lt.merge[ order(lt.merge$Year, lt.merge$Agei), ]
    lt.merge$px <- 1.0 - lt.merge$qx
    lt.merge$Year <- as.numeric( lt.merge$Year)
    ## output of .csv file
    ## map column order into expected old form order, since no headers are output; FromYear = ToYear
    csv.cols <-c("Year", "Year", "Agei","Exposures", "Deaths", "mx", "ax", "qx", "px",  "lx", "dx", "Lx", 
                 "Tx", "ex" )
    csv.out <- lt.merge[, csv.cols]
    write.table(file = csv.name, x = csv.out, row.names = FALSE, col.names = FALSE, sep = ", " )
  }
 
}
