#' @title \code{Write_Deaths_lexis} a function to prepare and save the file \code{deaths_lexis.txt}.
#'
#' @description This function is essentially a stand-alone that prepares \code{deaths_lexis.txt} straight from LexisDB files. At this time is makes no indispensable calls to other functions, although it could be made to use \code{getPeriodComponents} as other functions do, but this is not a high priority.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param OPENAGE the desired open age. Default value is 110.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param CountryLong the HMD country full name.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' 
#' @return function called for its side effect of creating the files \code{Population.txt} or \code{Population5.txt}. No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export
# Author: triffe
###############################################################################
Write_Deaths_lexis <- function(
  WORKING = getwd(),
  STATSFOLDER = "RSTATS",
  OPENAGE = 110,
  XXX = NULL,
  CountryLong = NULL,
  LDBPATH = NULL,
  MPVERSION # explicit, no default
  ){
  
  # MatlabRound() is for rounding output, should give same result as matlab, assuming that's important
  # by CB, updated by TR to take digits as arg.
  
  MatlabRoundFW <- function(x, digits = 0, pad = TRUE, Age = FALSE, totalL = 8){ 
    
    # this 1) rounds, and
    # 2) makes sure that the zeros stay on the end of the number string to the specified number of digits
    NAsmatch <- paste0( substr("                   ", 1, totalL - 2), "NA" )
    NAsreplace <- paste0( substr("                   ", 1, totalL - 1), "." )
    if (is.numeric(x)){
      fac     <- rep(10 ^ digits, length(x))
      x       <- sprintf(paste0("%.", digits, "f"), floor(x * fac + sign(x) * 0.5) / fac)
    }
    # strings will potentially vary in length due to integers: adds white space out to the
    # longest string (this is relevant e.g. for dx). ensures digit alignment.
    if (pad){
      maxw    <- max(nchar(x))  
      x       <- sprintf(paste0("%", ifelse(Age, {maxw - 1}, maxw), "s"), x)
      x       <- sprintf(paste0("%-", maxw, "s"), x) # this only had affect if Age == TRUE (dealing with '+')
    }
    # add optional left padding to specify total character width
    x         <- sprintf(paste0("%", totalL, "s"), x)
    x         <- ifelse( x == NAsmatch, NAsreplace, x)  # replace string NAs with string "."
    return(x)

  }
  
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  
  # for the metadata header: country long name
  if(length(CountryLong) == 0){
    warning("*** !!! Missing long country name; output will be affected")
  }
  
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
#--------------------------------------------------------
# read in data:
# LDB 
  ldb.path.m  <- file.path(LDBPATH, paste0("m", XXX, ".txt"))
  LDBobj.m    <- read.table(ldb.path.m, header = FALSE, sep = ",", 
                   col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  ldb.path.f  <- file.path(LDBPATH, paste0("f", XXX, ".txt"))
  LDBobj.f    <- read.table(ldb.path.f, header = FALSE, sep = ",", 
                   col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  
  # short vec of years
  yr          <- unique( LDBobj.m$Year)
  yr          <- yr[-length(yr)]
  # remove last year, only there for pop
  LDBobj.m    <- LDBobj.m[LDBobj.m$Year %in% yr, ]
  LDBobj.f    <- LDBobj.f[LDBobj.f$Year %in% yr, ]

  Cohort      <- LDBobj.m$Cohort
  
  Male        <- LDBobj.m$Deaths
  Female      <- LDBobj.f$Deaths
  
  ## CAB: fix BEL which has NAs for some years, e.g. 1917 and where LDB has -1 for Population, Deaths
  Male <- ifelse( Male == -1, NA, Male)
  Female <- ifelse( Female == -1, NA, Female)
  
  # --------------------------
  # silly thing to sum open age group:
  Male        <- matrix(Male, ncol = length(yr), dimnames = list(rep(0:130, each = 2), yr))
  Female      <- matrix(Female, ncol = length(yr), dimnames = list(rep(0:130, each = 2), yr))
  
  i.OPENAGE   <- OPENAGE * 2 + 1
  Male[i.OPENAGE, ]   <- colSums(Male[i.OPENAGE : nrow(Male), ])
  Female[i.OPENAGE, ] <- colSums(Female[i.OPENAGE : nrow(Female), ])

  Male        <- Male[1:i.OPENAGE, ]
  Female      <- Female[1:i.OPENAGE, ]
  Total       <- Male + Female
  # --------------------------
# labeling is tricky
  years       <- as.integer(rep(yr, each = i.OPENAGE))
  ages.offset <- rep(as.integer(rownames(Male)) + c(rep(c(0,1),OPENAGE),0), length(yr))
  cohorts     <- years - ages.offset 
  age         <- c(rep(0:(OPENAGE - 1), each = 2), paste0(OPENAGE, "+"))  # CAB, cast as string, then back cast in Matlab to number, etc.
  cohorts[age == paste0(OPENAGE, "+")] <- "   ."

  # set up paths
  STATS.path     <- file.path(WORKING, STATSFOLDER)
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #system(paste0("chgrp hmdcalc ", STATS.path))
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
  }
  write.out.file <- file.path(STATS.path, "Deaths_lexis.txt")
  
  DateMod        <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ",")
  # Methods Protocol version
  MPvers         <- ifelse(MPVERSION == 5, " MPv5 (May07)", " MPv6 (Nov17)\n")
  DataType       <- ",  Deaths (Lexis triangle)"
  
  # write it out!
  cat(
    # metadata header
    paste0(CountryLong, DataType, DateMod, MPvers,"\n"),
    # column headers, copied and pasted
    " Year   Age   Cohort      Female        Male       Total",
    # the data, rounded and formatted in place- no tabbing
    paste(MatlabRoundFW(years, totalL = 5), 
      MatlabRoundFW(age, Age = TRUE, totalL = 6),
      MatlabRoundFW(cohorts, totalL = 9),
      MatlabRoundFW(as.vector(Female), digits = 2, totalL = 12),
      MatlabRoundFW(as.vector(Male), digits = 2, totalL = 12),
      MatlabRoundFW(as.vector(Total), digits = 2, totalL = 12), 
      sep = ""), 
    file = write.out.file, sep = "\n")
  #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
  #system(paste0("chgrp hmdcalc ", write.out.file))
  
}
