#' @title \code{Write_Births} a function to prepare, format and write \code{Births.txt}.
#'
#' @description This function pulls birth count data from the Lexis Database (age = 0, lexis = 1) and sets up the standard output with columns \code{Year}, \code{Female}, \code{Male}, and \code{Total}. This function makes no consequential calls to other functions, and is not called by any other function, nor does it require another function to have been run previously.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param PVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param CountryLong the HMD country full name.
#' 
#' @return function called for its side effect of creating the file \code{Births.txt}. No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @import reshape2 
#' @export
# Author: triffe
###############################################################################
Write_Births <- function(
  WORKING = getwd(),
  STATSFOLDER = "RSTATS",
  LDBPATH = NULL,
  IDBPATH = NULL,
  MPVERSION , # explicit, no default
  XXX = NULL,
  CountryLong = NULL){
  
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  
  # for the metadata header: country long name
  
  if(length(CountryLong) == 0){
    warning("*** !!! Missing long country name; output will be affected")
  }
  
  if (is.null(LDBPATH)){
    LDBPATH       <- file.path(WORKING, "LexisDB")
  }
    if (is.null(IDBPATH)){
    IDBPATH       <- file.path(WORKING, "InputDB")
  }
  # Lexis DB object
  ldb.path.f      <- file.path(LDBPATH, paste0("f", XXX, ".txt"))
  LDBobj.f        <- read.table(ldb.path.f, header = FALSE, sep = ",", 
                        col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  ldb.path.m      <- file.path(LDBPATH, paste0("m", XXX, ".txt"))
  LDBobj.m        <- read.table(ldb.path.m, header = FALSE, sep = ",", 
                        col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  
  Bf.ldb              <- with(LDBobj.f, Population[Age == 0 & Lexis == 1 & Population > -1])
  Bm.ldb              <- with(LDBobj.m, Population[Age == 0 & Lexis == 1 & Population > -1])
  Bt.ldb              <- Bf.ldb + Bm.ldb
  Year.ldb            <- with(LDBobj.m, Year[Age == 0 & Lexis == 1 & Population > -1])


  ## CAB: override the above and take Births from InputDB, since that will give the full span of
  ## available birth data.  
  ## births from InputDB XXXbirths.txt file
  ## TODO: reconstruct LexisDB have "wide" structure, containing full span of data and BOY, EOY values (implicit
  ## territorial adjustment factors)
  
  idb.path      <- file.path(IDBPATH, paste0(XXX, "birth.txt"))
  idb.births <- read.csv(idb.path, header=TRUE)

  ## working subset
  idb.births <- idb.births[, c("Sex", "Year", "Births", "Access","LDB")]
  idb.births <- idb.births[ idb.births$LDB == 1, ]
  # for some aggregates, the Births file consists of a Header only or a Header and a line of missing values
  # in which case we skip IDB
  
  if( length(unique(idb.births$Year) ) <= 1 && unique(idb.births$Year) == "."){ # use LDB as Births source
    isEmptyIDBBirths <- TRUE
    warning("*** Empty InputDB Births file -- using LexisDB for Births")
    Bf <- Bf.ldb
    Bm <- Bm.ldb
    Bt <- Bt.ldb
    Year <- Year.ldb
    
    
  } else {  # use IDB as Births source
    
    isEmptyIDBBirths <- FALSE
    ## set year as factor with complete recorded year range, to impute any missing entries 
    idbbirthsExpectedYears <- seq(from=min(idb.births$Year), to=max(idb.births$Year), by=1)
    
    idb.births.m <- melt(idb.births, id.vars=c("Sex", "Year", "Access", "LDB") )
    idb.births.c <- dcast(idb.births.m, Year + Sex + Access + LDB  ~ variable , drop=FALSE, fill=NULL)
    
    ## remove duplicates where there are identical Year,Sex, Access entries   
    ## Warn here because this situation is very confusing
    idb.births.dups <- duplicated(idb.births.c[,1:3]) 
    if( any(idb.births.dups) ){
      warning(paste("*** ", XXX, ": There are duplicate/conflicting entries for Births; preferring cases where LDB==1"))
      idb.births.c <- idb.births.c[ !idb.births.dups,]
    }
    Bf <- idb.births.c$Births[ idb.births.c$Sex=='f' ]
    Bm <- idb.births.c$Births[ idb.births.c$Sex=='m' ]
    Year <-  as.character( idb.births.c$Year[ idb.births.c$Sex=='f' ] ) #was factor
    Year <- as.integer(Year)
    if( length(Bf) != length(Bm) || length(Bf) != length(Year) ){
      warning(paste("*** ", XXX, ": Unequal Male, Female births lengths in IDB") )
    }
    
    if( length( c( setdiff(Year, Year.ldb), setdiff(Year.ldb, Year) )) > 0 ){
      warning(paste("*** ", XXX, ": Different Birth Years between IDB, LDB:",  paste(setdiff(Year, Year.ldb), setdiff(Year.ldb, Year) , collapse=" ")) )
    }
    if( length(setdiff(idbbirthsExpectedYears, Year)) > 0 ){
      warning( paste("*** Fewer Birth years than expected in IDB -- missing data for years", setdiff(idbbirthsExpectedYears, Year), collapse=" ") )
    }
    Bt <- Bf + Bm  
    
  }
  
  
  DateMod         <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ",")
  # Methods Protocol version
  MPvers          <- ifelse(MPVERSION == 5, " MPv5 (May07)", "MPv6 (Nov17)\n")
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
      sprintf(paste0("%-", 10, "i"), Year),
      ifelse(is.na(Bf), "         .", sprintf(paste0("%", 10, "i"), Bf) ),
      ifelse(is.na(Bm), "         .", sprintf(paste0("%", 10, "i"), Bm) ),
      ifelse(is.na(Bt), "         .", sprintf(paste0("%", 10, "i"), Bt) ),
      sep = ""), 
    file = write.out.file, sep = "\n")
  #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
  #system(paste0("chgrp hmdcalc ", write.out.file));

} #end of function Write_Births()

  
  ## idb.births.out <- data.frame(Year = Years, Female = Bf,  Male = Bm, Total = Bt)
  ## fmt.births.out <- "%10s%10i%10i%10i"
  ## fmt.births.header <- "%10s%10s%10s%10s"
  ## idb.births.header <- sprintf(fmt.births.header, c("Year", "Female", "Male", "Total"))
  
  ## ## http://stackoverflow.com/questions/28058501/split-data-frame-for-passing-to-sprintf-in-r#28059117
  ## writeData <- function(DataSet,FirstLine, FmtString,fName){
  ##   outLines <- do.call("sprintf", c(FmtString, DataSet))
  ##   if (nchar(FirstLine) > 0){
  ##     writeLines(FirstLine, fName)
  ##   }
    
  ##   writeLines(outLines,fName)
  ##   return(0)
  ## }


  ## }


## mb<-data.frame(Year=c(1900,1900, 1901, 1903, 1903),
##                Sex= c("M","F","F","M","F"),
##                Births=c(100,101,102,103,104))
## year.f <- factor(mb$Year, levels=min(mb$Year):max(mb$Year))
## sex.f  <- factor(mb$Sex, levels=c("M","F"))
## mb$Year <- year.f
## mb$Sex  <- sex.f
## mb.m <- melt(mb, id.vars=c("Sex", "Year" ) )
## mb.c <- dcast(mb.m, Year + Sex ~ variable , drop=FALSE,fill=NULL)

