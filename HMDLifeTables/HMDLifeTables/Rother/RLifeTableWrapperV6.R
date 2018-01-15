## R Start file
options(showErrorCalls=TRUE)
options(showWarnCalls=TRUE)
options(show.error.locations="top")
options(warn=1)  # immediately print warnings

cat("HMD R Lifetable Processing....\n")


cat("Loading devtools and RLifeTables package \n")

library(roxygen2)
library(devtools)

load_all(pkg='/data/commons/boe/HMDLifeTables_R_cab.git/HMDLifeTables/HMDLifeTables')

#cat(paste("Loading...",file.path(Sys.getenv("PROJ"),"RunHMDCountry.R"),"\n"))

#source(file.path(Sys.getenv("PROJ"),"RunHMDCountry.R"), echo=FALSE)
       
cat("Pulling environment variables that may override defaults...\n")

## construct parameter list.  NB  parms$field = NULL will not stick, while list(field=NULL) will
## so use list()
parms <- list(
    WORKING = getwd(),
    OPENAGE = ifelse( nzchar(Sys.getenv("OPENAGE")), as.integer(Sys.getenv("OPENAGE")), 110),
    RADIX =  ifelse( nzchar(Sys.getenv("RADIX")), as.integer(Sys.getenv("RADIX")), 100000),
    CAGEEXTRP = ifelse( nzchar(Sys.getenv("CAGEEXTRP")), as.integer(Sys.getenv("CAGEEXTRP")), 90),
    MPVERSION =  ifelse( nzchar(Sys.getenv("MPVER")), as.integer(Sys.getenv("MPVER")), 6),
    STATSFOLDER = "STATS",
    XXX = ifelse( nzchar(Sys.getenv("XXX")), Sys.getenv("XXX"), 
                  ifelse(exists("XXX"), XXX, NULL) ),
    CountryName = ifelse( nzchar(Sys.getenv("CountryName")), Sys.getenv("CountryName"),
                  ifelse(exists("CountryName"), CountryName, NULL) ),
    IDBPATH = NULL,
    LDBPATH = NULL,
    SAVEBIN = TRUE
    );
    
if(  nzchar(Sys.getenv("LDBPATH")) )
    parms["LDBPATH"] <- Sys.getenv("LDBPATH")
##

if(  nzchar(Sys.getenv("INDBPATH")) )
    parms["IDBPATH"] <- Sys.getenv("INDBPATH")
##
cat("Running the lifetables...\n")

print(parms)
## run the top-level script from the RLifeTables package, which suitable parameter settings
RunHMDCountry(WORKING=parms$WORKING, OPENAGE=parms$OPENAGE, RADIX=parms$RADIX, CAGEEXTRP=parms$CAGEEXTRP,
              MPVERSION=parms$MPVERSION, STATSFOLDER=parms$STATSFOLDER, XXX=parms$XXX, CountryName=parms$CountryName,
              LDBPATH=parms$LDBPATH,
              IDBPATH=parms$IDBPATH, SAVEBIN=parms$SAVEBIN)


cat("HMD R Lifetable ... finished")
