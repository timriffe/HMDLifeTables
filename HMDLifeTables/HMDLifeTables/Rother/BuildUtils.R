
# Author: triffe
###############################################################################

# this builds the package, incrementing the version automatically
BuildRLifeTablePackage <- function(MPVERSION = 6, MPVERSION.origin = "2013-01-01", commit = TRUE){
  
  # get time difference from origin
  increment <-  as.numeric(
    Sys.time() - 
      as.POSIXct(as.Date(MPVERSION.origin), 
        origin = ISOdatetime(1960,1,1,0,0,0),tz="PST")
  )
  
  # determine new version number:
  pkg.vs    <- paste0("Version: ",MPVERSION, ".", as.character(round(increment / 365.25, digits = 4)))
  
  # update version in DESCRIPTION file:
  DESC      <- readLines("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/DESCRIPTION")
  Version.i <- grep(DESC, pattern = "Version:")
  DESC[Version.i] <- pkg.vs
  writeLines(DESC, "/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/DESCRIPTION")
  
  if (commit){
    # commit package changes automatically:
    system(paste("cd /data/commons/triffe/COMMONS/git/HMD_Rlifetables_git/RLifeTables \n git commit -m ","'Package rebuild ",pkg.vs, "'"))
  }
  # build
  devtools::build(pkg = "/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables",
    path = "/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTablesBuilds")
}

# what's the newest build available?
NewestRLifeTablePackage <- function(build.folder = "/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTablesBuilds"){
  all.tar <- grep(list.files(build.folder),
    pattern = "tar.gz",value=TRUE)
  
  
  newest.build <- all.tar[which.max(gsub("[^[:digit:]]", "",  all.tar))]
  newest.build
}

# what version do I have installed?
InstalledRLifeTablePackage <- function(){
  paste0("RLifeTable_",installed.packages()["RLifeTable","Version"])
}
