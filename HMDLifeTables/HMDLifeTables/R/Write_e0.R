#' @title \code{Write_e0} a function to write formatted .txt files for period life expectancy values.
#'
#' @description This function requires that all pertinent \code{ltper_AxN.Rdata} data objects be present in the folder \code{WORKING/Rbin/} for males, females and both-sex tables. Objects are selected by grepping, which works for the present HMD file naming scheme. \code{.Rdata} files are read in, rounded, formatted and written back out into \code{WORKING/RSTATS/} by default, although this can be changed with the \code{STATSFOLDER} argument. This function must be called after ltper_AxN() (for all time intervals and sexes), which must have been run with the argument \code{save.bin = TRUE}.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' 
#' @return function called for its side effect of creating the files \code{E0per.txt}, \code{E0per_1x5.txt} and \code{E0per_1x10.txt}. No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export 
# Author: triffe
###############################################################################
Write_e0 <- function(
  WORKING = getwd(), 
  STATSFOLDER = "RSTATS",
  MPVERSION = 5,
  XXX = NULL){
  
  # define rounded function
  MatlabRoundFW <- function(x, digits = 0, pad = TRUE, Age = FALSE, totalL = 8){ 
    
    # this 1) rounds, and
    # 2) makes sure that the zeros stay on the end of the number string to the specified number of digits
    if (is.numeric(x)){
      fac   <- rep(10 ^ digits, length(x))
      x     <- sprintf(paste0("%.", digits, "f"), floor(x * fac + sign(x) * 0.5) / fac)
    }
    # strings will potentially vary in length due to integers: adds white space out to the
    # longest string (this is relevant e.g. for dx). ensures digit alignment.
    if (pad){
      maxw  <- max(nchar(x))  
      x     <- sprintf(paste0("%", ifelse(Age, {maxw - 1}, maxw), "s"), x)
      x     <- sprintf(paste0("%-", maxw, "s"), x) # this only had affect if Age == TRUE (dealing with '+')
    }
    # add optional left padding to specify total character width
    x       <- sprintf(paste0("%", totalL, "s"), x)
    x
  }
  # get country long name
  if (is.null(XXX)){
    XXX         <- ExtractXXXfromWORKING(WORKING) 
  }
  CountryLong   <- country.lookup[country.lookup[,1] == XXX, 2]
  
  # define, create stats path if necessary (should never be necessary)
  STATS.path    <- file.path(WORKING, STATSFOLDER)
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", STATS.path))
  }
  
  # Rbin folder must exist
  Rbin.folder   <- file.path(WORKING, "Rbin")
  if (!file.exists(Rbin.folder)){
    stop("no Rbin folder, function stopped")
  }
  # all files in folder
  files.2.write <- list.files(Rbin.folder)
  # all the files that have an 'lt' in them, should only be lifetables
  files.2.write <- files.2.write[grep(files.2.write, pattern = "lt")]
  # just 1x1, 1x5 and 1x10
  files.2.write <- files.2.write[grep(files.2.write, pattern = "1x")]
  coh.yes       <- any(grep(files.2.write, pattern = "coh"))
  # should be male female and 'both sex' version for each
  if (length(files.2.write) %% 3 != 0 | length(files.2.write) == 0){
    stop("in finding e0 estimates to write out, it looks like something is missing")
  }
  # some header info
  DateMod       <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ",")
  MPvers        <- ifelse(MPVERSION == 5, " MPv5 (May07)", "MPv6 (in development)\n")
  # year groups
  
  # begin loop over N
  for (N in c(1, 5, 10)){ #
    # we take e0 from the pertinent lifetables
    f.in.path   <- file.path(Rbin.folder, paste0("fltper_1x", N, ".Rdata"))
    m.in.path   <- file.path(Rbin.folder, paste0("mltper_1x", N, ".Rdata"))
    b.in.path   <- file.path(Rbin.folder, paste0("bltper_1x", N, ".Rdata"))
    
    # load females all the way in
    fem         <- local(get(load(f.in.path)))
    f.e0        <- with(fem, ex[Age == "0"])
    # just extract necessary info from others
    m.e0        <- with(local(get(load(m.in.path))), ex[Age == "0"])
    t.e0        <- with(local(get(load(b.in.path))), ex[Age == "0"])
    
    # pad year, if necessary
    Year        <- unique(fem$Year)
    if (max(nchar(Year)) == 4){
      Year      <- paste0("  ", Year, "   ")
    }
    # the name of the file to write out
    file.name      <- ifelse(N == 1, "E0per.txt", ifelse(N == 5, "E0per_1x5.txt", "E0per_1x10.txt"))
    write.out.file <- file.path(STATS.path, file.name)
    DataType       <- paste0(", Life expectancy at birth (period, 1x", N, ")")
   
    cat(
      # metadata header
      paste0(CountryLong, DataType, DateMod, MPvers,"\n"),
      # column headers, copied and pasted
      "  Year       Female    Male     Total",
      # the data, rounded and formatted in place- no tabbing
      paste(
         MatlabRoundFW(Year, totalL = 9),
         MatlabRoundFW(f.e0, digits = 2, totalL = 10),
         MatlabRoundFW(m.e0, digits = 2, totalL = 9),
         MatlabRoundFW(t.e0, digits = 2, totalL = 9),
         sep = ""), 
      file = write.out.file, sep = "\n")
    #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
    # --------------------------------------------------------------
    # lame copy and paste if cohort lifetables were run..
    if (coh.yes){
      
      DataType  <- paste0(", Life expectancy at birth (cohort, 1x", N, ")")
      
      f.in.path <- file.path(Rbin.folder, paste0("fltcoh_1x", N, ".Rdata"))
      m.in.path <- file.path(Rbin.folder, paste0("mltcoh_1x", N, ".Rdata"))
      b.in.path <- file.path(Rbin.folder, paste0("bltcoh_1x", N, ".Rdata"))
      
      fem       <- local(get(load(f.in.path)))
      f.e0      <- with(fem, ex[Age == "0"])
      m.e0      <- with(local(get(load(m.in.path))), ex[Age == "0"])
      t.e0      <- with(local(get(load(b.in.path))), ex[Age == "0"])
      
      Year      <- unique(fem$Year)
      if (max(nchar(Year)) == 4){
        Year    <- paste0("  ", Year, "   ")
      }
      file.name      <- ifelse(N == 1, "E0coh.txt", ifelse(N == 5, "E0coh_1x5.txt", "E0coh_1x10.txt"))
      write.out.file <- file.path(STATS.path, file.name)
      # write out
      cat(
        # metadata header
        paste0(CountryLong, DataType, DateMod, MPvers,"\n"),
        # column headers, copied and pasted
        "  Year       Female    Male     Total",
        # the data, rounded and formatted in place- no tabbing
        paste(
          MatlabRoundFW(Year, totalL = 9),
          MatlabRoundFW(f.e0, digits = 2, totalL = 10),
          MatlabRoundFW(m.e0, digits = 2, totalL = 9),
          MatlabRoundFW(t.e0, digits = 2, totalL = 9),
          sep = ""), 
        file = write.out.file, sep = "\n")
      #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", write.out.file))
    } # end cohort 'if' statement
  } # end N loop
} # end function definition
