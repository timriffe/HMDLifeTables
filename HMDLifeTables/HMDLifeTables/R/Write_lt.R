#' @title \code{Write_lt} a function to write formatted .txt files for both period and cohort lifetables.
#'
#' @description This function requires that all pertinent \code{ltper_AxN.Rdata} and \code{ltcoh_AxN()} data objects be present in the folder \code{WORKING/Rbin/} for males, females and both-sex tables. Objects are selected by grepping, which works for the present HMD file naming scheme. \code{.Rdata} files are read in, rounded, formatted and written back out into \code{WORKING/RSTATS/} by default, although this can be changed with the \code{STATSFOLDER} argument. This function must be called after ltper_AxN() (for all time intervals and sexes), which must have been run with the argument \code{save.bin = TRUE}.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param CountryLong the HMD country full name.
#' @return function called for its side effect of creating the lifetable txt output files, e.g. \code{mltper_1x1.txt} and other time/sex configurations. No value returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export 
# Author: triffe
###############################################################################
Write_lt <- function(
  WORKING = getwd(), 
  STATSFOLDER = "RSTATS", 
  MPVERSION , # explicit, no default
  XXX = NULL,
  CountryLong = NULL){
  
  # MatlabRound() is for rounding output, should give same result as matlab, assuming that's important
  # by CB, updated by TR to take digits as arg.
  MatlabRoundFW <- function(x, digits = 0, pad = TRUE, Age = FALSE, totalL = 8){ 
    NAsmatch <- paste0( substr("                   ", 1, totalL - 2), "NA" )
    NAsreplace <- paste0( substr("                   ", 1, totalL - 1), "." )
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
    x         <- ifelse( x == NAsmatch, NAsreplace, x)  # replace string NAs with string "."

    return(x)
  }
  # get country abbrev
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  # for the metadata header: country long name
  
  if(length(CountryLong) == 0){
    warning("*** !!! Missing long country name; output will be affected")
  }
  # ltper_AxN() should have been run with save.bin = TRUE, 
  # such that there's an Rbin folder waiting there with stuff in it
  Rbin.path       <- file.path(WORKING, "Rbin")
  # save formatted .txt out to this folder, make sure exists
  STATS.path      <- file.path(WORKING, STATSFOLDER)
  
  
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", STATS.path))
  }
  # names of all the files in Rbin
  files.2.write   <- list.files(Rbin.path)
  # all the files that have an 'lt' in them, should only be lifetables
  files.2.write   <- files.2.write[grep(files.2.write, pattern = "lt")]
  # read in each file one at a time, and repeat
  for (this.file in files.2.write){
    # Rdata files can be read in with load()- each file is simply called 'output'
    load(file.path(Rbin.path, this.file))
    # will have same name, but with txt suffix
    write.out.file <- file.path(STATS.path, gsub(this.file, pattern = "Rdata", replacement = "txt"))
    # this gets the '1x1', etc
    dims         <- gsub(gsub(this.file, pattern = ".Rdata", replacement = ""),pattern = ".*_",replacement="") 
    # either period or cohort
    perCond      <- any(grep(this.file, pattern = "per"))
    DataType     <- paste0(", Life tables (", ifelse(perCond, "period ", "cohort "),
                    dims, "), ")
    # if the filename contains an 'f' it's for females..
    magic.letter <- unlist(strsplit(this.file,""))[1]
    sexlong      <- ifelse( magic.letter == "f", "Females", ifelse (magic.letter == "m", "Males", "Total"))
    
    # time stamp (perhaps make more precise sometime)
    DateMod      <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ";")
    # Methods Protocol version
    MPvers       <- ifelse(MPVERSION == 5, " MPv5 (May07)", "  Methods Protocol: v6 (2017)")
    
    # this is fancy character padding. 'Year' is either 4 or 9 characters long- 
    # make 9, spaced properly
    if (max(nchar(output$Year)) == 4){
      output$Year <- paste0("  ", output$Year, "   ")
    }    
    # same for Age- this is harder to deal with- this
    # works in combo with the ifelse() statements below 
    AgeCond      <- max(nchar(output$Age)) == 4
    if (AgeCond){
      output$Age <- paste0(output$Age, "    ")
    }  
      
    # replace NAs, NaNs with ".": this appears to only come up with cohort tables, since
    # period tables have smoothed final values
    # CAB: both Nan and NA can occur so replace both
    mx <- MatlabRoundFW(output$mx, digits = 5, totalL = ifelse(AgeCond, 8, 11))
    mx <- gsub(mx, pattern = "NA|NaN", replacement = " .")
    
    qx <- MatlabRoundFW(output$qx, digits = 5, totalL = 9)
    qx <- gsub(qx, pattern = "NA|NaN", replacement = " .")
    
    ax <- MatlabRoundFW(output$ax, digits = 2, totalL = 6)
    ax <- gsub(ax, pattern = "NA|NaN", replacement = " .")
    
    ex <- MatlabRoundFW(output$ex, digits = 2, totalL = 7)
    ex <- gsub(ex, pattern = "NA|NaN", replacement = "  .")
    # begin writing out
    cat(
      # metadata header
      paste0(CountryLong, DataType, sexlong, DateMod, MPvers,"\n"),
      # column headers, copied and pasted
      "  Year          Age         mx       qx    ax      lx      dx      Lx       Tx     ex",
      # the data, rounded and formatted in place- no tabbing
      paste(MatlabRoundFW(output$Year, totalL = 9), 
        MatlabRoundFW(output$Age, Age = TRUE, totalL = ifelse(AgeCond, 14, 11)), 
        mx,
        qx,
        ax,
        MatlabRoundFW(output$lx),
        MatlabRoundFW(output$dx),
        MatlabRoundFW(output$Lx),
        MatlabRoundFW(output$Tx, totalL = 9),
        ex, 
        sep = ""), 
      file = write.out.file, sep = "\n")
    #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ",write.out.file))
  } # end file loop
} # end function definition

