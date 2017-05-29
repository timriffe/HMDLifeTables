#' @title \code{RunHMDCountry} a function to calculate, prepare and write out a full \code{STATS} folder given up-to-date \code{InputDB} and \code{LexisDB} folders.
#' 
#' @description This function is the highest-level function, calling all other lifetable functions in the tree either directly or indirectly. Some defaults are not controllable via arguments, though these may be incorporated as this function comes into use, especially those arguments that control whether or not cohort lifetables are to be run- right now it's set to automatic detection. This function can take around 30 seconds to run completely, though this may get faster in future edits.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param OPENAGE the desired open age. Default value is 110.
#' @param RADIX the lifetable radix (l_0). Default 1e5.
#' @param CAGEEXTRP earliest age at which almost-extinct cohorts may begin extrapolation. Default 90
#' @param MPVERSION 5 or 6. See individual functions for descriptions of differences
#' @param STATSFOLDER the folder name where output is to be written to (not a full path). Default \code{"RSTATS"}.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' @param SAVEBIN logical. Default \code{TRUE}.  Should interim \code{R} output be saved as e.g., \code{Rbin/ltper_1x1.Rdata} as well? Appropriate name is derived systematically. In this case objects are saved separately. This option can speed calculations.
#' @param reproduce.matlab logical. Default \code{FALSE}. Should we reproduce all aspects of the matlab code?  
#'
#' @details reproduce.matlab only has effect in ltperBoth_AxN.R  All implementations of Matlab lifetable code had discrepancies with 
#' V4, V5 of the HMD Methods Protocol.  When \code{TRUE}, then the discrepancies will be reproduced.  Affects V5 and V6.
#' @return Running this function has many side-effects, most notably, filling up a folder \code{WORKING/Rbin/} with unrounded R binary versions of all output objects, and then filling up a folder \code{WORKING/RSTATS} (name can be changed) with final .txt formatted output. 
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export

RunHMDCountry <- function(WORKING = getwd(), 
  OPENAGE = 110,
  RADIX = 1e5,
  CAGEEXTRP = 90,
  MPVERSION = 5,
  STATSFOLDER = "RSTATS",
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL,
  SAVEBIN = TRUE,
  reproduce.matlab=FALSE){
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH       <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH       <- file.path(WORKING, "InputDB")
  }
  # fresh start XXX <- "SWE"
  CleanRstuff(WORKING = WORKING)
  # country folder:
  # building blocks: feed in via memory rather than loaded
  females.p <- getPeriodComponents(WORKING = WORKING, 
                                   sex = "f", 
                                   OPENAGE = 110, 
                                   save.bin = SAVEBIN, 
                                   XXX = XXX, 
                                   LDBPATH = LDBPATH, 
                                   IDBPATH = IDBPATH)
  males.p   <- getPeriodComponents(WORKING = WORKING, 
                                   sex = "m", 
                                   OPENAGE = 110, 
                                   save.bin = SAVEBIN, 
                                   XXX = XXX, 
                                   LDBPATH = LDBPATH, 
                                   IDBPATH = IDBPATH)
  females.c <- getCohortComponents(WORKING = WORKING, 
                                   sex = "f", 
                                   save.bin = SAVEBIN, 
                                   XXX = XXX, 
                                   LDBPATH = LDBPATH, 
                                   IDBPATH = IDBPATH)
  males.c   <- getCohortComponents(WORKING = WORKING, 
                                   sex = "m", 
                                   save.bin = SAVEBIN, 
                                   XXX = XXX, 
                                   LDBPATH = LDBPATH, 
                                   IDBPATH = IDBPATH)
 # ensure fresh check on cohorts
 checkRuncoh(WORKING = WORKING, 
              males = males.c,
              females = females.c,
              OPENAGE = OPENAGE,
              save.bin = SAVEBIN, 
              XXX = XXX,
              LDBPATH = LDBPATH,
              IDBPATH = IDBPATH)                               
  # begin lifetables N=5;.abridged = TRUE
  for (.abridged in c(FALSE, TRUE)){
    for (N in c(1, 5, 10)){
      # --------------------------------------------------------- 
      # MALES
        # period tables, single sex 
        test <- try(ltper_AxN(
                          WORKING = WORKING,
                          perComp = males.p,
                          sex = "m", 
                          OPENAGE = OPENAGE, 
                          RADIX = RADIX, 
                          N = N, 
                          abridged = .abridged, 
                          MPVERSION = MPVERSION, 
                          save.bin = SAVEBIN, 
                          XXX = XXX, 
                          LDBPATH = LDBPATH, 
                          IDBPATH = IDBPATH))
        if (class(test) == "try-error"){
          stop("problem in ltper_AxN()", "m", paste0(ifelse(.abridged, 5, 1), "x", N))
        }
        # cohort tables, single sex (checked live if need to be run)
        test <-  try(ltcoh_AxN(
                          WORKING = WORKING,
                          cohComp = males.c,
                          sex = "m", 
                          OPENAGE = OPENAGE, 
                          RADIX = RADIX, 
                          CAGEEXTRP = CAGEEXTRP,
                          N = N, 
                          abridged = .abridged, 
                          save.bin = SAVEBIN, 
                          run.condition = "both", 
                          warn.if.not.run = FALSE, 
                          XXX = XXX, 
                          LDBPATH = LDBPATH, 
                          IDBPATH = IDBPATH))
        if (class(test) == "try-error"){
          stop("problem in ltcoh_AxN()", "m", paste0(ifelse(.abridged, 5, 1), "x", N))
        }
       # FEMALES
       # period tables, single sex females
        test <- try(ltper_AxN(
                          WORKING = WORKING,
                          perComp = females.p, 
                          sex = "f", 
                          OPENAGE = OPENAGE, 
                          RADIX = RADIX, 
                          N = N, 
                          abridged = .abridged, 
                          MPVERSION = MPVERSION, 
                          save.bin = SAVEBIN, 
                          XXX = XXX, 
                          LDBPATH = LDBPATH, 
                          IDBPATH = IDBPATH))
        if (class(test) == "try-error"){
          stop("problem in ltper_AxN()", "f", paste0(ifelse(.abridged, 5, 1), "x", N))
        }
        # cohort tables, single sex (checked live if need to be run)
        test <-  try(ltcoh_AxN(
                          WORKING = WORKING,
                          cohComp = females.c, 
                          sex = "f", 
                          OPENAGE = OPENAGE, 
                          RADIX = RADIX, 
                          CAGEEXTRP = CAGEEXTRP,
                          N = N, 
                          abridged = .abridged, 
                          save.bin = SAVEBIN, 
                          run.condition = "both", 
                          warn.if.not.run = FALSE, 
                          XXX = XXX, 
                          LDBPATH = LDBPATH, 
                          IDBPATH = IDBPATH))
        if (class(test) == "try-error"){
          stop("problem in ltcoh_AxN()", "f", paste0(ifelse(.abridged, 5, 1), "x", N))
        }
       # end innermost male-female loop
      # the non sex-specific functions:
      # both sex period tables, need to be done once both single sex tables are done
      test <-  try(ltperBoth_AxN(
                        WORKING = WORKING, 
                        males = males.p, 
                        females = females.p, 
                        OPENAGE = OPENAGE, 
                        RADIX = RADIX, 
                        N = N, 
                        abridged = .abridged, 
                        MPVERSION = MPVERSION, 
                        save.bin = SAVEBIN, 
                        XXX = XXX, 
                        LDBPATH = LDBPATH, 
                        IDBPATH = IDBPATH,
                        reproduce.matlab=reproduce.matlab ))
      if (class(test) == "try-error"){
        stop("problem in ltperBoth_AxN()", paste0(ifelse(.abridged, 5, 1), "x", N))
      }
      # both sex cohort tables
      test <-  try(ltcohBoth_AxN(
                        WORKING = WORKING, 
                        males = males.c, 
                        females = females.c, 
                        OPENAGE = OPENAGE, 
                        CAGEEXTRP = CAGEEXTRP, 
                        RADIX = RADIX, 
                        N = N, 
                        abridged = .abridged,     
                        save.bin = SAVEBIN, 
                        run.condition = "both", 
                        warn.if.not.run = FALSE, 
                        XXX = XXX, 
                        LDBPATH = LDBPATH, 
                        IDBPATH = IDBPATH))
      if (class(test) == "try-error"){
        stop("problem in ltcohBoth_AxN()", paste0(ifelse(.abridged, 5, 1), "x", N))
      }
      # period Deaths, Exposures, Mx
      test <-  try(Exposures_Deaths_Mx_AxN(
                        WORKING = WORKING, 
                        males = males.p, 
                        females = females.p,
                        abridged = .abridged, 
                        N = N,
                        OPENAGE = OPENAGE, 
                        MPVERSION = MPVERSION,
                        save.bin = SAVEBIN, 
                        XXX = XXX,
                        LDBPATH = LDBPATH, 
                        IDBPATH = IDBPATH))
      if (class(test) == "try-error"){
        stop("problem in Exposures_Deaths_Mx_AxN()", paste0(ifelse(.abridged, 5, 1), "x", N))
      }
      # cohort Exposures, Mx
      test <- try(cExposures_and_cMx_AxN(
                        WORKING = WORKING, 
                        males = males.c, 
                        females = females.c,
                        abridged = .abridged, N = N,
                        OPENAGE = OPENAGE, 
                        min.years = 30, 
                        MPVERSION = MPVERSION,
                        save.bin = SAVEBIN, 
                        XXX = XXX, 
                        LDBPATH = LDBPATH, 
                        IDBPATH = IDBPATH))
      if (class(test) == "try-error"){
        stop("problem in cExposures_and_cMx_AxN()", paste0(ifelse(.abridged, 5, 1), "x", N))
      }
    } # end N loop
    test <- try(Population_A(
                        WORKING = WORKING, 
                        abridged = .abridged, 
                        OPENAGE = OPENAGE, 
                        save.bin = SAVEBIN, 
                        XXX = XXX, 
                        LDBPATH = LDBPATH, 
                        IDBPATH = IDBPATH))
    if (class(test) == "try-error"){
      stop("problem in Population_A()", ifelse(.abridged, 5, 1))
    }
  } # end abridge loop
  # write out Rbin files to txt, formatted
  Write_lt(                 WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            XXX = XXX, 
                            MPVERSION = MPVERSION)
  Write_Exposures_Deaths_Mx(WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            XXX = XXX, 
                            MPVERSION = MPVERSION)
  Write_Births(             WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            LDBPATH = LDBPATH,
                            IDBPATH = IDBPATH,
                            XXX = XXX, 
                            MPVERSION = MPVERSION)
  Write_e0(                 WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            XXX = XXX, 
                            MPVERSION = MPVERSION)
  Write_Deaths_lexis(       WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            OPENAGE = OPENAGE,
                            XXX = XXX, 
                            LDBPATH = LDBPATH,
                            MPVERSION = MPVERSION)
  
  Write_csv(                WORKING = WORKING, 
                            STATSFOLDER = STATSFOLDER, 
                            MPVERSION = MPVERSION,
                            OldStyle = TRUE,
                            XXX = XXX 
                            )  
  
 
} # end function definition

###############################################################################
# Package Documentation:
#' The full suite of lifetable functions used by the HMD
#'
#' RLifeTables provides the entirety of post LexisDB functionality required
#' by the HMD for country updates and testing. There are four kinds of functions 
#' in this package: 1) top-level functions, used for creating in-memory R \code{data.frame}
#' output for the various time and sex configurations of period and cohort lifetables.
#' These include \code{ltper_AxN()}, \code{ltcoh_AxN()}, their both-sex variants, and the omnibus 
#' \code{RunHMDCountry()}. 2) auxiliary functions that do minor things. These have been written 
#' in the name of modularity and include functions \code{Abrdige()}, \code{YearAgg()}, code{CDa0()}
#' and others. These can be useful, but need not be handled at the top level for a mere
#' country update. Top level functions call these as necessary. 3) data preparation functions
#' such, \code{getCohortComponents()} and \code{getPeriodComponents()}, which feed data into 
#' downstream functions in a standard and useful format (territorial adjustments in-place) and 
#' 4) output formatting functions for preparing the .txt output files in a STATS folder (here 
#' RSTATS, unless otherwise specified by the user).
#'
#' Generating lifetables.
#' For a country update, try: (if you're at UCB)
#' @examples \dontrun{
#'    RunHMDCountry(WORKING = "/data/commons/hmd/HMDWORK/ISL") # Does Iceland 
#' } 
#' Assuming the standard HMD folder structure, LexisDB and InputDB file formats, the 
#' entire body of HMD output will be generated. Intermediate files (unrounded) are also 
#' produced as R binary (\code{.Rdata}) files in the \code{WORKING/Rbin/} folder. These 
#' can be useful either as an alternative final product or for testing. \code{RunHMDCountry()} 
#' fully erases the specified \code{STATSFOLDER} and \code{Rbin} folder to ensure that 
#' it is always working with the most recent files. In this way, each run is a full and 
#' fresh run.
#' 
#' While this is the only necessary function for typical country updates, all other functions 
#' are directly available to users, such that one could dynamically change arguments to produce 
#' different output. Imagine for instance that you wanted cohort lifetables for a wider range of 
#' years than the default; you could specify:
#' @examples \dontrun{
#'    ISLfclt <- ltcoh_AxN(WORKING = "/data/commons/hmd/HMDWORK/ISL",
#'                         sex = "f",
#'                         OPENAGE = 95,   # default 110
#'                         CAGEEXTRP = 80, # default 90
#'                         save.bin = FALSE,
#'                         MPVERSION = 5) 
#' } 
#' This would give something like 25 more cohorts than the default, albeit under stronger assumptions.
#' One could also avoid reference to a hard working directory by using LexisDB data in-memory:
#' In this case, one would use the \code{cohComp} argument, which refers to the standard output
#' of \code{getCohortComponents()}
#'  @examples \dontrun{
#' ISLcohComp <- getCohortComponents(
#'                         WORKING = "/data/commons/hmd/HMDWORK/ISL", 
#'                         OPENAGE = 95, 
#'                         sex = "f")
#' head(ISLcohComp) # if you had data from some non-HMD source or test-source
#'                  # and wanted to use HMD lifetables on it, then it must be
#'                  # in a \code{data.frame} whose structure mirrors this one.
#' # lifetable using in-memory data:
#' ISLfclt <- ltcoh_AxN(cphComp = ISLcohComp,
#'                         sex = "f",
#'                         OPENAGE = 95,
#'                         CAGEEXTRP = 80,
#'                         save.bin = FALSE,
#'                         MPVERSION = 5) 
#' }
#' in this way you can use the same data object for repeated function calls, which speeds functions up.
#' It also opens up functions for use with data that did not necessarily originate in the HMD input files.
#' As long as, in this case, \code{cohComp} has the correct columns, classes, dimensions and data, it 
#' should work just fine (except for version 6 exposures, which draw on monthly births). Those will be 
#' moved in-memory at a later date.
#' 
#' @references Wilmoth, J.R. and Andreev, K. and Jdanov, D. and Glei, DA (2007). 
#' Methods protocol for the human mortality database. University of California, 
#' Berkeley, and Max Planck Institute for Demographic Research, Rostock, 
#' [version 31/05/2007] \url{http://mortality.org}.
#' @docType package
#' @name RLifeTables
NULL

###############################################################################
# country.lookup dataset Documentation:
#' A lookup table for HMD country abbreviations and full names.
#' 
#' A \code{data.frame} containing HMD country abbreviations and their full names. 
#'   The variables, both character class, are as follows:
#' 
#' \itemize{
#'   \item cntries.
#'   \item full.countries.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A \code{data.frame} with 47 rows and 2 variables
#' @name country.lookup
NULL

###############################################################################
# monthDurations dataset Documentation:
#' A lookup table for the number of days in each month
#' 
#' An integer \code{matrix} containing the number of days in each month from 1700 until 2100.
#' Months (1-12) are in rows and years (1700-2100) are in columns. This \code{matrix} is used in
#' MPv6 exposure calculations only.
#' 
#' @docType data
#' @keywords datasets
#' @format A \code{matrix} with 12 rows and 401 columns
#' @name monthDurations
NULL
