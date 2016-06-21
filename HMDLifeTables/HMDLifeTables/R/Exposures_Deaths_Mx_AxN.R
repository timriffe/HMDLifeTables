#' @title \code{Exposures_Deaths_Mx_AxN} a function to prepare unrounded Exposures, Deaths and Mx data.frames
#'
#' @description All countries produce Deaths_1x1, Deaths_1x5.. 1x10, 5x1, 5x5 and 5x10, and likewise for Exposures and Mx estimates. Since Mx requires the former two in order to be calculated, and all items have a similar output construction, they are all produced in one function. Output is returned as a list of (long) \code{data.frame}s, and optionally saved as appropriately named .Rdata objects in the /Rbin folder.
#'
#' @details This function is a top-level function, and can be called independently. If arguments \code{males}, \code{females} are missing, they are created by calls to \code{getPeriodComponents()}
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation.  Default \code{getwd()}.
#' @param males the cohort components \code{data.frame} object, as created by \code{getCohortComponents()} for males. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param females the cohort components \code{data.frame} object, as created by getCohortComponents() for females. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param abridged logical. Default \code{FALSE}.
#' @param N integer number of years to aggregate 1, 5 or 10. Default 1 (no aggregatation).
#' @param OPENAGE the desired open age. Default value is 110.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects exposure calculation. For both versions, exposure calculations are performed by \code{Exposures_per()}.
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/cExposures_1x1.Rdata} as well. Appropriate name is derived systematically. In this case both objects are saved separately.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a list of three \code{data.frame}s of the full  Deaths, Mx and Exposures output, in long format. These have 5 columns: \code{Year}, \code{Age}, \code{Female}, \code{Male} and \code{Total}. Age is a formatted character vector, while the other columns are numeric.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
Exposures_Deaths_Mx_AxN <- function(
  WORKING = getwd(), 
  males = NULL,
  females = NULL,
  abridged = FALSE, 
  N = 1,
  OPENAGE = 110,
  MPVERSION = 5,
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL
){

# ----------------------------------------------------
# reshape data into matrices
  # read in data
  # 'XXX' (country abbrev) is used here and there
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH <- file.path(WORKING, "InputDB")
  }
  if (is.null(males) | is.null(females)){
    percomp.path.f <- file.path(WORKING, "Rbin", "f_periodComponents.Rdata")
    percomp.path.m <- file.path(WORKING, "Rbin", "m_periodComponents.Rdata")
    if (file.exists(percomp.path.f) & file.exists(percomp.path.m)){
      females <- local(get(load(percomp.path.f)))
      males   <- local(get(load(percomp.path.m)))
    } else {
      females <- getPeriodComponents(
        WORKING = WORKING, 
        sex = "f",
        OPENAGE = OPENAGE,
        save.bin = TRUE,
        XXX = XXX,
        LDBPATH = LDBPATH,
        IDBPATH = IDBPATH)
      males <- getPeriodComponents(
        WORKING = WORKING, 
        sex = "m",
        OPENAGE = OPENAGE,
        # MPVERSION can only b
        save.bin = TRUE,
        XXX = XXX,
        LDBPATH = LDBPATH,
        IDBPATH = IDBPATH)
    }
  }
# now reshape to 4 matrices:
  pop1.f      <- acast(females, Age ~ Year, value.var = "pop1")
  pop2.f      <- acast(females, Age ~ Year, value.var = "pop2")
  dl.f        <- acast(females, Age ~ Year, value.var = "dl")
  du.f        <- acast(females, Age ~ Year, value.var = "du")
  pop1.m      <- acast(males, Age ~ Year, value.var = "pop1")
  pop2.m      <- acast(males, Age ~ Year, value.var = "pop2")
  dl.m        <- acast(males, Age ~ Year, value.var = "dl")
  du.m        <- acast(males, Age ~ Year, value.var = "du")
  
  # calculate exposure
  # externalize exposure calculation
  Exp.f       <- Exposures_per(WORKING = WORKING, 
                                pop1 = pop1.f,    
                                pop2 = pop2.f,
                                dl = dl.f,
                                du = du.f, 
                                sex = "f", 
                                OPENAGE = OPENAGE, 
                                save.bin = TRUE, 
                                MPVERSION = MPVERSION
                              )
  Exp.m       <- Exposures_per(WORKING = WORKING, 
                                pop1 = pop1.m,    
                                pop2 = pop2.m,
                                dl = dl.m,
                                du = du.m, 
                                sex = "m", 
                                OPENAGE = OPENAGE, 
                                save.bin = TRUE, 
                                MPVERSION = MPVERSION
                              )
  # ---------------------------------------------------------------------------------
  # Aggregate years if necessary: (if N = 1, does nothing)
  Exp.f        <- YearAgg(Exp.f, N = N)
  Exp.m        <- YearAgg(Exp.m, N = N)
  Dx.f         <- YearAgg(dl.f + du.f, N = N)
  Dx.m         <- YearAgg(dl.m + du.m, N = N)

  ages <- c(0:(OPENAGE - 1), paste0(OPENAGE, "+"))
  # abridge, if necessary
  if (abridged){
    # fac is a grouping factor for tapply()
    fac <- c(0, 1, 1, 1, 1, 5:OPENAGE - 5:OPENAGE %% 5)
    Exp.f <- apply(Exp.f, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    Exp.m <- apply(Exp.m, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    Dx.f <- apply(Dx.f, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    Dx.m <- apply(Dx.m, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    # redo ages if necessary: these are width-formatted
    ages <- c("    0  ", 
      paste(sprintf("%3s", c(1, seq(5, OPENAGE - 5, by = 5))),
        sprintf("%-3s", c(4, seq(9, OPENAGE - 1, by = 5))), 
        sep = "-"),
      paste0(OPENAGE, "+   ")
    )        
  }
# get Dx and Exp for total population
  Exp.t <- Exp.f + Exp.m
  Dx.t  <- Dx.f + Dx.m
# get Mx
  Mx.m  <- Dx.m / Exp.m
  Mx.f  <- Dx.f / Exp.f
  Mx.t  <- Dx.t / Exp.t
  
  years <- colnames(Exp.t)
  
  # create output object
  Exposures <- data.frame(
    Year = rep(years, each = length(ages)),
    Age = rep(ages, length(unique(years))),
    Female = as.vector(Exp.f),
    Male = as.vector(Exp.m),
    Total = as.vector(Exp.t),
    stringsAsFactors = FALSE)
  Deaths <- data.frame(
    Year = rep(years, each = length(ages)),
    Age = rep(ages, length(unique(years))),
    Female = as.vector(Dx.f),
    Male = as.vector(Dx.m),
    Total = as.vector(Dx.t),
    stringsAsFactors = FALSE)
  Mx <- data.frame(
    Year = rep(years, each = length(ages)),
    Age = rep(ages, length(unique(years))),
    Female = as.vector(Mx.f),
    Male = as.vector(Mx.m),
    Total = as.vector(Mx.t),
    stringsAsFactors = FALSE)
  # optionally save to R binary format (for instance for post production of both-sex tables)
  if (save.bin){
    # long file name to mark parameters
    out.name1 <- paste0("Exposures",
      ifelse(abridged, "_5", "_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    out.name2 <- paste0("Deaths",
      ifelse(abridged, "_5", "_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    out.name3 <- paste0("Mx",
      ifelse(abridged, "_5", "_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
    }
    out.path1 <- file.path(dir.path, out.name1)
    out.path2 <- file.path(dir.path, out.name2)
    out.path3 <- file.path(dir.path, out.name3)
    # saving happens here
    save(Exposures, file = out.path1)
    save(Deaths, file = out.path2)
    save(Mx, file = out.path3)
    #Sys.chmod(c(out.path1, out.path2, out.path3), mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path1))
    #system(paste0("chgrp hmdcalc ", out.path2))
    #system(paste0("chgrp hmdcalc ", out.path3))
  }
  # invisibly return output, so that the console isn't cluttered up
  output <- list(Exposures = Exposures, Deaths = Deaths, Mx = Mx)
  invisible(output)
}




