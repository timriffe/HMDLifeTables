#' @title \code{cExposures_and_cMx_AxN} a function to prepare unrounded cExposures and cMx data.frames
#'
#' @description IFF cohort measures are published for country \code{XXX}, this function will prepare the \code{data.frame}s for cMx and cExposures, saving them to \code{/Rbin/} in order to be used as such or else to be formatted for publishing later by WriteExposuresDeathsMx(). 
#'
#' @details This function has no native check for whether or not results should be run. That check needs to happen before this function runs. In case there aren't a sufficient number of cohorts with \code{min.years} worth of data, no warning is given at this time. If you run this function on a country with less than \code{min.years} worth of period data, it will throw an uninformative native R error, as no check is done.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param males the cohort components \code{data.frame} object, as created by \code{getCohortComponents()} for males. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param females the cohort components \code{data.frame} object, as created by getCohortComponents() for females. If left as \code{NULL}, \code{getCohortComponents()} is called.
#' @param abridged logical. Default \code{FALSE}.
#' @param N integer number of years to aggregate 1, 5 or 10. Default 1 (no aggregatation)
#' @param OPENAGE the desired open age. Default value is 110
#' @param min.years the minimum number of data years that we need in order to publish anything for a cohort. Default 30.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects the selection of which cohorts to show. For v5, matlab programs had an out-of-protocol but reasonable requirement of \code{min.years} continuous years of at most age 100. v6 will take the \code{min.years} consecutive years literally unless changed in the MP.
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/cExposures_1x1.Rdata} as well. Appropriate name is derived systematically. In this case both objects are saved separately.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the past path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a list of two \code{data.frame}s of the full cMx and cExposures output, in long format. These have 5 columns: \code{Year}, \code{Age}, \code{Female}, \code{Male} and \code{Total}. Age is a formatted character vector, while the other columns are numeric.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
cExposures_and_cMx_AxN <- function(
  WORKING = getwd(), 
  males = NULL,
  females = NULL,
  abridged = FALSE, 
  N = 1,
  OPENAGE = 110,
  min.years = 30,
  MPVERSION , # explicit, no default
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL
  ){
    if (is.null(XXX)){
      XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
    }
    if (is.null(LDBPATH)){
      LDBPATH <- file.path(WORKING, "LexisDB")
    }
    if (is.null(IDBPATH)){
      IDBPATH <- file.path(WORKING, "InputDB")
    }
# -----------------------------------------------------
    # read in data, or generate if not available:
  if (is.null(males) | is.null(females)){
    cohcomp.path.f <- file.path(WORKING, "Rbin", "f_cohortComponents.Rdata")
    cohcomp.path.m <- file.path(WORKING, "Rbin", "m_cohortComponents.Rdata")
    if (file.exists(cohcomp.path.f) & file.exists(cohcomp.path.m) ){
      females <- local(get(load(cohcomp.path.f)))
      males   <- local(get(load(cohcomp.path.m)))
    } else {
      # generate data
      females <- getCohortComponents(
        WORKING = WORKING, 
        sex = "f",
        save.bin = TRUE,
        XXX = XXX,
        LDBPATH = LDBPATH,
        IDBPATH = IDBPATH)
      males <- getCohortComponents(
        WORKING = WORKING, 
        sex = "m",
        save.bin = TRUE,
        XXX = XXX,
        LDBPATH = LDBPATH,
        IDBPATH = IDBPATH)
    }
  }
# ---------------------------------------------------------------------------------
# get matrices to work with
    pop.f       <- acast(females, Age ~ Year, value.var = "pop")
    dl.f        <- acast(females, Age ~ Year, value.var = "dl")
    du.f        <- acast(females, Age ~ Year, value.var = "du")
    pop.m       <- acast(males, Age ~ Year, value.var = "pop")
    dl.m        <- acast(males, Age ~ Year, value.var = "dl")
    du.m        <- acast(males, Age ~ Year, value.var = "du")
# ---------------------------------------------------------------------------------
    # produce exposures externally
    Exp.f       <- Exposures_coh(WORKING = WORKING, 
                      pop = pop.f,    
                      dl = dl.f,
                      du = du.f, 
                      sex = "f", 
                      OPENAGE = OPENAGE, 
                      save.bin = TRUE,   # CAB test
                      MPVERSION = MPVERSION,
                      XXX = NULL,
                      LDBPATH = LDBPATH,
                      IDBPATH = IDBPATH)
    Exp.m       <- Exposures_coh(WORKING = WORKING, 
                      pop = pop.m,    
                      dl = dl.m,
                      du = du.m, 
                      sex = "m", 
                      OPENAGE = OPENAGE, 
                      save.bin = TRUE, 
                      MPVERSION = MPVERSION,
                      XXX = NULL,
                      LDBPATH = LDBPATH,
                      IDBPATH = IDBPATH)
# ---------------------------------------------------------------------------------
# imputeNAsfrom du          # this incidentally is the non-intuitive step

    dl.f[is.na(du.f)]   <- NA     # that was necessary to replicate original
    pop.f[is.na(du.f)]  <- NA     # matlab results, and that is not spelled out 
    dl.m[is.na(du.m)]   <- NA     
    pop.m[is.na(du.m)]  <- NA     
    
    # makes sense after looking at MP diagram for a while; 
    # ?? cut back partial cohort parallelogram 
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# CAB revisit this logic and see if can be found in MP
# Roll everything back to open age at this step
    i.OPENAGE <- OPENAGE + 1
# sum to age 110+, in case there are later events   
    NA110                   <- is.na(dl.f[i.OPENAGE, ])
    dl.f[i.OPENAGE, ]       <- colSums(dl.f[i.OPENAGE:nrow(dl.f),], na.rm = TRUE)
    du.f[i.OPENAGE, ]       <- colSums(du.f[i.OPENAGE:nrow(du.f),], na.rm = TRUE)
    pop.f[i.OPENAGE, ]      <- colSums(pop.f[i.OPENAGE:nrow(pop.f),], na.rm = TRUE)
    dl.m[i.OPENAGE, ]       <- colSums(dl.m[i.OPENAGE:nrow(dl.m),], na.rm = TRUE)
    du.m[i.OPENAGE, ]       <- colSums(du.m[i.OPENAGE:nrow(du.m),], na.rm = TRUE)
    pop.m[i.OPENAGE, ]      <- colSums(pop.m[i.OPENAGE:nrow(pop.m),], na.rm = TRUE)  
# cut down
    dl.f                    <- dl.f[1:i.OPENAGE, ]
    du.f                    <- du.f[1:i.OPENAGE, ]
    pop.f                   <- pop.f[1:i.OPENAGE, ]
    dl.m                    <- dl.m[1:i.OPENAGE, ]
    du.m                    <- du.m[1:i.OPENAGE, ]
    pop.m                   <- pop.m[1:i.OPENAGE, ]
    
# keep age 110 NAs in tact   
## CAB: why is age 110 forced to NA iff dl.f was NA ?    
    dl.f[i.OPENAGE, NA110]  <- NA
    du.f[i.OPENAGE, NA110]  <- NA
    pop.f[i.OPENAGE, NA110] <- NA
    dl.m[i.OPENAGE, NA110]  <- NA
    du.m[i.OPENAGE, NA110]  <- NA
    pop.m[i.OPENAGE, NA110] <- NA
    
# ---------------------------------------------------------------------------------   
# keep only cohorts with 30 years of non-zero data (exclude NAs and 0s)    
# TODO: look for this in MP revison 6
    if (MPVERSION == 5){
      # this just uses the same criteria as the matlab: the 100+ thing is out-of-protocol
      # generally these give the same results, e.g. for Iceland. It's also reasonable, and 
      # perhaps worth keeping, since such ages would only be from early cohorts where old
      # ages are often composed of broken open age groups.
      # CAB: various changes here do not match up exactly with matlab, which has a slightly different 
      # criterion for choosing candidate cohorts (not just based on pop.f)
      c.keep                  <- colSums(!(is.na(pop.f[1:109, ])), na.rm = TRUE) >= min.years
      # c.keep                  <- colSums(!(is.na(pop.f[1:101, ])), na.rm = TRUE) >= min.years
    } else {
      c.keep                  <- colSums(!(is.na(pop.f) | (pop.f + pop.m + dl.f + dl.m + du.f + du.m) == 0), na.rm = TRUE) >= min.years
    }
    
    # offer stop break if there aren't any cohorts with 30 years of observations
    if (sum(c.keep) == 0){
      #cat("Not enough data to publish cohort exposures and deaths\nneed ", min.years, "worth of data")
      return(NULL)
    }
    
    pop.f     <- pop.f[, c.keep]
    dl.f      <- dl.f[, c.keep]
    du.f      <- du.f[, c.keep]
    pop.m     <- pop.m[, c.keep]
    dl.m      <- dl.m[, c.keep]
    du.m      <- du.m[, c.keep]
    Exp.f     <- Exp.f[, colnames(pop.f)]
    Exp.m     <- Exp.m[, colnames(pop.m)]
# ---------------------------------------------------------------------------------
# Aggregate years if necessary:
    dl.f      <- YearAgg(dl.f, N = N)
    du.f      <- YearAgg(du.f, N = N)
    pop.f     <- YearAgg(pop.f, N = N)
    dl.m      <- YearAgg(dl.m, N = N)
    du.m      <- YearAgg(du.m, N = N)
    pop.m     <- YearAgg(pop.m, N = N)
    Exp.f     <- YearAgg(Exp.f, N = N)
    Exp.m     <- YearAgg(Exp.m, N = N)
# ---------------------------------------------------------------------------------
# calculate exposures
    #Exp.m     <- pop.m + (1 / 3) * (dl.m - du.m)
    #Exp.f     <- pop.f + (1 / 3) * (dl.f - du.f)
    
    Dx.f      <- du.f + dl.f
    Dx.m      <- du.m + dl.m
    
    ages      <- c(0:(OPENAGE - 1), paste0(OPENAGE, "+"))
# ---------------------------------------------------------------------------------
# abridge if necessary
    if (abridged){
      # fac is a grouping factor for tapply()
      fac   <- c(0, 1, 1, 1, 1, 5:OPENAGE - 5:OPENAGE %% 5)
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
    # get total exposure
    Exp.t <- Exp.m + Exp.f
# ---------------------------------------------------------------------------------    
    # get Mx ready
    Mx.f  <- Dx.f / Exp.f
    Mx.m  <- Dx.m / Exp.m
    Mx.t  <- (Dx.f + Dx.m) / Exp.t
    
    years <- colnames(Exp.t)
# ---------------------------------------------------------------------------------    
    # create output object
    cExposures <- data.frame(
      Year = rep(years, each = length(ages)),
      Age = rep(ages, length(unique(years))),
      Female = as.vector(Exp.f),
      Male = as.vector(Exp.m),
      Total = as.vector(Exp.t),
      stringsAsFactors = FALSE
    )
    cMx <- data.frame(
      Year = rep(years, each = length(ages)),
      Age = rep(ages, length(unique(years))),
      Female = as.vector(Mx.f),
      Male = as.vector(Mx.m),
      Total = as.vector(Mx.t),
      stringsAsFactors = FALSE
    )
    
    # optionally save to R binary format 
    if (save.bin){
      # long file name to mark parameters
      out.name1 <- paste0("cExposures",
        ifelse(abridged, "_5", "_1"),"x",
        ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
      out.name2 <- paste0("cMx",
        ifelse(abridged, "_5", "_1"),"x",
        ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
      dir.path <- file.path(WORKING, "Rbin")
      if(!file.exists(dir.path)) {
        dir.create(dir.path)
        #system(paste0("chgrp hmdcalc ", dir.path))
        #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      }
      out.path1 <- file.path(dir.path, out.name1)
      out.path2 <- file.path(dir.path, out.name2)
      # saving happens here
      save(cExposures, file = out.path1)
      save(cMx, file = out.path2)

      #system(paste0("chgrp hmdcalc ", out.path1))
      #system(paste0("chgrp hmdcalc ",out.path2))
      #Sys.chmod(c(out.path1, out.path2), mode = "2775", use_umask = FALSE)

    }
    # invisibly return output, so that the console isn't cluttered up
    output <- list(cExposures = cExposures, cMx = cMx)
    invisible(output)
  }
