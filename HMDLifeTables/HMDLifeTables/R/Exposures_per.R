#' @title \code{Exposures_per} a function to prepare either v5 or v6 period exposures.
#'
#' @description This function is meant to be internal. It is called by
#'   \code{Exposures_Deaths_Mx_AxN()}, \code{ltper_AxN()} and
#'   \code{ltperBoth_AxN()}. If necessary, it calls
#'   \code{getPeriodComponents()}. Exposures from the appropriate MP version are
#'   returned in an 'age by year' matrix already cut down to \code{OPENAGE}.
#'
#' @details This function can take data arguments in-memory, from the \code{/Rbin/} folder, or else derive them by calling \code{getPeriodComponents()}. This function is separate because using monthly birth data for v6 exposures is a fair amount of code. Since a few functions need these exposures it's best to modularize and call where needed. This function uses a package dataset \code{monthDurations.rda} in order to save having to recompute it on each run- this is a matrix with the number of days in each month every year across several centuries, reaching way beyond the potential of the HMD to as to avoid future bugs. It gets cut to size as needed. This function also calls \code{AC2AP()}.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param pop1 optional. An 'age by year' matrix of Jan 1st population counts. Taken from \code{perComp} or derived if necessary.
#' #' @param pop2 optional. An 'age by year' matrix of Dec 31st population counts territorially adjusted back from the january 1st year t+1 back to year t. Taken from \code{perComp} or derived if necessary.
#' @param dl optional. An 'age by year' matrix of lower triangle death counts. Taken from \code{perComp} or derived if necessary.
#' @param du optional. An 'age by year' matrix of upper triangle death counts. Taken from \code{perComp} or derived if necessary.
#' @param perComp optional. Output from \code{getPeriodComponents()}. 
#' @param sex \code{"m"} or \code{"f"}.
#' @param OPENAGE the desired open age. Default value is 110
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/cohExposuresRaw.Rdata} as well.
#' @param MPVERSION 5 or 6. Version 5 exposures assume uniformity, v6 allows for non-uniformity across birthdays in the death distribution by taking information from monthly birth distributions. 
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#'  
#' @return a list of 'age by year' matrices of exposures: \code{Exp, Eu, El}, corresponding to total exposures, upper triangle exposures, and lower triangle exposures.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe

## TODO: CAB - restructure to handle empty XYZbirthbymonth.txt monthly birth files and 
##             require its presence under V6.  Change logic so that immutable parameter
##             MPVERSION does not need to be reassigned

###############################################################################
Exposures_per <- function(WORKING = getwd(), 
  pop1 = NULL,    
  pop2 = NULL,
  dl = NULL,
  du = NULL, 
  perComp = NULL,
  sex ,   # CAB: should not have a default
  OPENAGE = 110, 
  save.bin = TRUE, 
  MPVERSION,  # version must be explicit
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL
  ){
  # MPVERSION can only be 5 or 6
  if(!MPVERSION %in% c(5, 6)){
    stop("\nonly MP versions 5 and 6 are supported at this time\n")
  }
  # --------------------
  # get country short name
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH <- file.path(WORKING, "InputDB")
  }
  
  # --------------------
  # this chunk bypassed if pop1, pop2, dl, du supplied as arguments, otherwise derived, loaded or produced...
  if (is.null(pop1) | is.null(pop2) | is.null(dl) | is.null(du)){
    # typical period objects required   # and also allow load() of Rbin if produced previously
    if (is.null(perComp)){
      percomp.path <- file.path(WORKING, "Rbin", paste0(sex, "_periodComponents.Rdata"))
      if (file.exists(percomp.path)){
        perComp <- local(get(load(percomp.path)))
      } else {
        perComp <- getPeriodComponents(
          WORKING = WORKING, 
          sex = sex,
          OPENAGE = OPENAGE,
          save.bin = TRUE,
          XXX = XXX,
          LDBPATH = LDBPATH,
          IDBPATH = IDBPATH)
      }
    }
    # shape to AP matrices
    pop1          <- acast(perComp, Age ~ Year, sum, value.var = "pop1") 
    pop2          <- acast(perComp, Age ~ Year, sum, value.var = "pop2") 
    dl            <- acast(perComp, Age ~ Year, sum, value.var = "dl")
    du            <- acast(perComp, Age ~ Year, sum, value.var = "du")
  }

  # ---------------------------------------------------------------------------
  # TR: I intuit that this has been modified by CB. And that the monthly file name has changed??
  use.old.exposure.formula <- TRUE
  births.monthly.path <- file.path(IDBPATH, paste0(XXX, "birthbymonth.txt"))
  if (MPVERSION > 5){
    use.old.exposure.formula <- FALSE 
    
    if (!file.exists(births.monthly.path)){
      cat("\nMPVERSION was given as", MPVERSION, "but monthly births file was missing:\n", 
            births.monthly.path, "\nreverted to MPVERSION 5 exposures\n")
        
      use.old.exposure.formula <- TRUE
        
      } 
  }
# old exposures, considerably simpler :-)
  if (use.old.exposure.formula){
    Eu          <- (pop1 / 2) -  (du / 6)
    El          <- (pop2 / 2) +  (dl / 6)
    #Exp         <- El + Eu
    Exp         <- (pop1 + pop2) / 2 + (dl - du) / 6
    
    # optional save out
    if (save.bin){
      out.path0   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresRaw.Rdata"))
      save(Exp, file = out.path0)
      out.pathL   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresLRaw.Rdata"))
      save(El, file = out.pathL)
      out.pathU   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresURaw.Rdata"))
      save(Eu, file = out.pathU)
      # now save mx
      #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", out.path0))
    }
    return(list(Exp=Exp, Exp.u=Eu, Exp.l=El) )
  }
  
# -------------------------------------
  nameBasedRecode<- function(.x,.y){  #substitute y into x based on names attribute of x, y
    stopifnot( !is.null(attributes(.x)[["names"]]) & !is.null(attributes(.y)[["names"]]))
    .x[ names(.y)[names(.y) %in% names(.x)] ] <- .y[ names(.y) %in% names(.x)]
    return(.x)
  }
  # -------------------------------------  
# set up dims, names used throughout
  yrs     <- as.integer(colnames(pop1))
  cohs    <- (min(yrs) - nrow(pop1) - 1):(max(yrs) + 1)
  Ncohs   <- length(cohs)
  ages    <- 0:OPENAGE
  Nages   <- length(ages)
# TODO: check for monthly in Rbin, as LDB can save this,
# be sure to use optional LDBPATH
  
  BM              <- read.table(births.monthly.path,
                      header = TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE, 
                      na.strings = ".")

# remove TOT, convert Month to integer for cast()ing
  BM              <- BM[BM$Month != "TOT" & BM$Month != "UNK", ]
  
#  unks            <- FALSE
#  if (any(BM$Month == "UNK")) {
#    unks            <- TRUE
#    UNKs            <- BM[BM$Month == "UNK", ]
#    BM              <- BM[!BM$Month == "UNK", ]
#    UNKs            <- acast(UNKs, Month ~ Year, sum, value.var = "Births")
#    #colSums(UNKs)
#  } 
  
  BM$Month        <- as.integer(BM$Month)
  
# create Month x Year (cohort) matrix
  Bmat            <- acast(BM[BM$LDB == 1, ], Month ~ Year, sum, value.var = "Births")

# distribute births of unknown month proportionally
#  if (unks){
#    Bmat.w.UNKs   <- Bmat[, colnames(UNKs), drop = FALSE] # need to keep 2 dims for JPN..
#    Bmat.pdf      <- t(t(Bmat.w.UNKs) / colSums(Bmat.w.UNKs))
#    Bmat[, colnames(UNKs)] <- Bmat.w.UNKs + t(t(Bmat.pdf) * c(UNKs))
#  }
  
# --------------------------------------------------------------
# Edit: load() durations matrix, since it's static
  #durations       <- monthDurations # TODO: check if throw error, mabe data(monthDurations)
  monthDurations       <- monthDurations[, colnames(Bmat)] # cut down to Bmat years

# --------------------------------------------------------------
  # John's formula, per nov 12, 2012 email:        
  f.i                       <- t(t(Bmat) / colSums(Bmat))
  BM                        <- apply(monthDurations, 2, cumsum)
  b.i                       <- rbind(0, t(t(BM) / colSums(monthDurations)))
  b.bar                     <- colSums(f.i * (b.i[1:12, ] + b.i[2:13, ]) / 2)
  b.bar.full                <- rep(.5, Ncohs)
  names(b.bar.full)         <- cohs
 
  ## name-based substitution
  b.bar.full <- nameBasedRecode(b.bar.full, b.bar)
  #b.bar.full[ names(b.bar)[names(b.bar) %in% names(b.bar.full)] ] <- b.bar[names(b.bar) %in% names(b.bar.full) ]
  
  #an AC B.bar matrix
  b.bar.mat                 <- matrix(b.bar.full, 
    nrow = Nages, 
    ncol = Ncohs, 
    byrow = TRUE, 
    dimnames = list(ages, cohs))
  
  sigmasq                     <- colSums(f.i * ((b.i[1:12, ] ^ 2 + b.i[2:13, ] * b.i[1:12, ] + b.i[2:13, ] ^ 2) / 3)) - b.bar ^ 2
  sigmasq.full                <- rep(1 / 12, Ncohs) # default uniform distribution
  names(sigmasq.full)         <- cohs

  #sigmasq.full[names(sigmasq)]  <- sigmasq
  sigmasq.full <- nameBasedRecode(sigmasq.full, sigmasq)
  sigmasq.mat                 <- matrix(sigmasq.full, 
                                nrow = Nages, 
                                ncol = Ncohs, 
                                byrow = TRUE, 
                                dimnames = list(ages, cohs))
  
# --------------------------------------------------------------------
# next part stays same after sigma formula determined.
# --------------------------------------------------------------------
# convert to AP, large dimensions, will need to extract relevant years
# AP1 and AP2 are numerically equal, but have different Year indexing, used
# for correct alignment of triangles
  b.bar.AP1                  <- AC2AP(b.bar.mat, Lexis = 2)
  sigmasq.AP1                <- AC2AP(sigmasq.mat, 2)
  b.bar.AP2                  <- AC2AP(b.bar.mat, Lexis = 1)
  sigmasq.AP2                <- AC2AP(sigmasq.mat, 1)
# sigma1 is the time at birth for cohort t-x-1. b2.bar for cohort t-x. These are PERIOD matrices, yay
# NB: 1 is for upper triangle here, contra everywhere else in HMD, so be careful!

  sigmasq1                  <- sigmasq.AP1[1:nrow(pop1), colnames(pop1)]
  sigmasq2                  <- sigmasq.AP2[1:nrow(pop2), colnames(pop2)]
# b.bar is the mean time at birth (proportion of year completed) 
  b1.bar                    <- b.bar.AP1[1:nrow(pop1), colnames(pop1)] # for cohort t-x-1
  b2.bar                    <- b.bar.AP2[1:nrow(pop2), colnames(pop2)] # for cohort t-x
# -----------------------------------------------------------
# now the actual formulas, per the memo and Dima's correction
# l1,l2 changed to s1, s2
  s1            <- 1 - b2.bar                               # insert formulas
  s2            <- (1 - b2.bar) / 2 - sigmasq2 / (2 * (1 - b2.bar))
  u1            <- b1.bar
  u2            <- b1.bar / 2 - sigmasq1 / (2 * b1.bar)
  
  El            <- s1 * pop2 + s2 * dl
  Eu            <- u1 * pop1 - u2 * du

# version 6 period exposures:
  Exp           <- El + Eu
  
  # optional save out
  if (save.bin){
    out.path0   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresRaw.Rdata"))
    save(Exp, file = out.path0)
    out.pathL   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresLRaw.Rdata"))
    save(El, file = out.pathL)
    out.pathU   <- file.path(WORKING, "Rbin",  paste0(sex, "perExposuresURaw.Rdata"))
    save(Eu, file = out.pathU)
    # now save mx
    #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }

  invisible(return(list(Exp=Exp, Exp.u=Eu, Exp.l=El)) )
}


