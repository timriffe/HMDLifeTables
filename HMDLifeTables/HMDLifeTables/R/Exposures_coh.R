#' @title \code{Exposures_coh} a function to prepare either v5 or v6 cohort exposures.
#'
#' @description This function is meant to be internal. It is called by \code{cExposures_and_cMx_AxN()}, \code{ltcoh_AxN()} and \code{ltcohBoth_AxN()}. If necessary, it calls \code{getCohortComponents()}. Exposures from the appropriate MP version are returned in an 'age by year' matrix already cut down to \code{OPENAGE}, but with all cohorts present that are present in the uncut \code{pop}- i.e. cohort selection for lifetables or exposure writing needs to happen downstream.
#'
#' @details This function can take data arguments in-memory, from the \code{/Rbin/} folder, or else derive them by calling \code{getCohortComponents()}. This function is separate because using monthly birth data for v6 exposures is a fair amount of code. Since a few functions need these exposures it's best to modularize and call where needed. This function uses a package dataset \code{monthDurations.rda} in order to save having to recompute it on each run- this is a matrix with the number of days in each month every year across several centuries, reaching way beyond the potential of the HMD to as to avoid future bugs. It gets cut to size as needed.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param pop optional. An 'age by year' matrix of Jan 1st population counts. Taken from \code{cohComp} or derived if necessary.
#' @param dl optional. An 'age by year' matrix of year t lower triangle death counts. Taken from \code{cohComp} or derived if necessary.
#' @param du optional. An 'age by year' matrix of year t+1 upper triangle death counts, territorially adjusted back to the year t standard. Taken from \code{cohComp} or derived if necessary.
#' @param cohComp optional. Output from \code{getCohortComponents()}. 
#' @param sex \code{"m"} or \code{"f"}.
#' @param OPENAGE the desired open age. Default value is 110
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/cohExposuresRaw.Rdata} as well.
#' @param MPVERSION 5 or 6. Version 5 exposures assume uniformity, v6 allows for non-uniformity across birthdays in the death distribution by taking information from monthly birth distributions. Differences are negligible.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return an 'age by year' matrix of exposures. if in test mode, a list with matrices of both deaths and exposures is returned.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe

## TODO: CAB - restructure to handle empty XYZmonthly.txt monthly birth files and 
##             require its presence under V6.  Change logic so that immutable parameter
##             MPVERSION does not need to be reassigned

###############################################################################
Exposures_coh <- function(WORKING = getwd(),
  pop = NULL,    
  dl = NULL,
  du = NULL, 
  cohComp = NULL,
  sex,   #CB, no defaults or ambiguous arg
  OPENAGE = 110, 
  save.bin = TRUE, 
  MPVERSION , # explicit, no default
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
  if (is.null(pop) | is.null(dl) | is.null(du)){
    if (is.null(cohComp)){
      cohcomp.path <- file.path(WORKING, "Rbin", paste0(sex, "_cohortComponents.Rdata"))
      if (file.exists(cohcomp.path)){ # tip thanks to Charles on SO
        cohComp <- local(get(load(cohcomp.path)))
        
      } else {
        # load cohort components function
         # generate data
        cohComp <- getCohortComponents(
          WORKING = WORKING, 
          sex = sex,
          save.bin = TRUE,
          XXX = XXX,
          LDBPATH = LDBPATH,
          IDBPATH = IDBPATH)
      }
    }
    # require(reshape2)
  # now reshape to 3 matrices, age x year
  # dl and du are pre-staggered, though they have the same column names
  # pop refers to the parallelogram central population.
    pop       <- acast(cohComp, Age ~ Year, value.var = "pop")
    dl        <- acast(cohComp, Age ~ Year, value.var = "dl")
    du        <- acast(cohComp, Age ~ Year, value.var = "du")
  }
  # -------------------------------------------
  i.openage         <- OPENAGE + 1
  # aggregate open age interval, cut down to size:
  # pop
  openNas           <- is.na(pop[i.openage, ])
  pop[i.openage, ]  <- colSums(pop[i.openage:nrow(pop), ], na.rm = TRUE)
  pop               <- pop[1:i.openage, ]
  pop[i.openage, openNas] <- NA
  # dl
  openNas           <- is.na(dl[i.openage, ])
  dl[i.openage, ]   <- colSums(dl[i.openage:nrow(dl), ], na.rm = TRUE)
  dl                <- dl[1:i.openage, ]
  dl[i.openage, openNas] <- NA
  # du
  openNas           <- is.na(du[i.openage, ])
  du[i.openage, ]   <- colSums(du[i.openage:nrow(du), ], na.rm = TRUE)
  du                <- du[1:i.openage, ]
  du[i.openage, openNas] <- NA
  # -------------------------------------------

  use.old.exposure.formula <- TRUE
  births.monthly.path <- file.path(IDBPATH, paste0(XXX, "birthbymonth.txt"))
  if (MPVERSION > 5){
    use.old.exposure.formula <- FALSE 
    
    if (!file.exists(births.monthly.path)){
      cat("\nMPVERSION was given as", MPVERSION, "but monthly births file was missing:\n", births.monthly.path, "\nreverted to MPVERSION 5 exposures\n")
      use.old.exposure.formula <- TRUE
    } 
  }
  
  # old exposures, considerably simpler :-)
  if (use.old.exposure.formula){
    Exp       <- pop + (1 / 3) * (dl - du)
    # optional save out
    if (save.bin){
      out.path0   <- file.path(WORKING, "Rbin",  paste0(sex, "cohExposuresRaw.Rdata"))
      save(Exp, file = out.path0)
      # now save mx
      #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", out.path0))
    }
    # YES cut off after open age
    return(Exp[1:i.openage, ]) 
  }
  
  # -------------------------------------
  nameBasedRecode<- function(.x,.y){  #substitute y into x based on names attribute of x, y
    stopifnot( !is.null(attributes(.x)[["names"]]) & !is.null(attributes(.y)[["names"]]))
    .x[ names(.y)[names(.y) %in% names(.x)] ] <- .y[ names(.y) %in% names(.x)]
    return(.x)
  }
  # -------------------------------------  
  
#  cohs              <- as.integer(colnames(pop))
#  Ncohs             <- length(cohs)
  ages              <- 0:OPENAGE
  Nages             <- length(ages)
  
  BM                <- read.table(births.monthly.path,
                        header = TRUE, 
                        sep = ",", 
                        stringsAsFactors = FALSE, 
                        na.strings = ".")
  
# remove TOT, convert Month to integer for cast()ing
  BM                <- BM[BM$Month != "TOT" & BM$Month != "UNK" & BM$LDB == 1, ]
  
#  unks              <- FALSE
#  if (any(BM$Month == "UNK")) {
#    unks            <- TRUE
#    UNKs            <- BM[BM$Month == "UNK", ]
#    BM              <- BM[!BM$Month == "UNK", ]
#    UNKs            <- acast(UNKs, Month ~ Year, sum, value.var = "Births")
#    #colSums(UNKs)
#  } 
  
  BM$Month          <- as.integer(BM$Month)
 
  # create Month x Year (cohort) matrix
  Bmat              <- acast(BM[BM$LDB == 1, ], Month ~ Year, sum, value.var = "Births")
  
  cohs              <- sort(union(as.integer(colnames(pop)),as.integer(colnames(Bmat))))
  Ncohs             <- length(cohs)
# distribute births of unknown month proportionally (no need- pointless, just need distribution)
#  if (unks){
#    Bmat.w.UNKs             <- Bmat[, colnames(UNKs), drop = FALSE] # need to keep 2 dims for JPN..
#    Bmat.pdf                <- t(t(Bmat.w.UNKs) / colSums(Bmat.w.UNKs))
#    Bmat[, colnames(UNKs)]  <- Bmat.w.UNKs + t(t(Bmat.pdf) * c(UNKs))
#  }
  
  monthDurations            <- monthDurations[, colnames(Bmat)]
  
  # --------------------------------------------------------------
  # John's formula, per nov 12, 2012 email:        
  f.i                       <- t(t(Bmat) / colSums(Bmat))
  BM                        <- apply(monthDurations, 2, cumsum) # cum number at the end of each month
  b.i                       <- rbind(0, t(t(BM) / colSums(monthDurations)))
  b.bar                     <- colSums(f.i * (b.i[1:12, ] + b.i[2:13, ]) / 2)
  b.bar.full                <- rep(.5, Ncohs)
  names(b.bar.full)         <- cohs
  
  ## name-based substitution
  b.bar.full <- nameBasedRecode(b.bar.full, b.bar)
#  b.bar.full[ names(b.bar)[names(b.bar) %in% names(b.bar.full)] ] <- b.bar[names(b.bar) %in% names(b.bar.full) ]
  
# an AC B.bar matrix
  b.bar.mat                 <- matrix(b.bar.full, 
    nrow = Nages, 
    ncol = Ncohs, 
    byrow = TRUE, 
    dimnames = list(ages, cohs))

  sigmasq                       <- colSums(f.i * ((b.i[1:12, ] ^ 2 + b.i[2:13, ] * b.i[1:12, ] + b.i[2:13, ] ^ 2) / 3)) - b.bar ^ 2
  sigmasq.full                  <- rep(1 / 12, Ncohs) # default uniform distribution
  names(sigmasq.full)           <- cohs
  #sigmasq.full[names(sigmasq)]  <- sigmasq # broken
  sigmasq.full <- nameBasedRecode(sigmasq.full, sigmasq)
  sigmasq.mat                   <- matrix(sigmasq.full, 
                                    nrow = Nages, 
                                    ncol = Ncohs, 
                                    byrow = TRUE, 
                                    dimnames = list(ages, cohs))
 
  # now calculate zU and zL:
  zL            <- (1 - b.bar.mat) / 2 + sigmasq.mat / (2 * (1 - b.bar.mat))
  zL            <- zL[, colnames(dl)]
  zU            <- b.bar.mat / 2 + sigmasq.mat / (2 * b.bar.mat)
  zU            <- zU[, colnames(du)]
  
  # exposure calcs (cohorts NOT cut down to required size)
  Exp           <- pop + zL * dl - zU * du
  
  if (save.bin){
    out.path0   <- file.path(WORKING, "Rbin",  paste0(sex, "cohExposuresRaw.Rdata"))
    save(Exp, file = out.path0)
    # now save mx
    #Sys.chmod( out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }
  # if (test){
  #   return(list(Deaths = dl + du, Exp = Exp))
  # }
  invisible(Exp)
}

