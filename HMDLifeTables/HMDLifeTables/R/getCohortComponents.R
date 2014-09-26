#' @title \code{getCohortComponents} a function to read in and prepare basic data components for cohort calculations
#'
#' @description This function modularizes the reading in, staggering, and post-territorial adjusting of the 3 basic (v5) data components of cohort calculations: pop, dl and du, and sticks them in a long-format \code{data.frame} amenable to reshaping in downstream functions. 
#'
#' @details This function is called by \code{Exposures_coh()}, \code{ltcoh_Axn()}, \code{ltcohBoth_AxN()} and \code{cExposures_and_cMx_AxN()} when needed, though these functions can also take the output of \code{getCohortComponents()} as input for in-memory data object passing (faster). The function \code{RunHMDCountry()} uses it this way. Another alternative is to save off the output of this function and read it in using \code{load()}, which is another built-in feature of the above functions. This function calls \code{cohTadj()} when pertinent.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param sex either \code{"m"} or \code{"f"}
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/f_cohortComponents.Rdata.Rdata} / \code{Rbin/m_cohortComponents.Rdata.Rdata} as well? Appropriate name is derived automatically.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a long \code{data.frame} with 5 columns: \code{"Year"},\code{"Age"}, \code{"dl"}, \code{"du"} and \code{"pop"}. 
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' @importFrom compiler cmpfun
#' 
#' @export
# Author: triffe
###############################################################################
getCohortComponents <- function(
  WORKING = getwd(), 
  sex = "f",
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL){
#--------------------------------------------------------
# define internal cohTadj, only used in this function.
  cohTadj <- cmpfun(function(LDBobj, tadj, sex){
      # some tadj files also contain "Rb", limit to "Vx" type
      tadj        <- tadj[tadj$Sex == sex & tadj$Type == "Vx", ]
      # make sure properly ordered (LDBobj sorted prior to calling this)
      tadj        <- tadj[order(tadj$Year, tadj$Age), ]
      # this selects the LDBobj indices that are present in tadj
      tadj.i      <- LDBobj$Age %in% tadj$Age & LDBobj$Year %in% tadj$Year &  LDBobj$Lexis == 2
      # add tadj column, by default does nothing
      LDBobj$tadj <- 1
      # cheap warning, should never come up now
      if (sum(tadj.i) != length(tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year])){
        stop("tadj prob, what country is this?")
      }
      # the denominator selects only the tadj items present in LDBobj (match to tadj.i)
      # matching is guaranteed because tadj and LDBobj have been sorted using the same criteria
      # we take the inverse because we're UNDOING the tadj
      LDBobj$tadj[tadj.i] <- 1 / tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year]
      # now, simply mutliply into Deaths and Population
      LDBobj$Population   <- LDBobj$Population * LDBobj$tadj
      LDBobj$Deaths       <- LDBobj$Deaths * LDBobj$tadj
      # this later step is simply useful for indexing outside this function, which relies on position of NAs
      LDBobj$Population[is.na(LDBobj$Population)] <- 0
      LDBobj$Deaths[is.na(LDBobj$Deaths)] <- 0 
      LDBobj
    })
  
#--------------------------------------------------------
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH <- file.path(WORKING, "InputDB")
  }
#--------------------------------------------------------
# read in data:
# TODO- check for these objects as Rdata, produced by LDB R programs...
# LDB 
  ldb.path      <- file.path(LDBPATH, paste0(sex, XXX, ".txt"))
  LDBobj        <- read.table(ldb.path, header = FALSE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
# tadj
  tadj.path     <- file.path(IDBPATH, paste0(XXX, "tadj.txt"))
  if (file.exists(tadj.path)){
    tadj <- read.table(tadj.path, header = TRUE, sep = ",", na.strings = ".")
  } else {
    tadj <- NULL
  }
  
# if tadj specified, then do it. If not then skip
  if(!missing(tadj)){
    if (!is.null(tadj)){
      mytadj <- try(cohTadj(LDBobj = LDBobj, tadj = tadj, sex = sex))
      if (class(mytadj) == "try-error"){
        stop("Check tadj for ", sex, XXX)
      } else {
        LDBobj  <- mytadj
      }
    }
  }
  
# -------------------------------------------------------
# acast() transforms the data into a matrix with cohorts in the columns and ages in the rows.
# NAs are imputed for incomplete cohorts.
  LY        <- max(LDBobj$Year)
  pop       <- acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Year != LY, ], Age ~ Cohort, value.var = "Population")
  du        <- acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Year != LY, ], Age ~ Cohort, value.var = "Deaths")
  dl        <- acast(LDBobj[LDBobj$Lexis == 1 & LDBobj$Year != LY, ], Age ~ Cohort, value.var = "Deaths")
# --------------------------------------------------------
# added for case of BEL (or other potential cases of intervening years with missing d
  pop[pop == -1]  <- NA
  dl[dl == -1]    <- NA
  du[du == -1]    <- NA
# cohort range not the same: dimension to be equal
  cohs.sclice   <- sort(intersect(intersect(colnames(pop), colnames(dl)), colnames(du)))
  dl        <- dl[, cohs.sclice]
  du        <- du[, cohs.sclice]
  pop       <- pop[, cohs.sclice]
  
# ---------------------------------------------------------
# now we have the 3 items necessary:
  years   <- as.integer(colnames(dl))
  ages    <- as.integer(rownames(dl))
  cohComp <- data.frame(
                        Year = rep(years, each = nrow(dl)),
                        Age = rep(ages, length(unique(years))),
                        dl = as.vector(dl),
                        du = as.vector(du),
                        pop = as.vector(pop)
                        )
  if (save.bin){
    # file name to mark parameters
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
    }
    # prepare file names 4 items to save
    out.name0 <- paste0(sex, "_cohortComponents.Rdata")
    out.path0 <- file.path(dir.path, out.name0)
    # saving happens here
    save(cohComp, file = out.path0)
    #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  } 
  invisible(cohComp)
}
