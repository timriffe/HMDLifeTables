#' @title \code{getPeriodComponents} a function to read in and prepare basic data components for period calculations
#'
#' @description This function modularizes the reading in, staggering, and post-territorial adjusting of the 4 basic (v5) data components of period calculations: pop1, pop2, dl and du, and sticks them in a long-format \code{data.frame} amenable to reshaping in downstream functions. 
#'
#' @details This function is called by \code{ltper_Axn()}, \code{ltperBoth_AxN()} and \code{Exposures_Deaths_Mx_AxN()} when needed, though these functions can also take the output of \code{getPeriodComponents()} as input for in-memory data object passing (faster). The function \code{RunHMDCountry()} uses it this way. Another alternative is to save off the output of this function and read it in using \code{load()}, which is another built-in feature of the above functions. This function calls \code{perTadj()} when pertinent.
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param sex either \code{"m"} or \code{"f"}
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/f_periodComponents.Rdata.Rdata} / \code{Rbin/m_periodComponents.Rdata.Rdata} as well? Appropriate name is derived automatically.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a long \code{data.frame} with 5 columns: \code{"Year"},\code{"Age"}, \code{"dl"}, \code{"du"} and \code{"pop1"} (Jan 1st), \code{"pop2"} (Dec 31st). 
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
getPeriodComponents <- function(
  WORKING = getwd(), 
  sex ,  # no defaults, this is a required parameter
  OPENAGE = 110,
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL){
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
# LDB 
  ldb.path      <- file.path(LDBPATH, paste0(sex, XXX, ".txt"))
  LDBobj        <- read.table(ldb.path, header = TRUE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"), stringsAsFactors=FALSE)
# tadj
  tadj.path     <- file.path(IDBPATH, paste0(XXX, "tadj.txt"))
  if (file.exists(tadj.path)){
    tadj <- read.table(tadj.path, header = TRUE, sep = ",", na.strings = ".", stringsAsFactors=FALSE)
  } else {
    tadj <- NULL
  }
  
# ----------------------------------------------------
# reshape data into matrices
# acast() transforms the data into a matrix with years in the columns and ages in the rows.
# pop1 is Jan 1st; pop2 is Dec 31st; du deaths upper; dl deaths lower
  LY        <- max(LDBobj$Year)
  pop1      <- acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Year != LY, ], Age ~ Year, value.var = "Population")
  if(!missing(tadj)){
    if (!is.null(tadj)){
      mytadj <- try(perTadj(LDBobj = LDBobj, tadj = tadj, sex = sex))
      if (class(mytadj) == "try-error"){
        stop("Check tadj for ", sex, XXX)
      } else {
        LDBobj  <- mytadj
      }
    }
  }
# ----------------------------------------------------
# pop2 = end of year; goes 1 year further and might have tadj
 
  pop2      <- acast(LDBobj[LDBobj$Lexis == 2, ], Age ~ Year, value.var = "Population")
# Deaths >= 0 because LDBobj goes 1 year extra for end of year pop. Deaths entered as -1 for missing
  du        <- acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Year != LY, ], Age ~ Year, value.var = "Deaths")
  dl        <- acast(LDBobj[LDBobj$Lexis == 1 & LDBobj$Year != LY, ], Age ~ Year, value.var = "Deaths")
# ---------------------------------------------------------------------------------
# added for case of BEL (or other potential cases of intervening years with missing data)
  pop1[pop1 == -1]      <- NA
  pop2[pop2 == -1]      <- NA
  dl[dl == -1]          <- NA
  du[du == -1]          <- NA
# Roll everything back to open age at this step
  i.OPENAGE             <- OPENAGE + 1 # used later too (index of open age)
  dl[i.OPENAGE, ]       <- colSums(dl[i.OPENAGE:nrow(dl), ]) # NAs now NOT accounted for on purpose
  du[i.OPENAGE, ]       <- colSums(du[i.OPENAGE:nrow(du), ]) # shouldn't be in LDB except explicitly
  pop1[i.OPENAGE, ]     <- colSums(pop1[i.OPENAGE:nrow(pop1), ])
  pop2[i.OPENAGE, ]     <- colSums(pop2[i.OPENAGE:nrow(pop2), ])
  
  dl        <- dl[1:i.OPENAGE, ]
  du        <- du[1:i.OPENAGE, ]
  pop1      <- pop1[1:i.OPENAGE, ]
# also remove first year for end of year pops (already 1 col wider)
  pop2      <- pop2[1:i.OPENAGE, 2:ncol(pop2)]  
# ----------------------------------------------------
# now we have the 4 items necessary:
  years     <- as.integer(colnames(dl))
  ages      <- as.integer(rownames(dl))
  perComp   <- data.frame(
                          Year = rep(years, each = (OPENAGE + 1)),
                          Age = rep(ages, length(unique(years))),
                          dl = as.vector(dl),
                          du = as.vector(du),
                          pop1 = as.vector(pop1),
                          pop2 = as.vector(pop2)
                          )
    if (save.bin){
      # long file name to mark parameters
      dir.path <- file.path(WORKING, "Rbin")
      if(!file.exists(dir.path)) {
        dir.create(dir.path)
        #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
        #system(paste0("chgrp hmdcalc ", dir.path))
      }
      # prepare file names 4 items to save
      out.name0 <- paste0(sex, "_periodComponents.Rdata")
      out.path0 <- file.path(dir.path, out.name0)
      # saving happens here
      save(perComp, file = out.path0)
      #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", out.path0))
    } 
  invisible(perComp)
}
