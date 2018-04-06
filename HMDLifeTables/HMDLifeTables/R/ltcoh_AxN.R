#' @title \code{ltcoh_AxN} calculates single sex cohort lifetables in various time intervals given LDB inputs
#' 
#' @description \code{ltcoh_AxN()} is a top-level function intended to be called by users directly. The function uses a simple heuristic, \code{checkRuncoh()}, to determine whether cohort lifetables are to be calculated in the first place, and returns a warning in case insufficient data are available. Thus, the function can in fact be run on each HMD country (except at this time Belgium).
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param cohComp optional. Output from \code{getCohortComponents()}. If \code{cohComp} is already in memory it's faster to supply this argument, otherwise \code{getCohortComponents()} is called internally.
#' @param sex \code{"m"} or \code{"f"}.
#' @param OPENAGE the desired open age. Default value is 110
#' @param CAGEEXTRP earliest age at which almost-extinct cohorts may begin extrapolation. Default 90
#' @param RADIX the lifetable radix (l_0). Default 1e5
#' @param N year interval: 1, 5 or 10. Other intervals would also work in theory.
#' @param abridged logical. Default \code{FALSE}. Should the lifetable by in single ages or abridged ages (0,1,5,10...)
#' @param save.bin logical. Default \code{TRUE}.  Should the output be saved as e.g. \code{Rbin/ltcoh_1x1.Rdata} as well. Appropriate name is derived systematically. In this case both objects are saved separately.
#' @param MPVERSION 5 or 6. Version 5 exposures assume uniformity, v6 allows for non-uniformity across birthdays in the death distribution by taking information from monthly birth distributions. Differences are negligible.
#' @param run.condition \code{"m"}, \code{"f"}, \code{"both"} or \code{"either"}. Default = \code{sex} argument. According to the heuristic in \code{checkRuncoh()}, should the lifetable be computed or not? For example, if at least one male cohort is extinct and \code{run.condition} is \code{"m"}, then it computes. In production, this should be set to \code{"either"} or \code{"both"}, otherwise we could get a border case where males are run but not females.
#' @param warn.if.not.run logical. Default \code{TRUE}. If there isn't enough data (enough extinct cohorts), should the function return a warning?
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @details In the case of multiyear cohorts, the first year is always determined by year modulo \code{N} = 0. The minimum number of years for a cohort to be included for any N > 1 is 2. Thus the first and last cohorts of 5 or 10 year cohort data might be 2, 3 or 4-year cohorts, depending on the start and end dates. 
#' 
#' This function uses \code{checkRuncoh()} to automatically determine whether at least one cohort has been observed from birth until extinction for the given \code{run.condition}. The specific check is whether any cohort observed from birth ever has a population count of 0 (i.e. it could go extinct prior to or after \code{OPENAGE}.
#' 
#' This function calls several other functions, including \code{getCohortComponents()} (which itself calls \code{cohTadj()}), \code{checkRuncoh()}, \code{CDa0()} and \code{Abridge()}. It is either used directly or called by \code{RunHMDCountry()}.
#' 
#' @return a \code{data.frame} of output, unrounded in long (stacked) format. Columns for \code{"Year"}, \code{"Age"}, \code{"mx"}, \code{"qx"}, \code{"ax"}, \code{"lx"}, \code{"dx"}, \code{"Lx"}, \code{"Tx"}, \code{"ex"}
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#' @importFrom compiler cmpfun
#'  
#' @export
# Author: triffe
###############################################################################
ltcoh_AxN <- function(
  WORKING = getwd(), 
  cohComp = NULL,
  sex = "f",
  OPENAGE = 110, 
  CAGEEXTRP = 90,
  RADIX = 1e5,  
  N = 1, 
  abridged = FALSE, 
  save.bin = TRUE,
  MPVERSION ,  #  version must be explicit
  run.condition = sex, # can be this sex, "either" or "both"
  warn.if.not.run = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL){
# ------------------------------------------------------------------------
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
# ------------------------------------------------------------------------
# check to determine whether lifetables to be calculated happens first
# there may be a check sitting in the Rbin folder.
# if it's from today, we use it, if not, generate a new one.
  check.run.path <- file.path(WORKING, "Rbin", "checkRuncoh.Rdata")
  if (file.exists(check.run.path)){
    load(check.run.path)
    today.date   <- format(Sys.time(), "%d %b %Y")
    if (extinct["timestamp"] != today.date){
      extinct    <- checkRuncoh(WORKING, save.bin = TRUE, XXX = XXX, LDBPATH = LDBPATH, IDBPATH = IDBPATH)
    }
  } else {
    extinct      <- checkRuncoh(WORKING, save.bin = TRUE, XXX = XXX, LDBPATH = LDBPATH, IDBPATH = IDBPATH)
  }
 
  if (!extinct[[run.condition]]){
    if (warn.if.not.run){
      cat("\nCohort lifetables not calculated for", XXX, run.condition, "\n" )
    }
    return(NULL)
  }
  
# ------------------------------------------------------------------------
# Time Saver, implemented Aug 17th, 2012
# if the lifetable is abridged, we can save ourselves having to calculate everything
# if a 1x1 Rbinary version has been saved. only need the lx, Tx, ex columns, reshaped 
# to wide.
if (abridged){
  # if a single-age version of the lifetable exists, then we can exploit this:
  # just be careful not to have rounded files sitting around. I'm gonig to remove the option...
  path.name <- paste0(sex, "ltcoh_1x", ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
  path.name <- file.path(WORKING, "Rbin", path.name)
  if (file.exists(path.name)){
    # the object will be called 'output'
    output     <- local(get(load(path.name)))
    output$Age <- as.integer(gsub(output$Age, pattern = "\\+", replacement = ""))
    # now simply extract the 3 items needed for abridging
    lx <- acast(output, Age ~ Year, value.var = "lx")
    Tx <- acast(output, Age ~ Year, value.var = "Tx")
    ex <- acast(output, Age ~ Year, value.var = "ex")
    # call the Abridge() function
    output <- Abridge(lx, Tx, ex, OPENAGE)
    
    # write out, if desired
    if (save.bin){
      # file name to mark parameters
      dir.path <- file.path(WORKING, "Rbin")
      out.name0 <- paste0(sex, "ltcoh_5x",
        ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
      out.path0 <- file.path(dir.path, out.name0)
      
      # saving happens here:
      save(output, file = out.path0)
      # now save mx
      #Sys.chmod(c(dir.path, out.path0), mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
      #system(paste0("chgrp hmdcalc ", out.path0))
      
    }
    return(output)
  }
} 
# again, this if statement only run if table is both to be abridged, and 
# if the 1xN data were already stored as an unrounded .Rdata object in the /Rbin/ directory

#-------------------------------------------------------------------------------  
# define internal functions

# CohortExtender() looks 1:5 spaces to the left of the first NA and take the sum.
# this is the same as a weighted average if used for Dx and Pop
# if we change the extension process to not impute NAs fro 'du'
# then this will need to be changed to divide by 5.
# this is applied row-wise over any cohort column major matrix
# Described pages 42-44
  CohortExtender <- cmpfun(function(x){
      if (any(is.na(x))){
        x[is.na(x)] <- sum(rev(x[!is.na(x)])[1:5]) / 5
      }
      x
    }
  )

# end internal function definitions
# ---------------------------------------------------------------------------------
  if (is.null(cohComp)){
    cohcomp.path <- file.path(WORKING, "Rbin", paste0(sex, "_cohortComponents.Rdata"))
    if (file.exists(cohcomp.path)){ # tip thanks to Charles on SO
      cohComp <- local(get(load(cohcomp.path)))
      
    } else {
    
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
# now reshape to 4 matrices, age x year
  pop       <- acast(cohComp, Age ~ Year, value.var = "pop")
  dl        <- acast(cohComp, Age ~ Year, value.var = "dl")
  du        <- acast(cohComp, Age ~ Year, value.var = "du")
# ---------------------------------------------------------------------------------
# first want to include only cohorts observed from birth
  Exp       <- Exposures_coh(WORKING = WORKING, 
                  pop = pop,    
                  dl = dl,
                  du = du, 
                  sex = sex, 
                  OPENAGE = OPENAGE, 
                  save.bin = TRUE, 
                  MPVERSION = MPVERSION,
                  XXX = NULL,
                  LDBPATH = LDBPATH,
                  IDBPATH = IDBPATH
                )
# 1) limit to cohorts observed from birth:
  c.from.birth <- !is.na(dl[1, ])
  
  pop 			<- pop[, c.from.birth]
  dl 				<- dl[, c.from.birth]
  du 				<- du[, c.from.birth]
  
# ---------------------------------------------------------------------------------
# imputeNAsfrom du      	  # this incidentally is the non-intuitive step
  dl[is.na(du)] 	<- NA     # that was necessary to replicate original
  pop[is.na(du)] 	<- NA     # matlab results, and that is not spelled out 
                            # makes sense after looking at MP diagram for a while
# ---------------------------------------------------------------------------------
# TODO: CAGEEXTRP is not in the MP. Instead the MP suggests that we dynamically determine
#       the x* as the age where 1% of the exposure is above it. Either we need to change
#       the MP or the code. This was a carry-over from matlab, but now would be a good
#       chance to fix it.
  c.keep    <- colSums(!is.na(pop)) > CAGEEXTRP
  pop       <- pop[, c.keep]
  dl        <- dl[, c.keep]
  du        <- du[, c.keep]
  
  Exp       <- Exp[, colnames(pop)] 
# Extend non-extinct cohorts to age 110(+)
  dl 				<- t(apply(dl, 1, CohortExtender))
  du 				<- t(apply(du, 1, CohortExtender))
  pop 			<- t(apply(pop, 1, CohortExtender))
  Exp 			<- t(apply(Exp, 1, CohortExtender)) # should be identical
# this could be important if there is only 1 cohort observed from birth to extinction
  
# ---------------------------------------------------------------------------------
# Roll everything back to open age at this step
  i.openage <- OPENAGE + 1
 
  dl[i.openage, ]      <- colSums(dl[i.openage:nrow(dl),], na.rm = TRUE)
  du[i.openage, ]      <- colSums(du[i.openage:nrow(du),], na.rm = TRUE)
  pop[i.openage, ]     <- colSums(pop[i.openage:nrow(pop),], na.rm = TRUE)
  dl        <- dl[1:i.openage, ]
  du        <- du[1:i.openage, ]
  pop       <- pop[1:i.openage, ]
  
# ---------------------------------------------------------------------------------
# Aggregate years if necessary:
  dl        <- YearAgg(dl, N = N)
  du        <- YearAgg(du, N = N)
  pop       <- YearAgg(pop, N = N)
  Exp       <- YearAgg(Exp, N = N)
# ---------------------------------------------------------------------------------
  #Exp       <- pop + (1 / 3) * (dl - du)                            # Eq 52, again! MPv5
  Dx        <- dl + du    # total deaths
# ---------------------------------------------------------------------------------  
# begin lifetable calculations 
# ---------------------------------------------------------------------------------
# initial mx, qx, ax: these are big matrices
  mx 				<- Dx / Exp																						# Eq 50  MPv5
  qx 				<- Dx / (pop + dl)																		# Eq 70  MPv5 (same in V6)
  # ax 				<- ((1 / 3) * dl + (2 / 3) * du) / Dx  						  # Eq 71  MPv5
  # prefer solve for ax from identity formula- more flexible given version
  ax        <- ((mx + 1) * qx - mx) / (mx * qx)          

  # ---------------------------------------------------------------------------------
# Age 0 steps:
# swap out a0, q0, m0 using procedure described in protocol pg 40-41
  if (MPVERSION == 5){
  m0 		    <- qx[1, ]
  a0 				<- CDa0(m0, sex)																			# coale demeny eq 61,62 MPv5
  for( i in 1:30){		# 30 is overkill																	
    m0 				<- qx[1, ] / ((a0 - 1) *  qx[1, ] + 1)  						# eq 60, solved for m0 MPv5
    a0				<- CDa0(m0, sex)
  }
# replace first row, age 0, with new estimates
  mx[1, ]		<- m0
  
  # optionally return just mx earlier 
  ax[1, ]		<- a0
  }
  if (MPVERSION >= 6){
    ax[1, ] <- AKq02a0(qx[1, ], sex)
    mx[1, ]	<- qx[1, ] / (1 - (1 - ax[1, ]) * qx[1, ])
  }
# End special age 0 stuff
# ---------------------------------------------------------------------------------
# some out-of-protocol but rather necessary fixes
  ax[is.nan(ax)]  <- .5                                            # found in CB code
  qx.1.i          <- zapsmall(qx - 1) == 0                         # this is a machine precision issue
  qx[qx.1.i]      <- 1                                             # similar to CB
# the index of the age of extinction is the smaller of the last age with pop or death or the open age
  extinct.i       <- apply(Dx + Exp, 2, function(x, i.openage){
                               min(max(which(x > 0)), i.openage)
                             }, i.openage = i.openage
                           )   
  qx[cbind(extinct.i, 1:ncol(qx))] <- 1                           # impute 1 at these ages
  ax[i.openage, ]  <- 1 / mx[i.openage, ]                         # reciprocal of MP, same thing
 
# TODO: note to self and CB: why do we assume constant hazard for open age group
#       when imputed? we make it 1/ax. Idea: why not pool all HMD data >1980 for ages
#       105 to 115, say and estimate a rule of thumb for when we have to impute the 
#       open age mx, ax? Here we actually do it the other way around, getting ax from
#       mx, but the assumption is the same.
  
# impute NAs after the 1 (in other measures too)								# described pg 41 MPv5
  qx <- apply(qx, 2, function(x, i.openage){
      i.extinct <-  min(which(x == 1))
      x[1:i.openage > i.extinct] <- NA
      x
    }, i.openage = i.openage
  )
# make agree
  mx[is.na(qx)]   <- NA
  ax[is.na(qx)]   <- NA
  qx[is.na(qx)]   <- NA # otherwise NaN, not a big deal
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
  px 				      <- 1 - qx 																				# Eq 64 MPv5
  px[is.nan(px)]  <- 0
# lx needs to be done columnwise over px, argument 2 refers to the margin.
  lx 				<- apply(px, 2, function(px., RADIX, OPENAGE){ 					# Eq 65 MPv5
      c(RADIX, RADIX * cumprod(px.[1:OPENAGE]))
    }
    , RADIX = RADIX, OPENAGE = OPENAGE)
  rownames(lx)    <- 0:OPENAGE # these got thrown off because l0 imputed.
  lx[is.na(lx)]   <- 0         # want 0s if no person years lived
# multiplying 2 matrices using '*' does the hadamard product in R (elementwise).
  dx 				      <- lx * qx 																				# Eq 66 MPv5
  dx[is.na(dx)]   <- 0
  Lx 				      <- lx - (1 - ax) * dx 														# Eq 67 MPv5
  Lx[i.openage, ]	<- lx[i.openage, ] * ax[i.openage, ]
# we need to do operations on Lx, but taking its NAs to mean 0
  Lx[is.na(Lx)] 	<- 0
# Tx needs to be done columnwise over Lx, argument 2 refers to the margin.
  Tx 				<- apply(Lx, 2, function(Lx., i.openage, OPENAGE){
      c(rev(cumsum(rev(Lx.[1:OPENAGE]))),0) + Lx.[i.openage]		    # Eq 68 MPv5
    }, OPENAGE = OPENAGE, i.openage = i.openage
  )
  rownames(Tx)    <- 0:OPENAGE
  ex 				      <- Tx / lx 	                                      # Eq 69 MPv5
  
# ---------------------------------------------------------------------------------
# combine into output lifetable
  if (abridged){
    output <- Abridge(lx, Tx, ex, OPENAGE)
  } else {
    cohs <- colnames(Tx)
    # sorry, paste0() is a version 2.15 function, upgrade!
    ages <- paste0(0:OPENAGE, c(rep("", OPENAGE), "+"))
      output <- data.frame(Year = rep(cohs, each = 111),
        Age = rep(ages, length(unique(cohs))),
        mx = as.vector(mx),
        qx = as.vector(qx),
        ax = as.vector(ax),
        lx = as.vector(lx),
        dx = as.vector(dx),
        Lx = as.vector(Lx),
        Tx = as.vector(Tx),
        ex = as.vector(ex),
        stringsAsFactors = FALSE)
  }
  if (save.bin){
    # long file name to mark parameters
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
      
    }
    out.name0 <- paste0(sex, "ltcoh",
      ifelse(abridged, "_5","_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    out.path0 <- file.path(dir.path, out.name0)

    # saving happens here:
    save(output, file = out.path0)
    # now save mx
    #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }
  # invisibly return output, so that the console isn't cluttered up
  invisible(output)
} # close lifetable function
