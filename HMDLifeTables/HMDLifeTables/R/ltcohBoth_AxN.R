#' @title \code{ltcohBoth_AxN} calculates both sex cohort lifetables in various time intervals given LDB inputs
#' 
#' @description \code{ltcohBoth_AxN()} is a top-level function intended to be called by users directly. The function uses a simple heuristic, \code{checkRuncoh()}, to determine whether cohort lifetables are to be calculated in the first place, and returns a warning in case insufficient data are available. Thus, the function can in fact be run on each HMD country (except at this time Belgium).
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param males optional. Output from \code{getCohortComponents()} for males. If \code{cohComp} is already in memory it's faster to supply this argument, otherwise \code{getCohortComponents()} is called internally.
#' #' @param females optional. Output from \code{getCohortComponents()} for females. If \code{cohComp} is already in memory it's faster to supply this argument, otherwise \code{getCohortComponents()} is called internally.
#' @param OPENAGE the desired open age. Default value is 110
#' @param CAGEEXTRP earliest age at which almost-extinct cohorts may begin extrapolation. Default 90
#' @param RADIX the lifetable radix (l_0). Default 1e5
#' @param N year interval: 1, 5 or 10. Other intervals would also work in theory.
#' @param abridged logical. Default \code{FALSE}. Should the lifetable by in single ages or abridged ages (0,1,5,10...)
#' @param save.bin logical. Default \code{TRUE}.  Should the output be saved as e.g. \code{Rbin/ltcoh_1x1.Rdata} as well. Appropriate name is derived systematically. In this case both objects are saved separately.
#' @param MPVERSION 5 or 6. Version 5 exposures assume uniformity, v6 allows for non-uniformity across birthdays in the death distribution by taking information from monthly birth distributions. Differences are negligible.
#' @param run.condition \code{"both"} or \code{"either"}. Default = \code{"either"}. According to the heuristic in \code{checkRuncoh()}, should the lifetable be computed or not? For example, if the condition is either, then if either males or females are extinct, then both sex table will be run. \code{"both"} is slightly more restrictive, though these will only differ in cases such as NXL at the time of this writing. 
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
#' 
#' @export
# Author: triffe
###############################################################################
ltcohBoth_AxN <- function(
  WORKING = getwd(), 
  males = NULL,
  females = NULL,
  OPENAGE = 110, 
  CAGEEXTRP = 90,
  RADIX = 1e5, 
  N = 1,                
  abridged = FALSE,     
  save.bin = TRUE,
  MPVERSION , # explicit, no default
  run.condition = "either", # could also be "both" , for both
  warn.if.not.run = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL
){
# if use.bin TRUE, then read appropriate bin files, if they exist; if not, make them live
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
  if (file.exists(file.path(WORKING, "Rbin", "checkRuncoh.Rdata"))){
    extinct <- local(get(load(file.path(WORKING, "Rbin", "checkRuncoh.Rdata"))))
    today.date <- format(Sys.time(), "%d %b %Y")
    if (extinct[["timestamp"]] != today.date){
      extinct <- checkRuncoh(WORKING, save.bin = TRUE, XXX = XXX, LDBPATH = LDBPATH, IDBPATH = IDBPATH)
    }
  } else {
    extinct <- checkRuncoh(WORKING, save.bin = TRUE, XXX = XXX, LDBPATH = LDBPATH, IDBPATH = IDBPATH)
  }
  if (!extinct[[run.condition]]){
    if (warn.if.not.run){
      cat("\nCohort lifetables not calculated for", XXX, run.condition, "\n" )
    }
    return(NULL)
  }

# an internal function: require(compiler)
  CohortExtender <- cmpfun(function(x){
      if (any(is.na(x))){
        x[is.na(x)] <- sum(rev(x[!is.na(x)])[1:5]) / 5
      }
      x
    }
  )
# ------------------------------------------------------------------------
# Time Saver, implemented Aug 17th, 2012
# if the lifetable is abridged, we can save ourselves having to calculate everything
# if a 1x1 Rbinary version has been saved. only need the lx, Tx, ex columns, reshaped 
# to wide.
  if (abridged){
    # if a single-age version of the lifetable exists, then we can exploit this:
    # just be careful not to have rounded files sitting around. I'm gonig to remove the option...
    path.name <- paste0("bltcoh_1x", ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    path.name <- file.path(WORKING, "Rbin", path.name)
    if (file.exists(path.name)){
      # the object will be called 'output'
      output     <- local(get(load(path.name)))
      output$Age <- as.integer(gsub(output$Age, pattern = "\\+", replacement = ""))
      # now simply extract the 3 items needed for abridging
      lx <- acast(output, Age ~ Year, value.var = "lx")
      Tx <- acast(output, Age ~ Year, value.var = "Tx")
      ex <- acast(output, Age ~ Year, value.var = "ex")
      
      output <- Abridge(lx, Tx, ex, OPENAGE)
      
      # write out, if desired
      if (save.bin){
        # file name to mark parameters
        dir.path <- file.path(WORKING, "Rbin")
        out.name0 <- paste0("bltcoh_5x",
          ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
        out.path0 <- file.path(dir.path, out.name0)
        
        # saving happens here:
        save(output, file = out.path0)
        # now save mx
        #Sys.chmod(c(dir.path, out.path0), mode = "2775", use_umask = FALSE)
      }
      return(output)
    }
  } 
  
# again, this if statement only run if table is both to be abridged, and 
# if the 1xN data were already stored as an unrounded .Rdata object in the /Rbin/ directory
# ------------------------------------------------------------------------  
# otherwise, we can run a fresh table. read in data, or generate if not available:
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
  # require(reshape2)
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
                    save.bin = FALSE, 
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
                    save.bin = FALSE, 
                    MPVERSION = MPVERSION,
                    XXX = NULL,
                    LDBPATH = LDBPATH,
                    IDBPATH = IDBPATH)
 
# first want to include only cohorts observed from birth
  
# 1) limit to cohorts observed from birth:
   c.from.birth <- !is.na(dl.f[1, ])
  
   pop.f     <- pop.f[, c.from.birth]
   dl.f      <- dl.f[, c.from.birth]
   du.f      <- du.f[, c.from.birth]
   pop.m     <- pop.m[, c.from.birth]
   dl.m      <- dl.m[, c.from.birth]
   du.m      <- du.m[, c.from.birth]
# ---------------------------------------------------------------------------------
# imputeNAsfrom du          # this incidentally is the non-intuitive step
  dl.f[is.na(du.f)]   <- NA     # that was necessary to replicate original
  pop.f[is.na(du.f)]  <- NA     # matlab results, and that is not spelled out 
  dl.m[is.na(du.m)]   <- NA     
  pop.m[is.na(du.m)]  <- NA     
  
  # makes sense after looking at MP diagram for a while
# ---------------------------------------------------------------------------------
  c.keep    <- colSums(!is.na(pop.f)) > CAGEEXTRP
  # c.from.birth moved for testing reasons for now
  pop.f     <- pop.f[, c.keep]
  dl.f      <- dl.f[, c.keep]
  du.f      <- du.f[, c.keep]
  pop.m     <- pop.m[, c.keep]
  dl.m      <- dl.m[, c.keep]
  du.m      <- du.m[, c.keep]
  Exp.f     <- Exp.f[, colnames(pop.f)]
  Exp.m     <- Exp.m[, colnames(pop.m)]
# 
# Extend non-extinct cohorts to age 110(+)
  dl.f      <- t(apply(dl.f, 1, CohortExtender))
  du.f      <- t(apply(du.f, 1, CohortExtender))
  pop.f     <- t(apply(pop.f, 1, CohortExtender))
  dl.m      <- t(apply(dl.m, 1, CohortExtender))
  du.m      <- t(apply(du.m, 1, CohortExtender))
  pop.m     <- t(apply(pop.m, 1, CohortExtender))
  Exp.f 	  <- t(apply(Exp.f, 1, CohortExtender)) # should be identical
  Exp.m 		<- t(apply(Exp.m, 1, CohortExtender)) # should be identical
# ---------------------------------------------------------------------------------
# Roll everything back to open age at this step
  i.openage <- OPENAGE + 1
  
  dl.f[i.openage, ]      <- colSums(dl.f[i.openage:nrow(dl.f),], na.rm = TRUE)
  du.f[i.openage, ]      <- colSums(du.f[i.openage:nrow(du.f),], na.rm = TRUE)
  pop.f[i.openage, ]     <- colSums(pop.f[i.openage:nrow(pop.f),], na.rm = TRUE)
  dl.m[i.openage, ]      <- colSums(dl.m[i.openage:nrow(dl.m),], na.rm = TRUE)
  du.m[i.openage, ]      <- colSums(du.m[i.openage:nrow(du.m),], na.rm = TRUE)
  pop.m[i.openage, ]     <- colSums(pop.m[i.openage:nrow(pop.m),], na.rm = TRUE)  
  
  dl.f      <- dl.f[1:i.openage, ]
  du.f      <- du.f[1:i.openage, ]
  pop.f     <- pop.f[1:i.openage, ]
  dl.m      <- dl.m[1:i.openage, ]
  du.m      <- du.m[1:i.openage, ]
  pop.m     <- pop.m[1:i.openage, ]
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
# WEIGHTS NOT USED for cohort both sex tables, just add them together
  dl        <- dl.m + dl.f
  du        <- du.m + du.f
  pop       <- pop.m + pop.f
  
# ---------------------------------------------------------------------------------
# get Ex and Dx
# Exp       <- pop + (1 / 3) * (dl - du)                            # Eq 52 MPv5
  Exp.f[is.na(Exp.f)] <- 0
  Exp.m[is.na(Exp.m)] <- 0
  Exp                 <- Exp.f + Exp.m
  Dx                  <- dl + du    # total deaths
# ---------------------------------------------------------------------------------
# begin lifetable calculations 
# the both-sex primary lifetable elements:
  qx        <- Dx / (pop + dl)                                    # Eq 70  MPv5
  mx        <- Dx / Exp                                           # Eq 50  MPv5
# ---------------------------------------------------------------------------------
# first derive a0 for each sex   
# Dx will be used to weight a0s
  Dx.m                <- dl.m + du.m                    
  Dx.f                <- dl.f + du.f 
  # ------------------------------------------------------------------------------

  # q0 used to start iteration as first m0 in v5, to derive a0 in v6
  q0m                 <- Dx.m[1,] / (pop.m[1,] + dl.m[1,])    
  q0f                 <- Dx.f[1,] / (pop.f[1,] + dl.f[1,])  
  if (MPVERSION == 5){
    m0m              <- q0m
    m0f              <- q0f
    a0m              <- CDa0(m0m, "m")  
    a0f              <- CDa0(m0f, "f")  
    for( i in 1:30){    # 30 is overkill                                  
      m0m            <- q0m / ((a0m - 1) *  q0m + 1)              # eq 60 MPv5, solved for m0
      a0m            <- CDa0(m0m, "m")
      m0f            <- q0f / ((a0f - 1) *  q0f + 1)              # eq 60 MPv5, solved for m0
      a0f            <- CDa0(m0f, "f")
    }
# another similar method that yields different results
# derive D.prime as a product of adjusted m0 and original exposure:
    #Dxm.prime <- m0m * Exp.m[1, ]
    #Dxf.prime <- m0f * Exp.f[1, ]
    #Dx.prime  <- Dxm.prime + Dxf.prime 
    ## rearrive at adjusted both-sex m0 using the same trick in reverse:
    #m0        <- Dx.prime / Exp[1, ]

    #mx[1, ]   <- m0 # insert fancy 2-sex m0
    #qx        <- Dx / (pop + dl)                                    # Eq 70  MPv5
    # if we were doing the Dx.prime version..
    #a0        <- ((m0 + 1)* qx[1, ] - m0) / (m0 * qx[1, ])         # identity between a, q, m, solved for a
    #ax        <- ((1 / 3) * dl + (2 / 3) * du) / Dx                 # Eq 71  MPv5
  }
  if (MPVERSION >= 6){
    a0f <- AKq02a0(q0f, sex = "f")
    a0m <- AKq02a0(q0m, sex = "m")
  }
  # solve using identity to be version-robust:
  ax        <-  ((mx + 1) * qx - mx) / (mx * qx)
  
  # take deathth-weighted average
  ax[1, ]   <- (a0f * Dx.f[1, ] + a0m * Dx.m[1, ]) / (Dx.f[1, ] + Dx.m[1, ]) # Eq 63 MPv5
  mx[1, ]   <- qx[1, ] / ((ax[1, ] - 1) * qx[1, ] + 1)            # Eq 60, solved for m MPv5
# incidentally, the following also seems valid, but gives a different result
# not fully developed, this is test code
# Exp.m       <- pop.m + (1 / 3) * (dl.m - du.m)  
# Exp.f       <- pop.f + (1 / 3) * (dl.f - du.f)  
# m02 <- (m0m * Exp.m[1, ] + m0f * Exp.f[1, ]) / (Exp.m[1, ] + Exp.f[1, ]) # exposure weighted m0's
# years <- as.integer(colnames(Exp.m))
# plot(years, mx[1, ], type = 'l')
# lines(years, m02, col = "red")
# png("SWEbothsexcohortm0.png") # png emailed to John and Carl Nov 30 2012
# plot(years, (m02 - mx[1, ]) / mx[1, ], type = 'l', 
#  main = "Sweden both sex cohort m0\nrelative difference (m02-m01)/m01", ylab = "(m02-m01)/m01")
# dev.off()
# getwd()
# ---------------------------------------------------------------------------------
# End special age 0 stuff
# ---------------------------------------------------------------------------------
# some out-of-protocol but rather necessary fixes
  ax[is.nan(ax)]  <- .5                                            # found in CB code
  qx.1.i          <- zapsmall(qx - 1) == 0                         # this is a machine precision issue
  qx[qx.1.i]      <- 1                                             # similar to CB
  i.openage       <- OPENAGE + 1
# the index of the age of extinction is the smaller of the last age with pop or death or the open age
  extinct.i       <- apply(Dx + Exp, 2, function(x, i.openage){
                               min(max(which(x > 0)), i.openage)
                             }, i.openage = i.openage
                           )   
  qx[cbind(extinct.i, 1:ncol(qx))] <- 1                           # impute 1 at these ages
  ax[i.openage, ] <- 1 / mx[i.openage, ]                          # reciprocal of MPv5, same thing
  
# impute NAs after the 1 (in other measures too)                  # described pg 41 MPv5
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
  px              <- 1 - qx                                         # Eq 64 MPv5
  px[is.nan(px)]  <- 0
# lx needs to be done columnwise over px, argument 2 refers to the margin.
  lx              <- apply(px, 2, function(px., RADIX, OPENAGE){    # Eq 65 MPv5
                          c(RADIX, RADIX * cumprod(px.[1:OPENAGE]))
                        }, RADIX = RADIX, OPENAGE = OPENAGE
                      )
  rownames(lx)    <- 0:OPENAGE # these got thrown off because l0 imputed.
  lx[is.na(lx)]   <- 0         # want 0s if no person years lived
# multiplying 2 matrices using '*' does the hadamard product in R (elementwise).
  dx              <- lx * qx                                        # Eq 66 MPv5
  dx[is.na(dx)]   <- 0
  Lx              <- lx - (1 - ax) * dx                             # Eq 67 MPv5
  Lx[i.openage, ] <- lx[i.openage, ] * ax[i.openage, ]
# we need to do operations on Lx, but taking its NAs to mean 0
  Lx[is.na(Lx)]   <- 0
# Tx needs to be done columnwise over Lx, argument 2 refers to the margin.
  Tx        <- apply(Lx, 2, function(Lx., i.openage, OPENAGE){
      c(rev(cumsum(rev(Lx.[1:OPENAGE]))),0) + Lx.[i.openage]        # Eq 68 MPv5
    }, OPENAGE = OPENAGE, i.openage = i.openage
  )
  rownames(Tx)    <- 0:OPENAGE
  ex              <- Tx / lx                                        # Eq 69 MPv5
  
# ---------------------------------------------------------------------------------
# combine into output lifetable
  if (abridged){
    output <- Abridge(lx, Tx, ex, OPENAGE)
  } else {
    ages   <- paste0(0:OPENAGE, c(rep("", OPENAGE), "+"))
    years  <- colnames(pop)
    output <- data.frame(Year = rep(years, each = (OPENAGE + 1)),
      Age = rep(ages, length(unique(years))),
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
  # ---------------------------------------------------------------------------------
# optionally save to R binary format (for instance for post production of both-sex tables)
  if (save.bin){
    # long file name to mark parameters
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ", dir.path))
    }
    out.name0 <- paste0("bltcoh",
      ifelse(abridged, "_5", "_1"), "x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    
    out.path0 <- file.path(dir.path, out.name0)
    # saving happens here
    save(output, file = out.path0)
    #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }
  invisible(output)
}


