#' @title \code{ltperBoth_AxN} calculates both-sex period lifetables
#' in various time intervals given LDB inputs
#' 
#' @description \code{ltperBoth_AxN()} is a top-level function
#' intended to be called by users directly. All HMD countries
#' calculate period lifetables. This function has not yet been adapted
#' for Belgium.
#' 
#' @param WORKING path to working directory, which typically ends with
#' the HMD country abbreviation. Default \code{getwd()}.
#' @param males optional. Output from \code{getPeriodComponents()} for
#' males. If \code{males} is already in memory it's faster to supply
#' this argument, otherwise \code{getPeriodComponents()} is called
#' internally.
#' @param females optional. Output from \code{getPeriodComponents()}
#' for females. If \code{females} is already in memory it's faster to
#' supply this argument, otherwise \code{getPeriodComponents()} is
#' called internally.
#' @param sex \code{"m"} or \code{"f"}.
#' @param OPENAGE the desired open age. Default value is 110
#' @param RADIX the lifetable radix (l_0). Default 1e5
#' @param N year interval: 1, 5 or 10. Other intervals would also work in theory.
#' @param abridged logical. Default \code{FALSE}. Should the lifetable
#' by in single ages or abridged ages (0,1,5,10, ...)
#' @param save.bin logical. Default \code{TRUE}.  Should the output be
#' saved as e.g. \code{Rbin/ltper_1x1.Rdata} as well. Appropriate name
#' is derived systematically. In this case both objects are saved
#' separately.
#' @param MPVERSION 5 or 6 (or 7). Version 5 exposures assume
#' uniformity, v6 allows for non-uniformity across birthdays in the
#' death distribution by taking information from monthly birth
#' distributions. Differences are negligible. Version 7 uses a
#' different mx smoothing technique.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this
#' is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH LDBPATH in case the LexisDB is not in \code{WORKING}
#' (local testing), the full path to the LexisDB folder. If left as
#' \code{NULL} it is assumed to be \code{file.path(WORKING,
#' "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local
#' testing), the full path to the LexisDB folder. If left as
#' \code{NULL} it is assumed to be \code{file.path(WORKING,
#' "InputDB")}
#' 
#' @details In the case of multiyear cohorts, the first year is always
#' determined by year modulo \code{N} = 0. The minimum number of years
#' for a cohort to be included for any N > 1 is 2. Thus the first and
#' last cohorts of 5 or 10 year cohort data might be 2, 3 or 4-year
#' cohorts, depending on the start and end dates.
#' 
#' This function calls several functions, including
#' \code{getPeriodComponents()}, \code{ltper_getDx100()},
#' \code{ltper_mx_v5()} or \code{ltper_mx_v6()}, \code{Abrdige()},
#' \code{CDa0()} and indirectly \code{perTadj()}. It is either used
#' directly or called by \code{RunHMDCountry()}.
#' 
#' @return a \code{data.frame} of output, unrounded in long (stacked)
#' format. Columns for \code{"Year"}, \code{"Age"}, \code{"mx"},
#' \code{"qx"}, \code{"ax"}, \code{"lx"}, \code{"dx"}, \code{"Lx"},
#' \code{"Tx"}, \code{"ex"}
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 acast 
#' 
#' @export

ltperBoth_AxN <- function(
  WORKING = getwd(), 
  males = NULL,
  females = NULL,
  OPENAGE = 110, 
  RADIX = 1e5, 
  N = 1,                
  abridged = FALSE,     
  MPVERSION = 5,
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
  if (is.null(males) | is.null(females)){
    # if use.bin TRUE, then read appropriate bin files, if they exist; if not, make them live
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
# ---------------------------------------------------------------------------------  
# calculate exposure  
# Exp.f       <- (pop1.f + pop2.f) / 2 + (dl.f - du.f) / 6                    # Eq 49
# Exp.m       <- (pop1.m + pop2.m) / 2 + (dl.m - du.m) / 6     
  Exp.f       <- Exposures_per(WORKING = WORKING, 
                                pop1 = pop1.f,    
                                pop2 = pop2.f,
                                dl = dl.f,
                                du = du.f, 
                                sex = sex, 
                                OPENAGE = OPENAGE, 
                                save.bin = FALSE, 
                                MPVERSION = MPVERSION,
                                XXX = XXX,
                                LDBPATH = LDBPATH,
                                IDBPATH = IDBPATH
                               )
  Exp.m       <- Exposures_per(WORKING = WORKING, 
                                 pop1 = pop1.m,    
                                 pop2 = pop2.m,
                                 dl = dl.m,
                                 du = du.m, 
                                 sex = sex, 
                                 OPENAGE = OPENAGE, 
                                 save.bin = FALSE, 
                                 MPVERSION = MPVERSION,
                                 XXX = XXX,
                                 LDBPATH = LDBPATH,
                                 IDBPATH = IDBPATH
                               )
# ---------------------------------------------------------------------------------  
# Aggregate years if necessary: (if N = 1, does nothing)
  Exp.f       <- YearAgg(Exp.f, N = N)
  Exp.m       <- YearAgg(Exp.m, N = N)
  du.f        <- YearAgg(du.f, N = N)
  du.m        <- YearAgg(du.m, N = N)
  dl.f        <- YearAgg(dl.f, N = N)
  dl.m        <- YearAgg(dl.m, N = N)
# ---------------------------------------------------------------------------------  
# deaths
  Dx.f        <- du.f + dl.f
  Dx.m        <- du.m + dl.m
# ---------------------------------------------------------------------------------
# also need to read in lifetables for smoothed mx (no point re-smoothing)
  lt.path.f <- file.path(WORKING, "Rbin", paste0("fltper_1x", N, ".Rdata"))
  lt.path.m <- file.path(WORKING, "Rbin", paste0("mltper_1x", N, ".Rdata"))
  if (file.exists(lt.path.f) & file.exists(lt.path.m)){
    females <- local(get(load(lt.path.f)))
    males   <- local(get(load(lt.path.m)))
  } else {
    # in case these were not pre-calculated, we grab the data live
    females <- ltper_AxN(
                  WORKING = WORKING, 
                  sex = "f", 
                  OPENAGE = OPENAGE, 
                  RADIX = RADIX, 
                  N = N,                
                  abridged = FALSE,     
                  MPVERSION = MPVERSION,
                  save.bin = TRUE,
                  XXX = XXX,
                  LDBPATH = LDBPATH,
                  IDBPATH = IDBPATH)
    males <- ltper_AxN(
                  WORKING = WORKING, 
                  sex = "m", 
                  OPENAGE = OPENAGE, 
                  RADIX = RADIX, 
                  N = N,                
                  abridged = FALSE,     
                  MPVERSION = MPVERSION,
                  save.bin = TRUE,
                  XXX = XXX,
                  LDBPATH = LDBPATH,
                  IDBPATH = IDBPATH)
  }
  nyrs        <- length(unique(males$Year))
  mx.m        <- matrix(males$mx, ncol = nyrs)
  mx.f        <- matrix(females$mx, ncol = nyrs)

# ---------------------------------------------------------------------------------
# raw weights (proportion female of exposure)
  Exp.t        <- Exp.f + Exp.m
  w.f          <- Exp.f / Exp.t                               # Eq 56 
# ---------------------------------------------------------------------------------
# smooth weights
# age, for model fitting
  age.mid      <- 0:OPENAGE + .5
  i.openage    <- OPENAGE + 1
  years        <- colnames(Exp.f)
  # container to put 'smoothed' weights, pi
  pi.mat       <- w.f # replaced for MPv5, partly used for MPv6+
  
  # needed in either case. Same index as that used for old-age smoothing
  extrap.ages.i <- ltper_GetDx100(WORKING = WORKING, 
    OPENAGE = OPENAGE, 
    N = N, 
    XXX = XXX, 
    LDBPATH = LDBPATH, 
    IDBPATH = IDBPATH)
  
  if (MPVERSION == 5){
    # loop over years
    eps <- 1e-8 # eps from matlab
    for (i in 1:ncol (w.f)){ #
      if (!all(is.na(Exp.f[, i]))){ # clause added for BEL
        keep.i       <- !(Exp.f[, i] <= eps | Exp.m[, i] <= eps)
        #reg.weights  <- Exp.t[keep.i, i]
        w.f.i        <- w.f[keep.i, i]
        age.mid.i    <- age.mid[keep.i]
        #coefs.i      <- lm(log(w.f.i / (1 - w.f.i)) ~ age.mid.i + I(age.mid.i ^ 2), weights = reg.weights)$coef # Eq 57
        coefs.i      <- lm(log(w.f.i / (1 - w.f.i)) ~ age.mid.i + I(age.mid.i ^ 2))$coef # Eq 57
        
        z.i          <- coefs.i[1] + (coefs.i[2] * age.mid) + (coefs.i[3] * age.mid ^ 2) # Eq 58
##        pi.mat[, i]  <- exp(z.i) / (1 + exp(z.i))
        pi.mat[, i]  <- ifelse( age.mid < (extrap.ages.i[i] - 1), w.f[, i], exp(z.i) / (1 + exp(z.i)) ) #use obs values below extrap age
      }
    }
    # i.e. in the matlab code ALL ages are used to fit and ALL ages are blended together thusly
  }
  if (MPVERSION > 5){
    # first get raw:
    D.all      <- dl.f + dl.m + du.f + du.m
    mx         <-  D.all / Exp.t
    # loop over years
    eps <- 1e-8 # eps from matlab
    for (i in 1:ncol (w.f)){ #
      if (!all(is.na(Exp.f[, i]))){ # clause added for BEL
        # different from v5 because only ages 80+ are used to fit
        sm.ind       <- extrap.ages.i[i]:i.openage 
        keep.i       <- !(Exp.f[, i] <= eps | Exp.m[, i] <= eps) & (0:110) >= 80
        w.f.i        <- w.f[keep.i, i]
        age.mid.i    <- age.mid[keep.i]
        # also, approximate weights are used (log() would mean we need to discard exposures <= 1); CAB: change Exp.f to Exp.t
        coefs.i      <- lm(log(w.f.i / (1 - w.f.i)) ~ age.mid.i + I(age.mid.i ^ 2), weights = sqrt(Exp.t[keep.i, i]))$coef # Eq 57
        z.i          <- coefs.i[1] + (coefs.i[2] * age.mid[sm.ind]) + 
          (coefs.i[3] * age.mid[sm.ind] ^ 2) # Eq 58
        pi.vec       <- exp(z.i) / (1 + exp(z.i))   
        #
        pi.mat[sm.ind, i]  <- exp(z.i) / (1 + exp(z.i))
      }
    }
  }
  
  # pi.mat will vary depending on which MPVERSION is in use.
  mx        <- pi.mat * mx.f + (1 - pi.mat) * mx.m

# ---------------------------------------------------------------------------------
# mx now defined, begin the lifetable calcs colnames(mx)
# ---------------------------------------------------------------------------------
  ax        <- mx * 0 + .5                                          # ax = .5, pg 38
  if (MPVERSION == 5){
    a0f       <- CDa0(m0 = mx.f[1, ], sex = "f")                      # Eq 61 / 62
    a0m       <- CDa0(m0 = mx.m[1, ], sex = "m")
  }
  if (MPVERSION >= 6){
    a0f       <- AKm02a0(m0 = mx.f[1, ], sex = "f")                 
    a0m       <- AKm02a0(m0 = mx.m[1, ], sex = "m")
  }
  d0f       <- Dx.f[1, ]
  d0m       <- Dx.m[1, ]
  ax[1, ]   <- (a0f * d0f + a0m * d0m) / (d0f + d0m)                # Eq 63 
  qx        <- mx / (1 + (1 - ax) * mx)                             # Eq 60 
# ---------------------------------------------------------------------------------
# set open age qx to 1
  qx[i.openage, ] <- ifelse(is.na(qx[i.openage, ]), NA, 1)
  ax[i.openage, ] <- 1 / mx[i.openage, ]                            # reciprocal of MP, same thing
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
  px              <- 1 - qx                                         # Eq 64
  px[is.nan(px)]  <- 0 # skips BEL NAs, as these are distinct from NaNs
# lx needs to be done columnwise over px, argument 2 refers to the margin.
  lx 			        <- apply(px, 2, function(px., RADIX, OPENAGE){ 		# Eq 65 MPv5
                        if (all(is.na(px.))) { # added for BEL
                          px.
                        } else {
                          c(RADIX, RADIX * cumprod(px.[1:OPENAGE]))
                        }
                      }, RADIX = RADIX, OPENAGE = OPENAGE
                    )
  rownames(lx)    <- 0:OPENAGE # these got thrown off because l0 imputed.
  # NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
  # lx[is.na(lx)]   <- 0   # removed for BEL testing         
  dx              <- lx * qx                                        # Eq 66
  Lx              <- lx - (1 - ax) * dx                             # Eq 67
  Lx[i.openage, ] <- lx[i.openage, ] * ax[i.openage, ]
  # we need to do operations on Lx, but taking its NAs to mean 0
  # Lx[is.na(Lx)]   <- 0   # removed for BEL testing  
  # Tx needs to be done columnwise over Lx, argument 2 refers to the column margin.
  Tx              <- apply(Lx, 2, function(Lx., i.openage, OPENAGE){
                               c(rev(cumsum(rev(Lx.[1:OPENAGE]))), 0) + Lx.[i.openage]  # Eq 68
                             }, OPENAGE = OPENAGE, i.openage = i.openage
                           )
  rownames(Tx)    <- 0:OPENAGE
  ex              <- Tx / lx                                        # Eq 69 
  
# ---------------------------------------------------------------------------------
# combine into output lifetable
  if (abridged){
    output <- Abridge(lx, Tx, ex, OPENAGE)
  } else {
    # years defined earlier
    ages   <- paste0(0:OPENAGE, c(rep("", OPENAGE), "+"))
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
    # file name to mark parameters
    out.name <- paste0("bltper",
      ifelse(abridged, "_5","_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      #system(paste0("chgrp hmdcalc ",dir.path))
    }
    #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
    
    out.path <- file.path(dir.path, out.name)
    # saving happens here
    save(output, file = out.path)
    #Sys.chmod(out.path, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path))
  }
  invisible(output)
}

