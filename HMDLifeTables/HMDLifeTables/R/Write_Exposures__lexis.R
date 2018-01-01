#' @title \code{Write_Exposures_lexis} a function to prepare and save the file \code{exposures_lexis.txt}.
#'
#' @description This function is essentially a stand-alone that prepares \code{exposures_lexis.txt} straight from binary files. 
#'
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param BINFOLDER the folder name where binary output was written  (not a full path). Default \code{"Rbin"}.
#' @param OPENAGE the desired open age. Default value is 110.
#' @param MPVERSION 5 or 6. Default 5. Here this only affects file headers.
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' 
#' @return function called for its side effect of creating the files \code{exposures_lexis.txt}. No value returned.
#' 
#' @author Carl Boe \email{boe@@demog.berkeley.edu}
#'
#' @importFrom reshape2 acast
#'
#' @export
# Author: boe
###############################################################################
Write_Exposures_lexis <- function(
  WORKING = getwd(),
  STATSFOLDER = "RSTATS",
  BINFOLDER = "Rbin",
  BINPATH = NULL,  # optional override to location of binary files
  OPENAGE = 110,
  MPVERSION , # explicit, no default
  XXX = NULL
  ){
  
  # MatlabRound() is for rounding output, should give same result as matlab, assuming that's important
  # by CB, updated by TR to take digits as arg.
  MatlabRoundFW <- function(x, digits = 0, pad = TRUE, Age = FALSE, totalL = 8){ 
    
    # this 1) rounds, and
    # 2) makes sure that the zeros stay on the end of the number string to the specified number of digits
    if (is.numeric(x)){
      fac     <- rep(10 ^ digits, length(x))
      x       <- sprintf(paste0("%.", digits, "f"), floor(x * fac + sign(x) * 0.5) / fac)
    }
    # strings will potentially vary in length due to integers: adds white space out to the
    # longest string (this is relevant e.g. for dx). ensures digit alignment.
    if (pad){
      maxw    <- max(nchar(x))  
      x       <- sprintf(paste0("%", ifelse(Age, {maxw - 1}, maxw), "s"), x)
      x       <- sprintf(paste0("%-", maxw, "s"), x) # this only had affect if Age == TRUE (dealing with '+')
    }
    # add optional left padding to specify total character width
    x         <- sprintf(paste0("%", totalL, "s"), x)
    x
  }
  
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if ( is.null(BINPATH)){
    BINPATH <- file.path(WORKING, BINFOLDER)
  }
  if( !file.exists(BINPATH) ){
    stop("Bad location of R binary data files")
  }
#--------------------------------------------------------
# read in data: El, Eu, Et  from binary Rbin  
#   t.e0      <- with(local(get(load(b.in.path))), ex[Age == "0"])
#  Iceland,  Deaths (Lexis triangle)       Last modified: 20 Dec 2017,MPv6 (in development)
#  
#  
#  Year   Age   Cohort      Female        Male       Total
#  1838    0      1838      173.87      197.58      371.45
#  1838    0      1837      114.13      136.42      250.55

  bin.path.f.L  <- file.path(BINPATH, "fperExposuresLRaw.Rdata")
  bin.path.m.L  <- file.path(BINPATH, "mperExposuresLRaw.Rdata")
  bin.path.f.U  <- file.path(BINPATH, "fperExposuresURaw.Rdata")
  bin.path.m.U  <- file.path(BINPATH, "mperExposuresURaw.Rdata")
  bin.path.f.UL  <- file.path(BINPATH, "fperExposuresRaw.Rdata")
  bin.path.m.UL  <- file.path(BINPATH, "mperExposuresRaw.Rdata")
  
  if(!file.exists(bin.path.f.L)) stop(paste("Cannot access file: ",bin.path.f.L) )
  if(!file.exists(bin.path.m.L)) stop(paste("Cannot access file: ",bin.path.m.L) )
  if(!file.exists(bin.path.f.U)) stop(paste("Cannot access file: ",bin.path.f.U) )
  if(!file.exists(bin.path.m.U)) stop(paste("Cannot access file: ",bin.path.m.U) )
  if(!file.exists(bin.path.f.UL)) stop(paste("Cannot access file: ",bin.path.f.UL) )
  if(!file.exists(bin.path.m.UL)) stop(paste("Cannot access file: ",bin.path.m.UL) )
  
  Exp.L.f <- with(local(mget(load(bin.path.f.L))), El)
  Exp.L.m <- with(local(mget(load(bin.path.m.L))), El)
  Exp.U.f <- with(local(mget(load(bin.path.f.U))), Eu)
  Exp.U.m <- with(local(mget(load(bin.path.m.U))), Eu)
  Exp.UL.f <- with(local(mget(load(bin.path.f.UL))), Exp)
  Exp.UL.m <- with(local(mget(load(bin.path.m.UL))), Exp)
  
  stopifnot( max( abs(Exp.L.f + Exp.U.f - Exp.UL.f ) ) < 1e-8)  # sum tests L+U is total Lexis square, all conformable
  stopifnot( max( abs(Exp.L.m + Exp.U.m - Exp.UL.m ) ) < 1e-8)
  
  # not sure why the binary matrices were stored without good dimnames, fix here
  # the following function creates Year,Age,Triangle Exposure structures which can be merge()'ed
  
  f.meltAgeYearExposure <- function(x, valname="Exposure", triangle=c("L", "U")){ 
    dimnames(x)<- list(Age=dimnames(x)[[1]], Year=dimnames(x)[[2]])
    x.m <- melt(x, value.name=valname)
    #Cohort <- as.integer(x.m$Year) - as.integer(x.m$Age) - ifelse(triangle=="U", 1, 0)
    #x.m.res <- cbind(x.m[, c(1,2)], Cohort=Cohort, x.m[, -c(1,2)]   ) 
    #colnames(x.m.res)[4] <- valname
    x.m.res <- cbind(x.m, Triangle=triangle)
    return(x.m.res)
  }
  melt.Exp.L.f <- f.meltAgeYearExposure(Exp.L.f, valname="Female", triangle="L")
  melt.Exp.L.m <- f.meltAgeYearExposure(Exp.L.m, valname="Male", triangle="L") 
  melt.Exp.U.f <- f.meltAgeYearExposure(Exp.U.f, valname="Female", triangle="U")
  melt.Exp.U.m <- f.meltAgeYearExposure(Exp.U.m, valname="Male", triangle="U") 
 
  melt.Exp.L <- merge( melt.Exp.L.f, melt.Exp.L.m)
  melt.Exp.L$Total <- melt.Exp.L$Female + melt.Exp.L$Male
  melt.Exp.L$Cohort <- as.integer(melt.Exp.L$Year) - as.integer(melt.Exp.L$Age)
  
  melt.Exp.U <- merge( melt.Exp.U.f, melt.Exp.U.m)
  melt.Exp.U$Total <- melt.Exp.U$Female + melt.Exp.U$Male
  melt.Exp.U$Cohort <- as.integer(melt.Exp.U$Year) - as.integer(melt.Exp.U$Age) - 1
  
  melt.Exp <- rbind(melt.Exp.L,  melt.Exp.U)
  
  melt.Exp <- melt.Exp[ order(melt.Exp$Year, melt.Exp$Age), ]
  rownames(melt.Exp) <- NULL
  
  # LexisDB values are clipped at OPENAGE when first read in, so all subsequent arrays do not go beyond OPENAGE
  # This is a potential problem when OPENAGE is not set high enough since details get lost.  In theory, OPENAGE
  # would only come into play at the lifetable rates stage, and not truncate the death and exposure series.  So,
  # the consequence is that "Lexis triangles" need to be viewed as lifetable triangles and not as LexisDB triangles.
  
  # OPENAGE, e.g. age 110 has U,L entries.  The shape of the open age interval is really RR
  # in this case and so there is only a single UL
  isOpenAge <- (melt.Exp$Age == OPENAGE)
  isLowerTriangle <- (melt.Exp$Triangle =="L")
  isUpperTriangle <- (melt.Exp$Triangle =="U")
  Exp.OpenAgeUL <- melt.Exp[isOpenAge & isLowerTriangle, c("Female", "Male", "Total") ] +    # U+L for OPENAGE
                   melt.Exp[isOpenAge & !isLowerTriangle, c("Female", "Male", "Total") ]
  melt.Exp[isOpenAge & isLowerTriangle, c("Female", "Male", "Total") ] <- Exp.OpenAgeUL      # replace L with U+L
  melt.Exp$Cohort[isOpenAge & isLowerTriangle] <- NA     
  melt.Exp <- melt.Exp[ -seq(nrow(melt.Exp))[isOpenAge & isUpperTriangle] , ]  # drop if OPENAGE and U
  
  # set up paths
  STATS.path     <- file.path(WORKING, STATSFOLDER)
  if (!file.exists(STATS.path)){
    dir.create(STATS.path)
    #system(paste0("chgrp hmdcalc ", STATS.path))
    #Sys.chmod(STATS.path, mode = "2775", use_umask = FALSE)
  }
  write.out.file <- file.path(STATS.path, "Exposures_lexis.txt")
  
  CountryLong    <- country.lookup[country.lookup[,1] == XXX, 2]
  DateMod        <- paste0("\tLast modified: ", format(Sys.time(), "%d %b %Y"), ",")
  # Methods Protocol version
  MPvers         <- ifelse(MPVERSION == 5, " MPv5 (May07)", " MPv6 (Nov17)")
  DataType       <- ",  Exposures (Lexis triangle)"
  
  # special formatting for OPENAGE entries of Age and Cohort
  Age.string <- MatlabRoundFW(melt.Exp$Age, Age = TRUE, totalL = 6)
  Age.string <- ifelse( melt.Exp$Age == OPENAGE, paste0(Age.string,"+"), Age.string) 
  Cohort.string <- MatlabRoundFW(melt.Exp$Cohort, totalL = 9)
  Cohort.string <- ifelse( melt.Exp$Age == OPENAGE, "       .", Cohort.string)
  
  # write it out!
  cat(
    # metadata header
    paste0(CountryLong, DataType, DateMod, MPvers,"\n"),
    # column headers, copied and pasted
    " Year   Age   Cohort      Female        Male       Total",
    # the data, rounded and formatted in place- no tabbing
    paste(MatlabRoundFW( melt.Exp$Year, totalL = 5), 
      Age.string,
      Cohort.string,
      MatlabRoundFW(as.vector(melt.Exp$Female), digits = 2, totalL = 12),
      MatlabRoundFW(as.vector(melt.Exp$Male), digits = 2, totalL = 12),
      MatlabRoundFW(as.vector(melt.Exp$Total), digits = 2, totalL = 12), 
      sep = ""), 
    file = write.out.file, sep = "\n")
  #Sys.chmod(write.out.file, mode = "2775", use_umask = FALSE)
  #system(paste0("chgrp hmdcalc ", write.out.file))
  
}
