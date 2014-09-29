
# Author: triffe
###############################################################################

# helper function, called by Exposures_per()- just a matrix utility function:
# takes an AC matrix and changes it to an AP matrix. Lexis is 1 or 2 per LDB conventions.
# no need to use independently (though it's useful), as all arguments are supplied
# in Exposures_per() as necessary.

# where lexis is either 1 or 2:
# 1: lower triangle (TL) or Dec 31st populations
# 2: upper triangle (TU) or Jan 1st populations
# watch out, if unspecified and no attribute, then defaults to 2 with no warning.
# rownames must be ages (no '+'!), and colnames cohorts.

AC2AP <- function(ACmatrix, Lexis){
  if (missing(Lexis)){
    Lexis         <- ifelse(is.null(attr(ACmatrix, "Lexis")), 2, attr(ACmatrix, "Lexis"))
  }
  # back to LDB format, but lacking Year
  longform        <- reshape2:::melt(ACmatrix, varnames = c("Age", "Cohort"), value.name = "value")
  # assume Year is t+x+1. where 1 is the year.offset. Could also be 0. Test if not sure
  longform$Year   <- longform$Cohort + longform$Age + ifelse(Lexis == 2, 1, 0)
  longform        <- longform[!is.na(longform$value), ]
  # cast back to Age x Year matrix
  APmatrix        <- reshape2:::acast(longform, Age ~ Year, value.var = "value")
  attr(APmatrix, "Lexis") <- Lexis
  APmatrix
}

# ------------------------------------------------------------------
# getPeriodComponents(),
# function used to create the four basic data objects, matrices, used to calculate
# period lifetables: pop1, pop2, dl, du. This is called by several functions as a backup.
# Exposures_per can bypass this function if it's 4 output items are passed in as arguments,
# but here you'll need to produce them. This function calls perTadj(), for period tadj
# no need to call this function independenty- its arguments are supplied as needed by
# Exposures_per()

# country.folder is a file path to the folder whose name is an HMD country short name and contains
# InputDB, LexisDB. Other folder are created in there as needed, such as RSTATS, Rbin. Other arguments
# are self explanatory
getPeriodComponents <- function(
  country.folder = "/data/commons/hmd/HMDWORK/ISL", 
  sex = "f",
  openage = 110,
  save.bin = TRUE){
#--------------------------------------------------------
  parts         <- rev(unlist(strsplit(country.folder, split = .Platform$file.sep)))
  ctry          <- parts[!parts == ""][1] 
#--------------------------------------------------------
# read in data:
# LDB 
  ldb.path      <- file.path(country.folder, "LexisDB", paste0(sex, ctry, ".txt"))
  LDBobj        <- read.table(ldb.path, header = FALSE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
# tadj
  tadj.path     <- file.path(country.folder, "InputDB", paste0(ctry, "tadj.txt"))
  if (file.exists(tadj.path)){
    tadj <- read.table(tadj.path, header = TRUE, sep = ",", na.strings = ".")
  } else {
    tadj <- NULL
  }
  
# ----------------------------------------------------
# reshape data into matrices
# acast() transforms the data into a matrix with years in the columns and ages in the rows.
# pop1 is Jan 1st; pop2 is Dec 31st; du deaths upper; dl deaths lower
  pop1      <- reshape2:::acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Deaths >= 0, ], Age ~ Year, value.var = "Population")
  if(!missing(tadj)){
    if (!is.null(tadj)){
      #source("/data/commons/triffe/git/HMD_Rlifetables_git/R/perTadj.R")
      mytadj <- try(perTadj(LDBobj = LDBobj, tadj = tadj, sex = sex))
      if (class(mytadj) == "try-error"){
        stop("Check tadj for ", sex, ctry)
      } else {
        LDBobj  <- mytadj
      }
    }
  }
# ----------------------------------------------------
# pop2 = end of year; goes 1 year further and might have tadj
  pop2      <- reshape2:::acast(LDBobj[LDBobj$Lexis == 2, ], Age ~ Year, value.var = "Population")
# Deaths >= 0 because LDBobj goes 1 year extra for end of year pop. Deaths entered as -1 for missing
  du        <- reshape2:::acast(LDBobj[LDBobj$Lexis == 2 & LDBobj$Deaths >= 0, ], Age ~ Year, value.var = "Deaths")
  dl        <- reshape2:::acast(LDBobj[LDBobj$Lexis == 1 & LDBobj$Deaths >= 0, ], Age ~ Year, value.var = "Deaths")
# ---------------------------------------------------------------------------------
# Roll everything back to open age at this step
  i.openage             <- openage + 1 # used later too (index of open age)
  dl[i.openage, ]       <- colSums(dl[i.openage:nrow(dl), ], na.rm = TRUE)
  du[i.openage, ]       <- colSums(du[i.openage:nrow(du), ], na.rm = TRUE)
  pop1[i.openage, ]     <- colSums(pop1[i.openage:nrow(pop1), ], na.rm = TRUE)
  pop2[i.openage, ]     <- colSums(pop2[i.openage:nrow(pop2), ], na.rm = TRUE)
  
  dl        <- dl[1:i.openage, ]
  du        <- du[1:i.openage, ]
  pop1      <- pop1[1:i.openage, ]
# also remove first year for end of year pops (already 1 col wider)
  pop2      <- pop2[1:i.openage, 2:ncol(pop2)] 
# ----------------------------------------------------
# now we have the 4 items necessary:
  years <- as.integer(colnames(dl))
  ages  <- as.integer(rownames(dl))
  perComp <- data.frame(
    Year = rep(years, each = (openage + 1)),
    Age = rep(ages, length(unique(years))),
    dl = as.vector(dl),
    du = as.vector(du),
    pop1 = as.vector(pop1),
    pop2 = as.vector(pop2)
  )
  if (save.bin){
    # long file name to mark parameters
    dir.path <- file.path(country.folder, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
    }
    # prepare file names 4 items to save
    out.name0 <- paste0(sex, "_periodComponents.Rdata")
    out.path0 <- file.path(dir.path, out.name0)
    # saving happens here
    save(perComp, file = out.path0)
    Sys.chmod(c(dir.path, out.path0), mode = "2775", use_umask = FALSE)
  } 
  invisible(perComp)
}
# ------------------------------------
# perTadj()
# a tadj function used by  getPeriodComponents(). LDBobj is a data.frame of the LDB. 
# tadj is a data.frame of the tadj file from InputDB, if needed. These arguments
# are produced and supplied within getPeriodComponents(), so no need to worry about them.
perTadj <- compiler:::cmpfun(function(LDBobj, tadj, sex){
    # some tadj files also contain "Rb", limit to "Vx" type
    tadj        <- tadj[tadj$Sex == sex & tadj$Type == "Vx",]
    # make sure properly ordered (LDBobj sorted prior to calling this)
    tadj        <- tadj[order(tadj$Year, tadj$Age),]
    # this selects the LDBobj indices that are present in tadj
    tadj.i      <- LDBobj$Age %in% tadj$Age & LDBobj$Year %in% tadj$Year &  LDBobj$Lexis == 2
    # add tadj column, by default does nothing
    LDBobj$tadj <- 1
    # ITA has tadj years that aren't in the data... used in LDB?
    if (sum(tadj.i) != length(tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year])){
      stop("tadj prob, what country is this?")
    }
    # the denominator selects only the tadj items present in LDBobj (match to tadj.i)
    # matching is guaranteed because tadj and LDBobj have been sorted using the same criteria
    # we take the inverse because we're UNDOING the tadj
    LDBobj$tadj[tadj.i] <- 1 / tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year]
    # now, simply mutliply into Population
    LDBobj$Population   <- LDBobj$Population * LDBobj$tadj
    # this later step is simply useful for indexing outside this function, which relies on position of NAs
    LDBobj$Population[is.na(LDBobj$Population)] <- 0 # necessary so we can use NAs to index later
    LDBobj
  })
# ------------------------------------------------------------------------
# Exposures_per()
# this is the function you actually care about.
# arguments:
# 'country.folder' is the HMD shortname folder in the HMD working directory.
# however, the function expects XXXmonthly.txt to be found in the /InputDB/ folder,
# which it isn't yet for any country. For testing, create a copy of a country folder and
# include the XXXmonthly.txt file in the InputDB. These functions expect the standard 
# HMD folder layout.
# the next 4 arguments you can leave as NULL: pop1, pop2, dl, du are all produced
# by getPeriodComponents(), but this function (defined above) is also called internally 
# as a backstop, so can skip these arguments. 'perComp' you can also skip, leave as NULL,
# 'sex' is either "m" or "f"
# openage you leave as the default 110
# save.bin will save a copy of the output matrix (Exp) to the /Rbin/ folder in the country.folder.
# the function also invisibly returns this matrix, so if working in memory you can set it to false.
# MPversion can either be 5 or 6.

# be sure to have installed the packages:
# 'reshape2' (included in base R distributions these days)

# be sure to change the load() path inside
# output is an AP matrix of exposures, columns and rows labeled. 110+ is labeled as 110, not 110+.
Exposures_per <- function(country.folder = "/data/commons/hmd/HMDWORK/JPN", 
  pop1 = NULL,    
  pop2 = NULL,
  dl = NULL,
  du = NULL, 
  perComp = NULL,
  sex = "m", 
  openage = 110, 
  save.bin = TRUE, 
  MPversion = 5, #MPversion = 6
  test = FALSE, # test <- TRUE
  Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
  ){
  # MPversion can only be 5 or 6
  if(!MPversion %in% c(5, 6)){
    stop("\nonly MP versions 5 and 6 are supported at this time\n")
  }
  # --------------------
  # packages and functions we'll be needing:
  require(reshape2)
  #source("/data/commons/triffe/git/HMD_Rlifetables_git/R/getPeriodComponents.R")
  #source("/data/commons/triffe/git/HMD_Rlifetables_git/R/AC2AP.R")
  # --------------------
  # get country short name
  parts         <- rev(unlist(strsplit(country.folder, split = .Platform$file.sep)))
  ctry          <- parts[!parts == ""][1] 
  
  # --------------------
  # this chunk bypassed if pop1, pop2, dl, du supplied as arguments, otherwise derived, loaded or produced...
  if (is.null(pop1) | is.null(pop2) | is.null(dl) | is.null(du)){
    # typical period objects required   # and also allow load() of Rbin if produced previously
    if (is.null(perComp)){
      percomp.path <- file.path(country.folder, "Rbin", paste0(sex, "_periodComponents.Rdata"))
      if (file.exists(percomp.path)){
        perComp <- local(get(load(percomp.path)))
      } else {
        perComp <- getPeriodComponents(
          country.folder = country.folder, 
          sex = sex,
          openage = openage,
          save.bin = TRUE)
      }
    }
    # shape to AP matrices
    pop1          <- acast(perComp, Age ~ Year, sum, value.var = "pop1") 
    pop2          <- acast(perComp, Age ~ Year, sum, value.var = "pop2") 
    dl            <- acast(perComp, Age ~ Year, sum, value.var = "dl")
    du            <- acast(perComp, Age ~ Year, sum, value.var = "du")
  }

  # ---------------------------------------------------------------------------
  # try to get births monthly. If not available, jump to v5, warn test <- TRUE
  births.monthly.root <- ifelse(test, Monthly.folder, {file.path(country.folder, "InputDB")})

  if (MPversion == 6){
    births.monthly.path <- file.path(births.monthly.root, paste0(ctry, "monthly.txt"))
    if (!file.exists(births.monthly.path)){
      MPversion   <- 5
      cat("\nMPversion was given as", MPversion, "but necessary file was missing:\n", births.monthly.path, "\nreverted to MPversion 5 exposures\n")
    }
  }
# old exposures, considerably simpler :-)
  if (MPversion == 5){
    Exp         <- (pop1 + pop2) / 2 + (dl - du) / 6
    # optional save out
    if (save.bin){
      out.path0   <- file.path(country.folder, "Rbin",  paste0(sex, "perExposuresRaw.Rdata"))
      save(Exp, file = out.path0)
      # now save mx
      Sys.chmod(c(dir.path, out.path0), mode = "2775", use_umask = FALSE)
    }
    if (test){
      return(list(Deaths = dl+du, Exp = Exp))
    }
    return(Exp)
  }
  
# -------------------------------------
# set up dims, names used throughout
  yrs     <- as.integer(colnames(pop1))
  cohs    <- (min(yrs) - nrow(pop1) - 1):(max(yrs) + 1)
  Ncohs   <- length(cohs)
  ages    <- 0:openage
  Nages   <- length(ages)
  
# pop1 pop1 is AP
#births.monthly.path <- "/data/commons/triffe/git/HMD_CS/HMDwork/C_CHE/update_7.11.2012/UpdatedTXT/CHEmonthly.txt"
#births.monthly.path <- "/data/commons/triffe/Desktop/RUSmonthly.txt"
# read in CHEmonthly.txt
  BM              <- read.table(births.monthly.path,
                      header = TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE, 
                      na.strings = ".")
# unique(BM$Month)
# remove TOT, convert Month to integer for cast()ing
  BM              <- BM[!BM$Month == "TOT", ]
  
  unks            <- FALSE
  if (any(BM$Month == "UNK")) {
    unks            <- TRUE
    UNKs            <- BM[BM$Month == "UNK", ]
    BM              <- BM[!BM$Month == "UNK", ]
    UNKs            <- acast(UNKs, Month ~ Year, sum, value.var = "Births")
    #colSums(UNKs)
  } 
  
  BM$Month        <- as.integer(BM$Month)
  
# create Month x Year (cohort) matrix
  Bmat            <- acast(BM[BM$LDB == 1, ], Month ~ Year, sum, value.var = "Births")

# distribute births of unknown month proportionally
  if (unks){
    Bmat.w.UNKs   <- Bmat[, colnames(UNKs), drop = FALSE] # need to keep 2 dims for JPN..
    Bmat.pdf      <- t(t(Bmat.w.UNKs) / colSums(Bmat.w.UNKs))
    Bmat[, colnames(UNKs)] <- Bmat.w.UNKs + t(t(Bmat.pdf) * c(UNKs))
  }
# --------------------------------------------------------------
# get proportion of year passed for eah month midpoint
# replace simple with precise date code:
# month.mids                <- seq(1, 23, by = 2) / 24
# 1) Find rough midpoints (integers underneath)
# rough.mids.Date <- seq.Date(as.Date(paste(colnames(Bmat)[1], "01", "15", sep = "-")), 
#                          by = "month", 
#                          length.out = 12 * ncol(Bmat))
# 2) Round up, ceiling(), and down, floor(), and take the difference for exact month durations in days.
# durations       <- matrix(ceiling_date(rough.mids.Date, "month") - floor_date(rough.mids.Date, "month"), nrow = 12)
# --------------------------------------------------------------
# Edit: load() durations matrix, since it's static
# Pavel: change this, as usual...
  durations       <- local(get(load("/data/commons/triffe/git/HMD_Rlifetables_git/R/.durations.Rdata")))
  durations       <- durations[, colnames(Bmat)] # cut down to Bmat years
# --------------------------------------------------------------
# 3) Within each year the cumsum of the month durations, minus Jan (always 31) 
# plus 1/2 of the month durations (to get to day achieved as of midpoint), 
# then divided by days in given year. Yields proportion of year passed as of each month midpoint.
# month.mids.true <- t(t(apply(durations, 2, cumsum) - 31 + durations / 2) / colSums(durations))
  
# --------------------------------------------------------------
# B.bar is a vector, a single value for a each (precise) year. Monthly birth-weighted mean
# proportion of year passed at birth for years with monthly data.
# b.bar                     <- colSums(Bmat * month.mids.true) / colSums(Bmat) # same as John's formula
# pad with .5 for other cohorts lacking monthly births (uniform assumption)
#  b.bar.full                <- rep(.5, Ncohs)
#  names(b.bar.full)         <- cohs
#  b.bar.full[names(b.bar)]  <- b.bar # my fav way to do variable recoding...
## an AC B.bar matrix
#  b.bar.mat                 <- matrix(b.bar.full, 
#                                nrow = Nages, 
#                                ncol = Ncohs, 
#                                byrow = TRUE, 
#                                dimnames = list(ages, cohs))
## convert to AP, very large dimensions, will need to extract
#  b.bar.AP                  <- AC2AP(b.bar.mat, Lexis = 2)
#  
## b.bar is the mean time at birth (proportion of year completed) 
#  b1.bar                    <- B.bar.AP[1:nrow(pop1), colnames(pop1)] # for cohort t-x-1
#  b2.bar                    <- B.bar.AP[1:nrow(pop2), colnames(pop2)] # for cohort t-x
#  
# --------------------------------------------------------------
# now for sigma1, sigma2
# month.mids.mat             <- matrix(month.mids, nrow = nrow(Bmat), ncol = ncol(Bmat))
#  B.bar.small               <- matrix(B.bar, 
#                                nrow = nrow(Bmat), 
#                                ncol = ncol(Bmat), 
#                                byrow = TRUE, 
#                                dimnames = dimnames(Bmat))
## err, a formula for weighted standard deviation that I found on CrossValidated or somewhere...
#  sigma                     <- sqrt(colSums(Bmat * (month.mids.true - B.bar.small) ^ 2) / colSums(Bmat))
#  sigma.full                <- rep(sqrt(1/12), Ncohs) # default uniform distribution
#  names(sigma.full)         <- cohs
#  sigma.full[names(sigma)]  <- sigma
#  sigma.mat                 <- matrix(sigma.full, 
#                                nrow = Nages, 
#                                ncol = Ncohs, 
#                                byrow = TRUE, 
#                                dimnames = list(ages, cohs))
# -------------------------------------------------------------------------------------
  # John's formula, per nov 12, 2012 email:        
  f.i                       <- t(t(Bmat) / colSums(Bmat))
  BM                        <- apply(durations, 2, cumsum)
  b.i                       <- rbind(0, t(t(BM) / colSums(durations)))
  b.bar                     <- colSums(f.i * (b.i[1:12, ] + b.i[2:13, ]) / 2)
  b.bar.full                <- rep(.5, Ncohs)
  names(b.bar.full)         <- cohs
  b.bar.full[names(b.bar)]  <- b.bar # my fav way to do variable recoding...
# an AC B.bar matrix
  b.bar.mat                 <- matrix(b.bar.full, 
    nrow = Nages, 
    ncol = Ncohs, 
    byrow = TRUE, 
    dimnames = list(ages, cohs))
  
  sigmasq                     <- colSums(f.i * ((b.i[1:12, ] ^ 2 + b.i[2:13, ] * b.i[1:12, ] + b.i[2:13, ] ^ 2) / 3)) - b.bar ^ 2
  sigmasq.full                <- rep(1 / 12, Ncohs) # default uniform distribution
  names(sigmasq.full)         <- cohs
  sigmasq.full[names(sigmasq)]  <- sigmasq
  sigmasq.mat                 <- matrix(sigmasq.full, 
                                nrow = Nages, 
                                ncol = Ncohs, 
                                byrow = TRUE, 
                                dimnames = list(ages, cohs))
  
# --------------------------------------------------------------------
# next part stays same after sigma formula determined.
# --------------------------------------------------------------------
# convert to AP, large dimensions, will need to extract relevant years
  b.bar.AP                  <- AC2AP(b.bar.mat, Lexis = 2)
  sigmasq.AP                  <- AC2AP(sigmasq.mat, 2)
# sigma1 is the time at birth for cohort t-x-1. b2.bar for cohort t-x. These are PERIOD matrices, yay
  sigmasq1                    <- sigmasq.AP[1:nrow(pop1), colnames(pop1)]
  sigmasq2                    <- sigmasq.AP[1:nrow(pop2), colnames(pop2)]
# b.bar is the mean time at birth (proportion of year completed) 
  b1.bar                    <- b.bar.AP[1:nrow(pop1), colnames(pop1)] # for cohort t-x-1
  b2.bar                    <- b.bar.AP[1:nrow(pop2), colnames(pop2)] # for cohort t-x
# -----------------------------------------------------------
# now the actual formulas, per the memo and Dima's correction
# l1,l2 changed to s1, s2
  s1            <- 1 - b2.bar                               # insert formulas
  s2            <- (1 - b2.bar) / 2 - sigmasq2 / (2 * (1 - b2.bar))
  u1            <- b1.bar
  u2            <- b1.bar / 2 - sigmasq1 / (2 * b1.bar)
  
  El            <- s1 * pop2 + s2 * dl
  Eu            <- u1 * pop1 - u2 * du

  Exp         <- (pop1 + pop2) / 2 + (dl - du) / 6
# version 6 period exposures:
  Exp           <- El + Eu
  
  # optional save out
  if (save.bin){
    out.path0   <- file.path(country.folder, "Rbin",  paste0(sex, "perExposuresRaw.Rdata"))
    save(Exp, file = out.path0)
    # now save mx
    Sys.chmod(c(dir.path, out.path0), mode = "2775", use_umask = FALSE)
  }
  if (test){
    return(list(Deaths = dl+du, Exp = Exp))
  }
  invisible(Exp)
}


# notice in the above code at times there are chunks commented out that
# were used in the first iteration, which assumed that all months were
# of equal duration. The present code uses exact month midpoints.

# example test call:
Expv5 <- Exposures_per(country.folder = "path/to/pavels/test/folder/XXX", 
  sex = "f",
  save.bin = FALSE,
  MPversion = 5)
Expv6 <- Exposures_per(country.folder = "path/to/pavels/test/folder/XXX", 
  sex = "f",
  save.bin = FALSE,
  MPversion = 6)
#Expv6         <- El + Eu

# previous testing, comparing.
#Expv6.precise <- Expv6
#Expv6.simple <- Expv6
#hist(Expv6.precise - Expv6.simple)
#hist(Expv6.simple - Exp)
## old >= v5 exposures:
#Exp       <- (pop1 + pop2) / 2 + (dl - du) / 6
#colnames(Expv6)
#library(fields)
#library(RColorBrewer)
#library(grDevices)
#colorfun  <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space = "Lab")
#maxdiff   <- max(pretty(Expv6 - Exp, n = 20)) # left off here- work out breaks
#brks      <- seq(-maxdiff, maxdiff, by = 20)
#image.plot(x = 1876:2011 + .5, y = 0:110+.5, t(Expv6 - Exp), 
#  ylim = c(0, 111), xlim = c(1876, 2012), useRaster = TRUE,
#  breaks = brks, col = colorfun(length(brks) - 1), asp = 1)
## makes a big difference for infants
#
##install.packages("lubridate")
#colorfun  <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space = "Lab")
#maxdiff   <- max(pretty(Expv6.precise - Expv6.simple, n = 20)) # left off here- work out breaks
#brks      <- seq(-maxdiff, maxdiff, by = 2)
#image.plot(x = 1876:2011 + .5, y = 0:110+.5, t(Expv6.precise - Expv6.simple), 
#  ylim = c(0, 111), xlim = c(1876, 2012), useRaster = TRUE,
#  breaks = brks, col = colorfun(length(brks) - 1), asp = 1)


