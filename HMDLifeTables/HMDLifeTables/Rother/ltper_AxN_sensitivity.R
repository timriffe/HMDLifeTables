# for sensitivity testing only
# Author: triffe
###############################################################################
ltper_AxN_sensitivity <- function(
  WORKING = getwd(), 
  perComp = NULL,
  sex = "m", 
  OPENAGE = 110, 
  RADIX = 1e5, 
  N = 1,                
  abridged = FALSE,     
  MPVERSION = 5,
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL,
  exposureVersion = 5,
  a0Version = 5){ # complete v5 before going too far with 6.
  
  if (is.null(XXX)){
    XXX          <- ExtractXXXfromWORKING(WORKING) # not sourced!
  }
  if (is.null(LDBPATH)){
    LDBPATH <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH <- file.path(WORKING, "InputDB")
  }
# ------------------------------------------------------  
# Time Saver, implemented Aug 17th, 2012
# if the lifetable is abridged, we can save ourselves having to calculate everything
# if a 1x1 Rbinary version has been saved. only need the lx, Tx, ex columns, reshaped 
# to wide.
  if (abridged){
    # if a single-age version of the lifetable exists, then we can exploit this:
    # just be careful not to have rounded files sitting around. I'm gonig to remove the option...
    path.name <- paste0(sex, "ltper_1x", ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    path.name <- file.path(WORKING, "Rbin", path.name)
    if (file.exists(path.name)){
      # the object will be called 'output'
      output      <- local(get(load(path.name)))
      output$Age  <- as.integer(gsub(output$Age, pattern = "\\+", replacement = ""))
      # now simply extract the 3 items needed for abridging
      lx          <- acast(output, Age ~ Year, value.var = "lx")
      Tx          <- acast(output, Age ~ Year, value.var = "Tx")
      ex          <- acast(output, Age ~ Year, value.var = "ex")
      
      output <- Abridge(lx, Tx, ex, OPENAGE)
      
      # write out, if desired
      if (save.bin){
        # file name to mark parameters
        dir.path <- file.path(WORKING, "Rbin")
        out.name0 <- paste0(sex, "ltper_5x",
          ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
        out.path0 <- file.path(dir.path, out.name0)
        # saving happens here:
        save(output, file = out.path0)
        # now save mx
        #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
        #system(paste0("chgrp hmdcalc ", out.path0))
      }
      return(output)
    }
  } 
# again, this if statement only run if table is both to be abridged, and 
# if the 1xN data were already stored as an unrounded .Rdata object in the /Rbin/ directory
# if the file doesn't exist, then the function calculates everything live.
# ------------------------------------------------------
# now, we check to see if a period Components object is sitting around for us.
# if so, we read it in and use it as such. If not, we call that function and make them
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
# now reshape to 4 matrices:
  pop1      <- acast(perComp, Age ~ Year, value.var = "pop1")
  pop2      <- acast(perComp, Age ~ Year, value.var = "pop2")
  dl        <- acast(perComp, Age ~ Year, value.var = "dl")
  du        <- acast(perComp, Age ~ Year, value.var = "du")
  
# ---------------------------------------------------------------------------------  
# calculate exposure:  
# Exp       <- (pop1 + pop2[, 2:ncol(pop2)]) / 2 + (dl - du) / 6
# Exp       <- (pop1 + pop2) / 2 + (dl - du) / 6                    # Eq 49 MPv5
# Exposure calculations now generalized 
  Exp       <- Exposures_per(WORKING = WORKING, 
    pop1 = pop1,    
    pop2 = pop2,
    dl = dl,
    du = du, 
    sex = sex, 
    OPENAGE = OPENAGE, 
    save.bin = FALSE, 
    MPVERSION = exposureVersion,
    XXX = NULL,
    LDBPATH = LDBPATH,
    IDBPATH = IDBPATH
  )
# ---------------------------------------------------------------------------------
# Aggregate years if necessary: (if N = 1, does nothing)
  dl        <- YearAgg(dl, N = N)
  du        <- YearAgg(du, N = N)
  Exp       <- YearAgg(Exp, N = N)
  
  Dx        <- dl + du                                              # Eq 48 MPv5
  
  # begin lifetable calculations 
  # ---------------------------------------------------------------------------------
  # initial mx, qx, ax: these are big matrices
  # ---------------------------------------------------------------------------------
  # find age at which to begin imputing, based on male (or female) Dx <= 100
  # MP says male, but matlab apparently uses male or female
  # source the function until packaged
  extrap.ages.i <- ltper_GetDx100(WORKING = WORKING, 
    OPENAGE = OPENAGE, 
    N = N, 
    XXX = XXX, 
    LDBPATH = LDBPATH, 
    IDBPATH = IDBPATH)
  # ---------------------------------------------------------------------------------
  # call mx functions (smoothing outsourced to reduce clutter in ltper_AxN() 
  # so far just need Dx, Exp and extrap.ages.i to get back mx
  if (MPVERSION == 5 | MPVERSION == 6){
    mx <- ltper_mx_v5(Dx = Dx, Exp = Exp, extrap.ages.i = extrap.ages.i)
  }
  if (MPVERSION == 7){
    mx <- ltper_mx_v7(Dx = Dx, Exp = Exp, extrap.ages.i = extrap.ages.i)
  }
  
# ---------------------------------------------------------------------------------
# mx now defined, begin the lifetable calcs colnames(mx)
# ---------------------------------------------------------------------------------
  ax        <- Dx * 0 + .5                                          # ax = .5, pg 38 MPv5
  if (a0Version == 5){
    ax[1, ]   <- CDa0(m0 = mx[1, ], sex = sex)                        # Eq 61 / 62 MPv5
  }
  if (a0Version == 6){
    ax[1, ]   <- AKm02a0(m0 = mx[1, ], sex = sex)
  }
#  if (testa0){
#    ax[1, ]   <- AKm02a0_direct(m0 = mx[1, ], sex = sex)
#  }
# multiplying 2 matrices using '*' does the hadamard product in R (elementwise).
  qx        <- mx / (1 + (1 - ax) * mx)                             # Eq 60 MPv5
# ---------------------------------------------------------------------------------
# set open age qx to 1
  i.openage             <- OPENAGE + 1
  qx[i.openage, ]       <- ifelse(is.na(qx[i.openage, ]), NA, 1)
  ax[i.openage, ]       <- 1 / mx[i.openage, ]                   
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
  px 				      <- 1 - qx 																				# Eq 64 MPv5
  px[is.nan(px)]  <- 0 # skips BEL NAs, as these are distinct from NaNs
# lx needs to be done columnwise over px, argument 2 refers to the margin.
  lx 			        <- apply(px, 2, function(px., RADIX, OPENAGE){ 		# Eq 65 MPv5
      if (all(is.na(px.))) {
        px.
      } else {
        c(RADIX, RADIX * cumprod(px.[1:OPENAGE]))
      }
    }, RADIX = RADIX, OPENAGE = OPENAGE
  )
  rownames(lx)    <- 0:OPENAGE # these got thrown off because l0 imputed.
  # NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
  # lx[is.na(lx)]   <- 0 # removed for BEL testing        
  dx 				      <- lx * qx 																				# Eq 66 MPv5
  Lx 				      <- lx - (1 - ax) * dx 														# Eq 67 MPv5
  
  Lx[i.openage, ]	<- lx[i.openage, ] * ax[i.openage, ]
  # we need to do operations on Lx, but taking its NAs to mean 0
  # Lx[is.na(Lx)] 	<- 0 # removed for BEL testing
  # Tx needs to be done columnwise over Lx, argument 2 refers to the column margin.
  Tx 				      <- apply(Lx, 2, function(Lx., i.openage, OPENAGE){
      c(rev(cumsum(rev(Lx.[1:OPENAGE]))),0) + Lx.[i.openage]	# Eq 68 MPv5
    }, OPENAGE = OPENAGE, i.openage = i.openage
  )
  rownames(Tx)    <- 0:OPENAGE
  ex 				      <- Tx / lx 	                                      # Eq 69 MPv5
# ---------------------------------------------------------------------------------
# combine into output lifetable
  if (abridged){
    output <- Abridge(lx, Tx, ex, OPENAGE)
  } else {
    years  <- colnames(Exp)
    ages   <- paste0(0:OPENAGE, c(rep("", OPENAGE), "+"))
    output <- data.frame(
      Year = rep(years, each = (OPENAGE + 1)),
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
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #system(paste0("chgrp hmdcalc ", dir.path))
      #system(paste0("chgrp hmdcalc ", dir.path))
    }
    # prepare file names 4 items to save
    out.name0 <- paste0(sex, "ltper",
      ifelse(abridged, "_5","_1"),"x",
      ifelse(N == 1, "1.Rdata", ifelse(N == 5, "5.Rdata", "10.Rdata")))
    out.path0 <- file.path(dir.path, out.name0)
    # saving happens here
    save(output, file = out.path0)
    #Sys.chmod(out.path0, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path0))
  }
  # invisibly return output, so that the console isn't cluttered up
  invisible(output)
  
} # close lifetable function


