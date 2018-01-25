#' @title \code{Population_A} a function to prepare data for the \code{Population} and \code{Population5} data objects
#' 
#' @description This function is top level. It prepares population counts for males, females and Total, in the typical HMD output format (long, stacked). One peculiarity is that tadj years are stacked, with '-' and '+' appended to the year names. That makes this function slightly harder to integrate into the modularity of the other functions, a job which is not yet complete. For instance, it doesn't use \code{getPeriodExposures()}, but ideally would do so.
#' 
#' @param WORKING path to working directory, which typically ends with the HMD country abbreviation. Default \code{getwd()}.
#' @param abridged logical. Default \code{FALSE}.
#' @param OPENAGE the desired open age. Default value is 110.
#' @param save.bin logical. Default \code{TRUE}. Should the output be saved as e.g. \code{Rbin/Population.Rdata} as well? Appropriate name is derived systematically. 
#' @param XXX the HMD country abbreviation. If left \code{NULL}, this is extracted from \code{WORKING} as the last path part.
#' @param LDBPATH in case the LexisDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "LexisDB")}
#' @param IDBPATH in case the InputDB is not in \code{WORKING} (local testing), the full path to the LexisDB folder. If left as \code{NULL} it is assumed to be \code{file.path(WORKING, "InputDB")}
#' 
#' @return a \code{data.frame} in long (staced) format with \code{"Year"}, \code{"Age"}, \code{"Female"}, \code{"Male"} and \code{"Total"} columns. Values are unrounded.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 acast
#' 
#' @export
# Author: triffe
###############################################################################
Population_A <- function(
  WORKING = getwd(), 
  abridged = FALSE, 
  OPENAGE = 110,
  save.bin = TRUE,
  XXX = NULL,
  LDBPATH = NULL,
  IDBPATH = NULL
){
  # ------------------------------------------------------
  # read in data
  # 'XXX' (country abbrev) is used here and there
  if (is.null(XXX)){
    XXX           <- ExtractXXXfromWORKING(WORKING)
  }
  if (is.null(LDBPATH)){
    LDBPATH       <- file.path(WORKING, "LexisDB")
  }
  if (is.null(IDBPATH)){
    IDBPATH       <- file.path(WORKING, "InputDB")
  }
  # Lexis DB object
  ldb.path.f      <- file.path(LDBPATH, paste0("f", XXX, ".txt"))
  LDBobj.f        <- read.table(ldb.path.f, header = FALSE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  ldb.path.m      <- file.path(LDBPATH, paste0("m", XXX, ".txt"))
  LDBobj.m        <- read.table(ldb.path.m, header = FALSE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  
  # pop1 is Jan 1st; pop2 is Dec 31st; du deaths upper; dl deaths lower
  pop1.f          <- pop.f      <- acast(LDBobj.f[LDBobj.f$Lexis == 2, ], Age ~ Year, value.var = "Population")
  pop1.m          <- pop.m      <- acast(LDBobj.m[LDBobj.m$Lexis == 2, ], Age ~ Year, value.var = "Population")
  
  # territorial adjustment
  tadj.path       <- file.path(IDBPATH, paste0(XXX, "tadj.txt"))
  FY              <- min(LDBobj.f$Year) # bug fix for ITA (tadj in 1st year makes no sense)
  
  if (file.exists(tadj.path)){
    tadj                      <- read.table(tadj.path, header = TRUE, sep = ",", na.strings = ".")
    tadj.years                <- unique(tadj$Year[tadj$Type == "Vx" & tadj$Year > FY])
    LDBobj2.f                 <- try(perTadj(LDBobj = LDBobj.f, tadj = tadj, sex = "f"))
    LDBobj2.m                 <- try(perTadj(LDBobj = LDBobj.m, tadj = tadj, sex = "m"))
    pop2.f                    <- acast(LDBobj2.f[LDBobj2.f$Lexis == 2, ], Age ~ Year, value.var = "Population")
    pop2.m                    <- acast(LDBobj2.m[LDBobj2.m$Lexis == 2, ], Age ~ Year, value.var = "Population") 
    # Luckily order recognizes "-" before "+"..
    # this line necessary for ITA, which has a tadj before pop data actually start...
    tadj.years                <- tadj.years[as.character(tadj.years) %in% colnames(pop.f)]
    # order(c("1994+","1994-"))
    yr.index                  <- colnames(pop.f) %in% tadj.years
    # orig pop is already tadj-adjusted, (it's post = "+")
    colnames(pop.f)[yr.index] <- paste0(tadj.years, "+")
    colnames(pop.m)[yr.index] <- paste0(tadj.years, "+")
    # our tadj above 'undid' the tadj, to make it a Dec 31st population, give it a "-"
    tadj.add.f                <- pop2.f[ , as.character(tadj.years), drop = FALSE]
    tadj.add.m                <- pop2.m[ , as.character(tadj.years), drop = FALSE]
    colnames(tadj.add.f)      <- paste0(tadj.years, "-")
    colnames(tadj.add.m)      <- paste0(tadj.years, "-")
    # stick together order correctly
    pop.f                     <- cbind(pop.f, tadj.add.f)
    pop.m                     <- cbind(pop.m, tadj.add.m)
    pop.f                     <- pop.f[, order(colnames(pop.f))]
    pop.m                     <- pop.m[, order(colnames(pop.m))]
  } else {
    # otherwise still need pop2's
    ## CAB ??  Revisit.  NB: pop2 assignment does not look correct, looks identical to pop1, so if no tadj file, pop2 == pop1
    ## in addition pop2 is not used subsequent to this assignment.  Should just remove this else{} clause
    pop2.f            <- acast(LDBobj.f[LDBobj.f$Lexis == 2, ], Age ~ Year, value.var = "Population")
    pop2.m            <- acast(LDBobj.m[LDBobj.m$Lexis == 2, ], Age ~ Year, value.var = "Population") 
  }
  
  ## NAs are coded as -1 for pop and for deaths, recode
  pop.f[ pop.f == -1] <- NA
  pop.m[ pop.m == -1] <- NA
  
  # sum open age, cut down
  i.OPENAGE           <- OPENAGE + 1
  ## CAB: fix for BEL case of entire missing years of data, in which case i.OPENAGE gets NA as well, rather than 0
  isallNA.f <- apply(pop.f,2, function(x){ all( is.na(x) ) })
  isallNA.m <- apply(pop.m,2, function(x){ all( is.na(x) ) })
  
  pop.f[i.OPENAGE, ]  <- colSums(pop.f[i.OPENAGE:nrow(pop.f), ], na.rm = TRUE)
  pop.m[i.OPENAGE, ]  <- colSums(pop.m[i.OPENAGE:nrow(pop.m), ], na.rm = TRUE)
  
  ## CAB: apply fix
  pop.f[i.OPENAGE, ] <- ifelse(isallNA.f, NA, pop.f[i.OPENAGE, ] )
  pop.m[i.OPENAGE, ] <- ifelse(isallNA.m, NA, pop.m[i.OPENAGE, ] )
  
  pop.f               <- pop.f[1:i.OPENAGE, ]
  pop.m               <- pop.m[1:i.OPENAGE, ]
  
  ages <- c(0:(OPENAGE - 1), paste0(OPENAGE, "+"))
  
  # abridge, if necessary
  if (abridged){
    # fac is a grouping factor for tapply()
    fac <- c(0,1,1,1,1,5:OPENAGE - 5:OPENAGE %% 5)
    pop.f <- apply(pop.f, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    pop.m <- apply(pop.m, 2, function(x, fac){
        tapply(x, fac, sum)
      }, fac = fac 
    )
    # redo ages if necessary: these are width-formatted
    ages <- c("    0  ", 
      paste(sprintf("%3s", c(1, seq(5, OPENAGE - 5, by = 5))),
        sprintf("%-3s", c(4, seq(9, OPENAGE - 1, by = 5))), 
        sep = "-"),
      paste0(OPENAGE, "+   ")
    )        
  }
  
  # total population
  pop.t <- pop.f + pop.m
  years <- colnames(pop.t)
  
  # create output object
  output <- data.frame(
    Year = rep(years, each = length(ages)),
    Age = rep(ages, length(unique(years))),
    Female = as.vector(pop.f),
    Male = as.vector(pop.m),
    Total = as.vector(pop.t),
    stringsAsFactors = FALSE)
  
  # optionally save to R binary format (for instance for post production of both-sex tables)
  if (save.bin){
    # long file name to mark parameters
    out.name <- paste0("Population", ifelse(abridged, "5.Rdata",".Rdata"))
    dir.path <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      #system(paste0("chgrp hmdcalc ", dir.path))
      #Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
    }
    
    out.path <- file.path(dir.path, out.name)
    # saving happens here
    save(output, file = out.path)
    #Sys.chmod(out.path, mode = "2775", use_umask = FALSE)
    #system(paste0("chgrp hmdcalc ", out.path))
  }
  
  invisible(output)
} # end function
