# Author: triffe
###############################################################################

# script for comparing output from ltper_AxN() output with output from current HMD
# problem: smoothing will differ...

source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltper_AxN.R")

# function to get data to compare
GetTestControl <- function(XXX = "ISL", sex = "f", A = 1, N = 1){
  # get R-generated lifetable
  test <- ltper_AxN(
    country.folder = paste0("/hdir/0/hmd/HMDWORK/",XXX),
    sex = sex, 
    openage = 110, 
    radix = 1e+05, 
    N = N, 
    abridged = ifelse(A == 1, FALSE, TRUE), 
    round.output = TRUE,
    MPversion = 5,
    write.log = FALSE,
    save.bin = FALSE)
  
  # get current HMD table
  fname <- paste0(sex,"ltper_",A,"x",N,".txt")
  fpath <- file.path("/hdir/0/hmd/HMDWORK",XXX,"STATS",fname)
  control <- read.table(fpath, skip = 2, header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
  
 # slice off '+'
  control$Age <- as.integer(gsub( pattern = "\\+",replacement = "", x=control$Age))
  test$Age <- as.integer(gsub( pattern = "\\+",replacement = "", x=test$Age))
  
  list(test = test, control = control)
}

# difference of mx in smoothed ages
mx.sm.diff80p <- function(test, control){
  if (all(dim(test) == dim(control))){
    test.mx <- reshape2:::acast(test, Age ~ Year, value.var = "mx")
    control.mx <- reshape2:::acast(control, Age ~ Year, value.var = "mx")
    return(test.mx[81:111, ] - control.mx[81:111, ])
  } else {
    stop("dimensions of test and control didn't match- check them")
  }
}
# difference of mx in unsmoothed ages
mx.sm.diffunder80check <- function(test, control){
  if (all(dim(test) == dim(control))){
    test.mx <- reshape2:::acast(test, Age ~ Year, value.var = "mx")
    control.mx <- reshape2:::acast(control, Age ~ Year, value.var = "mx")
    return(test.mx[1:80, ] - control.mx[1:80, ])
  } else {
    stop("dimensions of test and control didn't match- check them")
  }
}

#TC <- GetTestControl("ISL","f",1,1)
#test <- TC$test
#control <- TC$control
# simple check: in theory these should always agree, diffs should be 0
#Diffs <- mx.sm.diffunder80check(test, control)
#all(Diffs == 0)

Make1x1mxcomparison <- function(XXX = "ISL", sex = "f"){
  TC      <- GetTestControl(XXX, sex, 1, 1)
  test    <- TC$test
  control <- TC$control
# now Diffs for 80+, we know these will differ often
# plot as surface to see by how much
  Diffs <- t(mx.sm.diff80p(test, control)) # transposed by image, so we match
# we'd like for absolute differences to be < 5e-5 (5/100000), but things may differ
# primarily because we're using a different optimization technique, but also 
# simply because it's a different statistical packages, and there may be rounding errors,
# implementation differences
# all(abs(Diffs) < 5e-5)
  dmax <- max(abs(Diffs))
# make table for margin histogram 
  Etab <- table(factor(round(Diffs, 6),levels = seq(-5e-5, 5e-5, length = 11)))
# for plotting, we want 0s to appear white, impute NA:
  Diffs[Diffs == 0] <- NA
# reds where test > control, blues for <
  cols <- c("blue",rev(colorRampPalette(RColorBrewer:::brewer.pal(9,"RdBu"), space = "Lab")(10)),"red")
  breaks <- c(-.5, seq(-5e-5, 5e-5, length = 11), .5)
  
  years <- as.integer(rownames(Diffs))
  
#graphics.off()
# dev.new(height=3,width = 10)
  par(mai = c(.3, .5, .7, 1), xpd = TRUE)
  image(
    x = years + .5,
    y = 80.5:110.5,
    Diffs, 
    col = cols, 
    breaks = breaks, 
    asp = 1,
    xlim = c(1750, 2011),
    ylim = c(80, 111),
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE)
  rect(min(years), 80, max(years) + 1, 111, border = gray(.4))
# y axis
  segments(min(years), c(80, 90, 100, 110), min(years) - 1, c(80, 90, 100, 110), col = gray(.4))
  text(min(years), c(80, 90, 100, 110), c(80, 90, 100, 110), pos = 2, cex = .8)
# x axis
  yrt <- c(years,max(years) + 1)[c(years, max(years) + 1) %% 10 == 0]
  segments(yrt, 79, yrt, 80, col = gray(.4))
  text(yrt, 79, yrt, cex = .8, pos = 1)
# legend
  cols2 <- c(cols[1:6],"white",cols[7:12])
  ypos  <- seq(80, 111, length = 14)
  rect(max(years) + 5, ypos[1:(length(ypos) - 1)], max(years) + 10, ypos[2:length(ypos)],
    col = cols2, border = NA)
  rect(max(years) + 5, 80, max(years) + 10, 111, border = gray(.4))
  labs <- breaks
  labs[1] <- "< -5e-5" ; labs[length(labs)] <- "> 5e-5"
  text(max(years) + 10, ypos[1:(length(ypos) - 1)]+diff(ypos)[1]/2, labs, pos = 4, cex = .7)
# title
  Title <- paste(XXX, ifelse(sex == "f", "females,","males,"),"m(x) test - control;  ","Max absolute difference =", round(dmax,5))
  text(1800, 115, Title, pos = 4)
  
# hist of errors, top right
  cols3 <- cols2[2:(length(cols2) - 1)]
  rgt  <- max(years) + 15
  l    <- length(Etab) * 4
  lft  <- rgt - l
  xpos <- seq(lft,rgt,length=length(Etab)+1)
  rect(xpos[1:(length(xpos) - 1)], 115, xpos[2:length(xpos)], 115 + 15 * Etab / sum(Etab),
    col = cols3, border = NA, lwd = .5)
  rect(xpos[5:7],115,xpos[6:8],115 + 15 * Etab[5:7] / sum(Etab), border = gray(.4), lwd = .5)
}

#dev.new(height = 6, width = 10)
# compare all mx, for 80+
## "FRACNP","GBRTENW","ISR","LUX","RUS","TWN" removed due to folder write permissions
#cntries <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE","DNK","DEUTE","DEUTW","DEUTNP",
#  "ESP","EST","FIN","FRATNP","GBR","GBRCENW","GBR_NIR","GBR_NP","GBR_SCO",
#  "HUN","IRL","ISL","ITA","JPN","LTU","LVA","NLD","NOR","POL","PRT",
#  "NZL_MA","NZL_NM","NZL_NP","SWE","USA","SVK","SVN","UKR","FRACNP","GBRTENW","ISR","LUX","RUS","TWN")

#fp <- "/data/commons/triffe/git/HMD_Rlifetables_git/src/algorithms/lifetable.period/lifetable.period.tim/testdata"
#pdf(file = file.path(fp,"Comparemx80plus.pdf"), height = 6, width = 11)
#par(mfrow = c(2,1))
#for (XXX in cntries){
#  Make1x1mxcomparison(XXX,"f")
#  Make1x1mxcomparison(XXX,"m")
#}
#dev.off()

# this just test all option combos to make sure the function actually runs:
# will take about 20 minutes or so:
# XXX 46 countries, 
# sex = m,f, 
# N = 1,5,10 
# abridged, unabridged, 
# v5,v6

# start selecting function to be used directly. 
# These functions call other functions as necessary
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/CleanRstuff.R")
## runs all runable-at-this-time tables: 
## single-sex programs:
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltcoh_AxN.R")
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltper_AxN.R")
## both-sex programs:
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltperBoth_AxN.R")
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltcohBoth_AxN.R")
## output writing program:
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/Write_lt.R")
## period and cohort components objects
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/getCohortComponents.R")
#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/getPeriodComponents.R")

# begin loop XXX <- "ISL" ; abridged = TRUE ; N = 1
#time.start <- Sys.time()
# make sure we're doing a fresh run, and not accidentally using old
# pieces floating around. Once these programs are in use, this would
# be a risky, bas function to use, although we'd likely wouldn't keep
# using directory names that start with R
#CleanRstuff()
# begin country loop
#for (XXX in cntries){
#  cat("\n\n", paste0(XXX, " lifetables running:\n"))
#  cat("\nproducing period data objects..")
#  test <- getPeriodComponents(country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#    sex = "f", openage = 110, save.bin = TRUE)
#  test <- getPeriodComponents(country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#    sex = "m", openage = 110, save.bin = TRUE)
#  cat("\nproducing cohort data objects..")
#  test <- getCohortComponents(country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#    sex = "f", openage = 110, save.bin = TRUE)
#  test <- getCohortComponents(country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#    sex = "m", openage = 110, save.bin = TRUE)
#  for (N in c(1, 5, 10)){
#    for (abridged in c(FALSE, TRUE)){
#      for (sex in c("m","f")){
#        # period tables, single sex 
#        test <- ltper_AxN(
#          country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX),
#          sex = sex, 
#          openage = 110, 
#          radix = 1e+05, 
#          N = N, 
#          abridged = abridged, 
#          MPversion = 5,        # just reproduce for now
#          save.bin = TRUE)
#        cat("\n", paste0(XXX, "_", sex, "ltper_", ifelse(abridged, 5, 1), "x", N))
#        # cohort tables, single sex (checked live if need to be run)
#        test <- ltcoh_AxN(
#          country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX),
#          sex = sex, 
#          openage = 110, 
#          radix = 1e+05,
#          CAGEEXTRP = 90,
#          N = N, 
#          abridged = abridged, 
#          save.bin = TRUE,
#          run.condition = "both",
#          warn.if.not.run = FALSE)
#        if (!is.null(test)){
#          cat("\n", paste0(XXX,"_", sex, "ltcoh_", ifelse(abridged, 5, 1), "x", N))
#        }
#      }
#      # both sex period tables, need to be done once both single sex tables are done
#      test <- ltperBoth_AxN(
#        country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#        openage = 110, 
#        radix = 1e+05, 
#        N = N, 
#        abridged = abridged, 
#        MPversion = 5,        # just reproduce for now
#        use.bin = TRUE,
#        save.bin = TRUE)
#      cat("\n", paste0(XXX, "_bltper_", ifelse(abridged, 5, 1), "x", N))
#      # both sex cohort tables
#      test <- ltcohBoth_AxN(
#        country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX), 
#        openage = 110, 
#        CAGEEXTRP = 90,
#        radix = 1e5, 
#        N = N,                
#        abridged = abridged,     
#        save.bin = TRUE,
#        run.condition = "both",
#        warn.if.not.run = FALSE
#      )
#      if (!is.null(test)){
#        cat("\n", paste0(XXX,"_bltcoh_", ifelse(abridged, 5, 1), "x", N))
#      }
#    }
#  }
#  Write_lt(country.folder = paste0("/hdir/0/hmd/HMDWORK/", XXX))
#  cat("\n------------------------------------------------")
#  cat("\n", paste0(XXX, " lifetables printed to txt in ",XXX,"/RSTATS/"))
#  cat("\n------------------------------------------------\n\n")
#}
#time.end <- Sys.time()
#gc()

#source("/data/commons/triffe/git/HMD_Rlifetables_git/R/RunHMDCountry.R")
## a few cntries at a time,
## begun Thu Aug 23, 5 pm, to see if particular countries throw errors
## it would appears 
#
#}