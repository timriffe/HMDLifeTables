
# Author: triffe
###############################################################################
devtools::load_all("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables", TRUE)
source("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables/Rother/ltper_AxN_sensitivity.R")
library(reshape2)
CTRIES <- DemogBerkeley::getHMDcountries()

#  females:
XXX <- "AUT"
a0testsf <- do.call(rbind,parallel::mclapply(CTRIES, function(XXX){
      cat(XXX,"\n")
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    # no v6 exposure differences unless monthly present
    monthlyTF <- any(grepl("monthly",list.files(file.path(WRKNG,"InputDB"))))
    E5A5       <- ltper_AxN_sensitivity(
                          sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5,        # only toggles mx smoothing
                          exposureVersion = 5, 
                          a0Version = 5,
                          save.bin = FALSE)
    E5A6       <- ltper_AxN_sensitivity(
                          sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5,        # only toggles mx smoothing
                          exposureVersion = 5, 
                          a0Version = 6,
                          save.bin = FALSE)       
    E6A5       <- ltper_AxN_sensitivity(
                          sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5,        # only toggles mx smoothing
                          exposureVersion = 6, 
                          a0Version = 5,
                          save.bin = FALSE)
    E6A6       <- ltper_AxN_sensitivity(
                          sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5,        # only toggles mx smoothing
                          exposureVersion = 6, 
                          a0Version = 6,
                          save.bin = FALSE)
    id <- E5A5$Age == 0
    data.frame(
      PopName = XXX, 
      Year = E5A5$Year[id],
      monthly = monthlyTF,
      E5A5a0 = E5A5$ax[id],
      E5A6a0 = E5A6$ax[id],
      E6A5a0 = E6A5$ax[id],
      E6A6a0 = E6A6$ax[id],
      E5A5e0 = E5A5$ex[id],
      E5A6e0 = E5A6$ex[id],
      E6A5e0 = E6A5$ex[id],
      E6A6e0 = E6A6$ex[id],
      stringsAsFactors = FALSE)
  }, mc.cores = 8))
a0testsm <- do.call(rbind,parallel::mclapply(CTRIES, function(XXX){
      cat(XXX,"\n")
      WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
      # no v6 exposure differences unless monthly present
      monthlyTF <- any(grepl("monthly",list.files(file.path(WRKNG,"InputDB"))))
      E5A5       <- ltper_AxN_sensitivity(
        sex = "m", 
        WORKING = WRKNG, 
        MPVERSION = 5,        # only toggles mx smoothing
        exposureVersion = 5, 
        a0Version = 5,
        save.bin = FALSE)
      E5A6       <- ltper_AxN_sensitivity(
        sex = "m", 
        WORKING = WRKNG, 
        MPVERSION = 5,        # only toggles mx smoothing
        exposureVersion = 5, 
        a0Version = 6,
        save.bin = FALSE)       
      E6A5       <- ltper_AxN_sensitivity(
        sex = "m", 
        WORKING = WRKNG, 
        MPVERSION = 5,        # only toggles mx smoothing
        exposureVersion = 6, 
        a0Version = 5,
        save.bin = FALSE)
      E6A6       <- ltper_AxN_sensitivity(
        sex = "m", 
        WORKING = WRKNG, 
        MPVERSION = 5,        # only toggles mx smoothing
        exposureVersion = 6, 
        a0Version = 6,
        save.bin = FALSE)
      id <- E5A5$Age == 0
      data.frame(
        PopName = XXX, 
        Year = E5A5$Year[id],
        monthly = monthlyTF,
        E5A5a0 = E5A5$ax[id],
        E5A6a0 = E5A6$ax[id],
        E6A5a0 = E6A5$ax[id],
        E6A6a0 = E6A6$ax[id],
        E5A5e0 = E5A5$ex[id],
        E5A6e0 = E5A6$ex[id],
        E6A5e0 = E6A5$ex[id],
        E6A6e0 = E6A6$ex[id],
        stringsAsFactors = FALSE)
    }, mc.cores = 8))



    hist((a0testsf$E5A5e0 -  a0testsf$E6A6e0) -
      ((a0testsf$E5A5e0 -  a0testsf$E6A5e0) + (a0testsf$E5A5e0 -  a0testsf$E5A6e0)))
    
  a0testsf$Sex <- "f"
  a0testsm$Sex <- "m"
  a0tests <- rbind(a0testsf,a0testsm)
  setwd("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables/")
  save(a0tests,file="Rother/a0tests.Rdata")
NAind <- !is.na(a0tests$E6A6e0)

  png("Rother/E6A6e0-E5A5e0.png")
  plot(density(a0tests$E6A6e0[NAind] - a0tests$E5A5e0[NAind]), xlim=c(-.4,.4),
    main = "e0, v6 - v5, both exposures and a0 changed")
  dev.off()
  png("Rother/E6A5e0-E5A5e0.png")
  plot(density(a0tests$E6A5e0[NAind]  - a0tests$E5A5e0[NAind]), xlim=c(-.4,.4),
    main = "e0, v6 - v5, only exposures changed")
  dev.off()
  png("Rother/E5A6e0-E5A5e0.png")
  plot(density(a0tests$E5A6e0[NAind] - a0tests$E5A5e0[NAind]), xlim=c(-.4,.4),
    main = "e0, v6 - v5, only a0 changed")
  dev.off()
  
  # -------------------------------------------------------------
# OLD TESTS

#e0diff <- parallel::mclapply(CTRIES, function(XXX){
#    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
#    v5       <- ltper_AxN(sex = "m", WORKING = WRKNG, 
#                          MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
#    v5a0vv6  <- ltper_AxN(sex = "m", WORKING = WRKNG, 
#      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
#    e0diff <- with(v5a0vv6,ex[Age == "0"]) - with(v5,ex[Age == "0"])
#    cbind(x=as.integer(sort(unique(v5$Year))),y=e0diff)
#  }, mc.cores = 4)
#
#names(e0diff) <- CTRIES
#args(png)
#png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/periode0males.png")
#plot(NULL,type = 'n', xlim=c(1750,2012),ylim = c(-.05,.35), 
#  xlab = "year", ylab = "e0 difference (years)", 
#  main = "Change in e0 due to change in a0 method only", sub = "(1 month is about 0.083)")
#abline(h=0,col="red",lwd=2)
#m <- lapply(CTRIES, function(XXX,.e0diff){
#    lines(e0diff[[XXX]], col = "#00000050")
#  },.e0diff = e0diff)
#lines(e0diff[["ISL"]], col = "royalblue")
#lines(e0diff[["NLD"]], col = "green3")
#lines(e0diff[["CHE"]], col = "orange")
#legend("topright",col = c("royalblue","green3","orange"),lty=1,legend=c("Iceland","Netherlands","Switzerland"))
#dev.off()
#CTRIES[order(unlist(lapply(e0diff,function(x){
#      max(x[,2])
#    })))]

#names(e0difff) <- CTRIES
#CTRIES[order(unlist(lapply(e0difff,function(x){
#          max(x[,2])
#        })))]
#png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/periode0females.png")
#plot(NULL,type = 'n', xlim=c(1750,2012),ylim = c(-.05,.35), 
#  xlab = "year", ylab = "e0 difference (years)", 
#  main = "Change in e0 due to change in a0 method only", sub = "(1 month is about 0.083)")
#abline(h=0,col="red",lwd=2)
#m <- lapply(CTRIES, function(XXX,.e0difff){
#    lines(e0difff[[XXX]], col = "#00000050")
#  },.e0difff = e0difff)
#lines(e0difff[["ISL"]], col = "royalblue")
#lines(e0difff[["NLD"]], col = "green3")
#lines(e0difff[["SWE"]], col = "orange")
#legend("topright",col = c("royalblue","green3","orange"),lty=1,legend=c("Iceland","Netherlands","Sweden"))
#dev.off()
#
#
#a0compm <- parallel::mclapply(CTRIES, function(XXX){
#    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
#    v5       <- ltper_AxN(sex = "m", WORKING = WRKNG, 
#      MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
#    v5a0vv6  <- ltper_AxN(sex = "m", WORKING = WRKNG, 
#      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
#    
#    cbind(x=as.integer(sort(unique(v5$Year))),
#      a0v5=with(v5,ax[Age == "0"]),
#      a0v6=with(v5a0vv6,ax[Age == "0"]))
#  }, mc.cores = 4)
#names(a0compm) <- CTRIES
#a0compf <- parallel::mclapply(CTRIES, function(XXX){
#    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
#    v5       <- ltper_AxN(sex = "f", WORKING = WRKNG, 
#      MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
#    v5a0vv6  <- ltper_AxN(sex = "f", WORKING = WRKNG, 
#      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
#    
#    cbind(x=as.integer(sort(unique(v5$Year))),
#      a0v5=with(v5,ax[Age == "0"]),
#      a0v6=with(v5a0vv6,ax[Age == "0"]))
#  }, mc.cores = 4)
#names(a0compf) <- CTRIES
#
#
#png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/Swedena0v5v6.png")
#plot(a0compm[["SWE"]][,1],a0compm[["SWE"]][,2],type = "l", col = "blue",lty=2, ylim = c(0,.35),
#  main= "Sweden a0 v5 'Coale-Demeny' versus v6 'Andreev-Kingkade'", ylab = "a0",xlab = "Year")
#lines(a0compm[["SWE"]][,1],a0compm[["SWE"]][,3], col = "blue")
#lines(a0compf[["SWE"]][,1],a0compf[["SWE"]][,2], col = "red",lty=2)
#lines(a0compf[["SWE"]][,1],a0compf[["SWE"]][,3], col = "red")
#legend("bottomleft",col = c("red","red","blue","blue"),lty=c(1,2,1,2),
#  legend=c("v6 females","v5 females","v6 males","v5 males"))
#dev.off()
