
# Author: triffe
###############################################################################
devtools::load_all("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables", TRUE)


CTRIES <- DemogBerkeley::getHMDcountries()

#  females:
XXX <- "AUT"
a0tests <- do.call(rbind,parallel::mclapply(CTRIES, function(XXX){
      cat(XXX,"\n")
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    v5       <- ltper_AxN(sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5, 
                          testa0 = FALSE, 
                          save.bin = FALSE)
    v5.d       <- ltper_AxN(sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 5, 
                          testa0 = TRUE, 
                          save.bin = FALSE)
    v6         <- ltper_AxN(sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 6, 
                          testa0 = FALSE, 
                          save.bin = FALSE)
    v6.d       <- ltper_AxN(sex = "f", 
                          WORKING = WRKNG, 
                          MPVERSION = 6, 
                          testa0 = TRUE, 
                          save.bin = FALSE)                  
    data.frame(PopName = XXX,
       Sex = "f",
       Year = with(v5,Year[Age == "0"]),
       a0v5 = with(v5,ax[Age == "0"]),
       a0v5.d = with(v5.d,ax[Age == "0"]),
       a0v6 = with(v6,ax[Age == "0"]),
       a0v6.d = with(v6.d,ax[Age == "0"]),
       L0v5 = with(v5,Lx[Age == "0"]),
       L0v5.d = with(v5.d,Lx[Age == "0"]),
       L0v6 = with(v6,Lx[Age == "0"]),
       L0v6.d = with(v6.d,Lx[Age == "0"]),
       stringsAsFactors = FALSE
       )
  }, mc.cores = 8))


hist((a0tests$L0v5- a0tests$L0v6)/1e5)




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
