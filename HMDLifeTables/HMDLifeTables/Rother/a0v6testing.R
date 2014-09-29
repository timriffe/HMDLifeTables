
# Author: triffe
###############################################################################
devtools::load_all("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables", TRUE)

args(ltper_AxN)
setwd("/data/commons/hmd/HMDWORK/SWE")
SWEv5    <- ltper_AxN(sex = "m", WORKING = "/data/wilmoth0/HMD/HMDWORK/SWE", 
  MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
SWEv5a0v6 <- ltper_AxN(sex = "m", WORKING = "/data/wilmoth0/HMD/HMDWORK/SWE", 
  MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
head(SWEv5a0v6)
head(SWEv5)
years <- as.integer(sort(unique(SWEv5a0v6$Year)))
plot(years, with(SWEv5,ax[Age == "0"]))
lines(years, with(SWEv5a0v6,ax[Age == "0"]))

plot(years, with(SWEv5,ex[Age == "0"]))
lines(years, with(SWEv5a0v6,ex[Age == "0"]))

plot(years, with(SWEv5a0v6,ex[Age == "0"]) - with(SWEv5,ex[Age == "0"]), type = 'l')

CTRIES <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE","DNK","DEUTE","DEUTW","DEUTNP",
  "ESP","EST","FIN","FRATNP","GBR","GBRCENW","GBR_NIR","GBR_NP","GBR_SCO",
  "HUN","IRL","ISL","ITA","JPN","LTU","LVA","NLD","NOR","POL","PRT",
  "NZL_MA","NZL_NM","NZL_NP","SWE","USA","SVK","SVN","UKR","FRACNP","GBRTENW","ISR","LUX","RUS","TWN")

e0diff <- parallel::mclapply(CTRIES, function(XXX){
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    v5       <- ltper_AxN(sex = "m", WORKING = WRKNG, 
                          MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
    v5a0vv6  <- ltper_AxN(sex = "m", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
    e0diff <- with(v5a0vv6,ex[Age == "0"]) - with(v5,ex[Age == "0"])
    cbind(x=as.integer(sort(unique(v5$Year))),y=e0diff)
  }, mc.cores = 4)

names(e0diff) <- CTRIES
args(png)
png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/periode0males.png")
plot(NULL,type = 'n', xlim=c(1750,2012),ylim = c(-.05,.35), 
  xlab = "year", ylab = "e0 difference (years)", 
  main = "Change in e0 due to change in a0 method only", sub = "(1 month is about 0.083)")
abline(h=0,col="red",lwd=2)
m <- lapply(CTRIES, function(XXX,.e0diff){
    lines(e0diff[[XXX]], col = "#00000050")
  },.e0diff = e0diff)
lines(e0diff[["ISL"]], col = "royalblue")
lines(e0diff[["NLD"]], col = "green3")
lines(e0diff[["CHE"]], col = "orange")
legend("topright",col = c("royalblue","green3","orange"),lty=1,legend=c("Iceland","Netherlands","Switzerland"))
dev.off()
CTRIES[order(unlist(lapply(e0diff,function(x){
      max(x[,2])
    })))]




#  females:

e0difff <- parallel::mclapply(CTRIES, function(XXX){
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    v5       <- ltper_AxN(sex = "f", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
    v5a0vv6  <- ltper_AxN(sex = "f", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
    e0diff <- with(v5a0vv6,ex[Age == "0"]) - with(v5,ex[Age == "0"])
    cbind(x=as.integer(sort(unique(v5$Year))),y=e0diff)
  }, mc.cores = 4)

names(e0difff) <- CTRIES
CTRIES[order(unlist(lapply(e0difff,function(x){
          max(x[,2])
        })))]
png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/periode0females.png")
plot(NULL,type = 'n', xlim=c(1750,2012),ylim = c(-.05,.35), 
  xlab = "year", ylab = "e0 difference (years)", 
  main = "Change in e0 due to change in a0 method only", sub = "(1 month is about 0.083)")
abline(h=0,col="red",lwd=2)
m <- lapply(CTRIES, function(XXX,.e0difff){
    lines(e0difff[[XXX]], col = "#00000050")
  },.e0difff = e0difff)
lines(e0difff[["ISL"]], col = "royalblue")
lines(e0difff[["NLD"]], col = "green3")
lines(e0difff[["SWE"]], col = "orange")
legend("topright",col = c("royalblue","green3","orange"),lty=1,legend=c("Iceland","Netherlands","Sweden"))
dev.off()


a0compm <- parallel::mclapply(CTRIES, function(XXX){
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    v5       <- ltper_AxN(sex = "m", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
    v5a0vv6  <- ltper_AxN(sex = "m", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
    
    cbind(x=as.integer(sort(unique(v5$Year))),
      a0v5=with(v5,ax[Age == "0"]),
      a0v6=with(v5a0vv6,ax[Age == "0"]))
  }, mc.cores = 4)
names(a0compm) <- CTRIES
a0compf <- parallel::mclapply(CTRIES, function(XXX){
    WRKNG    <- paste0( "/data/wilmoth0/HMD/HMDWORK/",XXX)
    v5       <- ltper_AxN(sex = "f", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = FALSE, save.bin = FALSE)
    v5a0vv6  <- ltper_AxN(sex = "f", WORKING = WRKNG, 
      MPVERSION = 5, testa0v6 = TRUE, save.bin = FALSE)
    
    cbind(x=as.integer(sort(unique(v5$Year))),
      a0v5=with(v5,ax[Age == "0"]),
      a0v6=with(v5a0vv6,ax[Age == "0"]))
  }, mc.cores = 4)
names(a0compf) <- CTRIES


png("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/Swedena0v5v6.png")
plot(a0compm[["SWE"]][,1],a0compm[["SWE"]][,2],type = "l", col = "blue",lty=2, ylim = c(0,.35),
  main= "Sweden a0 v5 'Coale-Demeny' versus v6 'Andreev-Kingkade'", ylab = "a0",xlab = "Year")
lines(a0compm[["SWE"]][,1],a0compm[["SWE"]][,3], col = "blue")
lines(a0compf[["SWE"]][,1],a0compf[["SWE"]][,2], col = "red",lty=2)
lines(a0compf[["SWE"]][,1],a0compf[["SWE"]][,3], col = "red")
legend("bottomleft",col = c("red","red","blue","blue"),lty=c(1,2,1,2),
  legend=c("v6 females","v5 females","v6 males","v5 males"))
dev.off()
