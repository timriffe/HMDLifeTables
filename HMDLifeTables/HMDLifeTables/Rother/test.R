
devtools::load_all("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables")
#devtools::document("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables")
WORKING <- "/hdir/0/triffe/HMDWORK/USA"
XXX <- "USA"

#BuildRLifeTablePackage()
args(RunHMDCountry)
RunHMDCountry(WORKING = WORKING)
args(Exposures_per)

test <- ltper_AxN(WORKING, save.bin=FALSE)
v5m <- Exposures_per(WORKING = "/hdir/0/triffe/Desktop/FRATNP",MPVERSION=5)
v6m <- Exposures_per(WORKING = "/hdir/0/triffe/Desktop/FRATNP",MPVERSION=6)
v5f <- Exposures_per(WORKING = "/hdir/0/triffe/Desktop/FRATNP",MPVERSION=5,sex="f")
v6f <- Exposures_per(WORKING = "/hdir/0/triffe/Desktop/FRATNP",MPVERSION=6,sex="f")

library(LexisUtils)
#pdf("/home/tim/workspace/HMDJan2016/Figures/LexisAPCohortEffect.pdf")
LexisMap(v6[51:71, as.character(1968:1988)] / v5[50:70, as.character(1968:1988)]  ,log=FALSE )
#dev.off()

save(v5m,file="/hdir/0/triffe/workspace/other/DATA/FRATNPExposuresv5m.Rdata")
save(v6m,file="/hdir/0/triffe/workspace/other/DATA/FRATNPExposuresv6m.Rdata")
save(v5f,file="/hdir/0/triffe/workspace/other/DATA/FRATNPExposuresv5f.Rdata")
save(v6f,file="/hdir/0/triffe/workspace/other/DATA/FRATNPExposuresv6f.Rdata")
