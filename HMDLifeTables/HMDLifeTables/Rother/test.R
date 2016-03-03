
devtools::load_all("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables")
#devtools::document("/data/commons/triffe/git/HMDLifeTables/HMDLifeTables/HMDLifeTables")
WORKING <- "/hdir/0/triffe/HMDWORK/USA"
XXX <- "USA"

#BuildRLifeTablePackage()
args(RunHMDCountry)
RunHMDCountry(WORKING = WORKING)

test <- ltper_AxN(WORKING)
test[test$Year == 2009, ]


library(devtools)
library(TimUtils)
IncrementVersion("/home/tim/git/HMDLifeTables/HMDLifeTables/HMDLifeTables", 
		major.version = "6", package.origin = "2013-01-01")
build("/home/tim/git/HMDLifeTables/HMDLifeTables/HMDLifeTables", "/home/tim/git/HMDLifeTables/HMDLifeTables")
?build

#######################################
#
