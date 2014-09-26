
# Author: triffe
###############################################################################

# dynamically reload package, without prior build
devtools::load_all("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables", TRUE)

setwd("/hdir/0/hmd/HMDWORK/ITA")
RunHMDCountry()
args(RunHMDCountry)


# Run on a Country
setwd("/hdir/0/hmd/HMDWORK/BEL")
RunHMDCountry()
getwd()
# compare against current installed version
RLifeTable::RunHMDCountry()

# in case in-line testing is necessary:
library(reshape2)
library(compiler)

# this builds the doc files, if they have been modified 
devtools::document("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables")

# custom build functions
source("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/BuildUtils.R")
BuildRLifeTablePackage()     # definitive rebuild 

NewestRLifeTablePackage()    # newest built, but not necessarily installed, build

InstalledRLifeTablePackage() # present installed version

#library(RLifeTable) # updated less often
RunHMDCountry(WORKING = "/hdir/0/hmd/HMDWORK/ITA")
# testthat perhaps not necessary:
# test output in STATS vs RSTATS

source("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/Rother/TestFunctions.R")
test_Contents("/hdir/0/hmd/HMDWORK/BEL")
test_ltper(WORKING = "/hdir/0/hmd/HMDWORK/SWE")
test_ltcoh(WORKING = "/hdir/0/hmd/HMDWORK/ITA")
test_file_nr_diff(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "ltcoh")
test_file_nr_diff(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "Deaths_1")
test_file_nr_diff(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "Deaths_5")
test_file_nr_diff(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "Deaths_l")
test_file_nr_diff(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "Population.txt")

# causes hangup:
# (fix = allow for dimension mismatch by filling tadj object out to 130 years..)
# A <- Population_A(WORKING = "/hdir/0/hmd/HMDWORK/ITA", OPENAGE = 111, save.bin = FALSE)
# works, cuz tadj goes to age 110.
#A <- Population_A(WORKING = "/hdir/0/hmd/HMDWORK/SWE", OPENAGE = 111, save.bin = FALSE)
#tail(A)
#A[A$Year == 2011, ]
#

list_files("/hdir/0/hmd/HMDWORK/SWE")
#tkdiff_file(WORKING = "/hdir/0/hmd/HMDWORK/ITA", "Population.txt")


#----------------------------------------------------------------------------|
# TODO: ltperBoth_AxN() has been matched to matlab, there are likely 2       |
#                       bugs that this entails, for which correct            |
#                       code has been commented out. Awaiting review         |
#                       verified for ITA and BEL, couldn't update SWE files. |                         
#----------------------------------------------------------------------------|
# TODO: YearAgg() 2 year minimum for N-year age groups not enforced for BEL  |
#                 intermediate missing years in matlab code. It IS enforced  |
#                 in current R implementation, meaning that 1915-1919 5-year |
#                 cohort is an NA in R, but is based solely on 1919 data for |
#                 matlab. Just bear this in mind when calibrating. Perhaps   |
#                 not worth                                                  |
#----------------------------------------------------------------------------|
# TODO: matlab bug; age 110+ Population.txt doesn't always aggregate         |
#----------------------------------------------------------------------------|
# TODO: matlab bug; both sex cohort lifetables don't weight ax properly      |
#----------------------------------------------------------------------------|


