 
# Author: triffe
###############################################################################

source("/data/commons/triffe/git/HMD_Rlifetables_git/R/RunHMDCountry.R")

cntries <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE","DNK","DEUTE","DEUTW","DEUTNP",
  "ESP","EST","FIN","FRATNP","GBR","GBRCENW","GBR_NIR","GBR_NP","GBR_SCO",
  "HUN","IRL","ISL","ITA","JPN","LTU","LVA","NLD","NOR","POL","PRT",
  "NZL_MA","NZL_NM","NZL_NP","SWE","USA","SVK","SVN","UKR","FRACNP","GBRTENW","ISR","LUX","RUS","TWN")

library(parallel)
(NC <- detectCores()) # how many cores ar available?
Sys.time()
null.list <- mclapply(cntries, RunHMDCountry, mc.cores = NC)
Sys.time()
1345831426-1345831118 #

# 6 or 7 minutes on Galton

# A <- system.time(source("/data/commons/triffe/git/HMD_Rlifetables_git/R/GaltonParallelScript.R"))

# now organize for balanced load
cntries <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE","DNK","DEUTE","DEUTW","DEUTNP",
  "ESP","EST","FIN","FRATNP","GBR","GBRCENW","GBR_NIR","GBR_NP","GBR_SCO",
  "HUN","IRL","ISL","ITA","JPN","LTU","LVA","NLD","NOR","POL","PRT",
  "NZL_MA","NZL_NM","NZL_NP","SWE","USA","SVK","SVN","UKR","FRACNP","GBRTENW","ISR","LUX","RUS","TWN")

# cohort lifetables in:
clt.list <- c("DNK","FIN","FRACNP","FRATNP","ISL","NOR","SWE","GBRCENW","GBRTENW","ITA","CHE","GBR_SCO","NLD")
short.list <- c("TWN", "SVN","ISR","CHL","DEUTNP","LUX", "LVA","LTU" )
med.list <- c("NZL_NM", "ESP","USA","CAN","AUS")
med2.list <- cntries[!cntries %in% clt.list & !cntries %in% short.list & !cntries %in% med.list]
# very roughly sorted by data size
cntries.b <- c(clt.list, med.list, med2.list, short.list)
cntries.b <- c("DNK", "FIN", "FRACNP", "FRATNP", "ISL", "NOR", "SWE", "GBRCENW", 
  "GBRTENW", "ITA", "CHE", "GBR_SCO", "NLD", "NZL_NM", "ESP", "USA", 
  "CAN", "AUS", "AUT", "BGR", "BLR", "CZE", "DEUTE", "DEUTW", "EST", 
  "GBR", "GBR_NIR", "GBR_NP", "HUN", "IRL", "JPN", "POL", "PRT", 
  "NZL_MA", "NZL_NP", "SVK", "UKR", "Sys.time(RUS", "TWN", "SVN", "ISR", 
  "CHL", "DEUTNP", "LUX", "LVA", "LTU")

Sys.time()
null.list <- mclapply(cntries.b, RunHMDCountry, mc.cores = NC)
Sys.time()
# brings back down to 4 minutes

# lesson: country order makes a big difference




