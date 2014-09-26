
# library(RLifeTables) # will work on Quigley.

# loads development version (this is what we do)
library(devtools)
load_all("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables")

# Requires this folder structure:

#XYZ/LexisDB
#   /InputDB
#   /STATS

#   /STATS
#   /RSTATS (optional, in case you want to compare matlab output with R output)
#   /Rbin (created by RunHMDCountry, contains unrounded R data.frames of each output file)

#   /CHECKS is where diagnostic plots live. Right now these are produced by SAS, but this
#           functionality will also move to R.

# for R CHECKS, see "/data/commons/triffe/git/HMD_Rlifetables_git/Diagnostics"
# for the current state of progress on that. The LexisDB scripts are the highest priority
# at this time, and we'll finish the diagnostics in R once that's done. 

# There is one main runner function, which assumes the above folder structure 
# (only needs the first two items in order to do its thing, and it will create a STATS
#  and Rbin folder)
# see:
?RunHMDCountry

# if you look inside it, all other functions in the package are documented in the same way.
# of course, any specific questions, just ask Tim or Carl.

# it's run like this:
RunHMDCountry(WORKING = "/data/commons/triffe/git/HMD_CS/HMDwork/C_BEL/BEL", 
  OPENAGE = 110,            # default
  RADIX = 1e5,              # default
  CAGEEXTRP = 90,           # default, also not relevant now for LAHMD
  MPVERSION = 5,            # default, can also be 6 (switching this year)
  STATSFOLDER = "RSTATS",   # could also just be "STATS" for primary use
  XXX = NULL,               # default, country code, inferred from WORKING otherwise
  LDBPATH = NULL,           # default, inferred from WORKING, but could be different
  IDBPATH = NULL)           # default, inferred from WORKING, but could be different
# running this produces all output as text files. If you want to distribute results in 
# spreadsheets, then it seems that translation can be done in a single pass, so an extra
# function can be written for that. If you want to write one, it can be worked into 
# RunHMDCountry() and toggled with an optional argument, such as:
# , makeSpreadsheets = FALSE, 

# It is possible to run most functions called by RunHMDCountry()
# without reference to the standard folder structure, feeding the function
# data.frames. That's good if you want to do live testing, however it's not
# necessary for routine operation.



