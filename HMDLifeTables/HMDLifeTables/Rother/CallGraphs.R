
# Author: triffe
###############################################################################
#install.packages("mvbutils")
CallGraph <- function(LTfunction = "ltper_AxN"){
  mvbutils::foodweb(where = "package:RLifeTable", prune = LTfunction, generics = c("ExtractXXXfromWORKING"),
    boxcolor = gray(.95))
  title(paste0(LTfunction, "() call tree"))
}

# these may need resizing to be legible
CallGraph("ltcoh_AxN")
CallGraph("ltcohBoth_AxN")
CallGraph("ltper_AxN")
CallGraph("ltperBoth_AxN")
CallGraph("RunHMDCountry")