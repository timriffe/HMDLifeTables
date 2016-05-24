
# Author: triffe
###############################################################################
#?install.packages
#install.packages("/home/tim/git/HMDLifeTables/HMDLifeTables/HMDLifeTables_6.2.4752.tar.gz",
#type="source",repos=NULL)
CallGraph <- function(LTfunction = "ltper_AxN"){
  mvbutils::foodweb(where = "package:HMDLifeTables", prune = LTfunction, generics = c("ExtractXXXfromWORKING"),
    boxcolor = gray(.95))
  title(paste0(LTfunction, "() call tree"))
}
rownames(installed.packages())
#devtools::install_github("mvbutils")
# these may need resizing to be legible\library(HMDLifeTables)
CallGraph("ltcoh_AxN")
CallGraph("ltcohBoth_AxN")
CallGraph("ltper_AxN")
CallGraph("ltperBoth_AxN")
CallGraph("RunHMDCountry")
CallGraph("perTadj")


