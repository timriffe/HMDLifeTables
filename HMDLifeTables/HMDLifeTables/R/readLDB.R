#'
#' @title function to read Lexis database files into R
#' 
#' @description This function reads in a given LexisDB file, imputing NAs properly (but not removing them), converting the Lexis column to character, and appending a Sex column, inferred from the first letter of the file name, per HMD conventions.
#' 
#' @param ldb.path path to the standard Lexis database file, including file extension \code{.txt}. File name should be of the form \code{"mESP.txt"} or \code{"fESP"}, where ESP could be any country code.
#' 
#' @return data.frame LexisDB object, with columns in a useable form and a column for Sex added
#' 
#' @export
#' 

readLDB <- function(ldb.path){
  LDBobj        <- read.table(ldb.path, header = FALSE, sep = ",", 
    col.names = c("Year", "Age", "Lexis", "Cohort", "Population", "Deaths"))
  LDBobj[LDBobj == -1] <- NA
  LDBobj$Lexis  <- ifelse(LDBobj$Lexis == 1, "TL","TU")
  sex           <- unlist(strsplit(rev(unlist(strsplit(ldb.path,split="/")))[1],split=""))[1]
  LDBobj$Sex    <- sex
  LDBobj
}



