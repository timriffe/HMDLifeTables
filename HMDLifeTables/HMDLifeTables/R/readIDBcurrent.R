#'
#' @title \code{readIDBcurrent()} read in any IDB file to a data.frame
#' 
#' @description Character columns are read in as character. Missing values are assumed to be indicated with \code{"."}, due to convention. We assume there are column headers. File must be comma separated. All of things will hold for standard InputDB files.
#' 
#' @param XXX character string; HMD country short code
#' @param item character string; one of \code{c("death","pop","birth","tadj","monthly")}
#' @param old logical should the file come from XXX_Old?
#' @param full.path optional. *If you're not at UCB, or if you're no working from the main HMD working directory \code{"~/HMDWORK"}, then give the full path to the file, including \code{".txt"}.
#' 
#' @export
#' 
#' @return data.frame for the given input file.
#' 
#' @examples
#' #bo <- readIDBcurrent("SWE","birth")
#' #bo <- resortBirths(bo)
#' 
readIDBcurrent <- function(XXX , item , old = FALSE, full.path = NULL){
  
  stopifnot(item %in% c(c("death","pop","birth","tadj","birthbymonth")))
  
  this.path <- ifelse(is.null(full.path), 
                      ifelse(old,
                        file.path("~/HMDWORK", paste0(XXX,"_Old"), "InputDB", paste0(XXX, item, ".txt")), 
                        file.path("~/HMDWORK", XXX, "InputDB", paste0(XXX, item, ".txt"))), 
                      full.path)
  ## CAB: separate call for tadj file will go here (developed in Rother/ for now)
  ## CAB: should have general tests for duplicate items once LDB is taken into account
  ## CAB: replace below with item-specific reads and implement tests for column headers - we have no checks
  ## CAB: that the correct columns are being read in and labeled correctly
  obj <- read.table(this.path, 
             sep = ",", 
             na.strings = ".", 
             header = TRUE, 
             stringsAsFactors = FALSE,
             strip.white = TRUE)
  invisible(obj)
}
