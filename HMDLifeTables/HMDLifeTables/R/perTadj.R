#' @title \code{perTadj} a function for doing lifetable territorial adjustments for period data
#' 
#' @description period territorial adjustments 'un'adjust population counts bring January 1 population estimates from year t + 1 to December 31  of year t, since both populations are needed to calculate year t exposures. This step should in the future be moved back to the LDB programs, which should return a more comprehensive Lexis Database object. This function is called only by \code{getPeriodComponents()}.
#' 
#' @param LDBobj Lexis database object, with proper column names, to be passed in by \code{getPeriodComponents()}.
#' @param tadj a territorial adjustment object, as read straight from the Input Database.
#' @param sex \code{"m"} or \code{"f"}
#' 
#' @return returns LDBobj back again, having 'un'-adjusted deaths and pop counts for pertinent years.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export
#' 
#' @importFrom compiler cmpfun

perTadj <- cmpfun(function(LDBobj, tadj, sex){
  # some tadj files also contain "Rb", limit to "Vx" type
  FY          <- min(LDBobj$Year) # bug fix for ITA (tadj in 1st year makes no sense)
  tadj        <- tadj[tadj$Sex == sex & tadj$Type == "Vx" & tadj$Year > FY,]
  # make sure properly ordered (LDBobj sorted prior to calling this)
  tadj        <- tadj[order(tadj$Year, tadj$Age),]
  # this selects the LDBobj indices that are present in tadj
  tadj.i      <- LDBobj$Age %in% tadj$Age & LDBobj$Year %in% tadj$Year & LDBobj$Lexis == 2
  # add tadj column, by default does nothing
  LDBobj$tadj <- 1
  # ITA has tadj years that aren't in the data... used in LDB?
  if (sum(tadj.i) != length(tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year])){
    stop("tadj prob, what country is this?")
  }
  # the denominator selects only the tadj items present in LDBobj (match to tadj.i)
  # matching is guaranteed because tadj and LDBobj have been sorted using the same criteria
  # we take the inverse because we're UNDOING the tadj
  LDBobj$tadj[tadj.i] <- 1 / tadj$Value[tadj$Age %in% LDBobj$Age & tadj$Year %in% LDBobj$Year]
  # now, simply mutliply into Population
  LDBobj$Population   <- LDBobj$Population * LDBobj$tadj
  # this later step is simply useful for indexing outside this function, which relies on position of NAs
  LDBobj$Population[is.na(LDBobj$Population)] <- 0 # necessary so we can use NAs to index later
  LDBobj
})
