#' @title \code{YearAgg} Aggregate an 'age by year' matrix into N-year groups
#' 
#' @description this function takes any 'age by year' matrix, ideally with columns and row appropriately labeled, and aggregates it into year groups. The first year of each group are years where \code{year %% N == 0}. if \code{N = 1}, nothing happens, for any \code{N > 1} the minimum number of years for a year group is 2. That is to say, for 5-year groupings, the first and last year-groups might not have divided evenly into 5. These are included if they have at least 2 years in them.
#' 
#' @param mat an 'age by year' numeric matrix, ideally with dimensions labeled.
#' @param years vector of years (columns) in \code{mat}. Defaults to the column names of \code{mat}
#' @param N the desired year-group width.
#' 
#' @return an 'age by year' matrix with columns. Column labls will now specify hyphen-separated year ranges.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export 
#' @importFrom compiler cmpfun

YearAgg <- cmpfun(function(mat, years = as.integer(colnames(mat)), N = 1){
    if (N == 1){                            # no need to aggregate single years
      return(mat)
    }
    fac       <- years - years %% N         # grouping factor
    fac[fac < years[1]] <- years[1]
    out <- t(apply(mat, 1, function(x, fac){    # sum by factor
                       tapply(x, fac, function(.x){ # modified to be BEL-friendly
                                                    # sum(!is.na(.x)) < 2 is for consistentcy
                                                    # with idea of min 2-year groups if N > 1
                                                    # all(is.na(.x)) is consistent with matlab...
                                 ifelse(sum(!is.na(.x)) < 2, NA, sum(.x, na.rm = TRUE))
                               }
                       )
                     }, fac = fac
                   )
             )
           
    yr0       <- as.integer(colnames(out))  # lable columns to indicate ranges
    yrlast    <- sort(unique(c(years[years %% N == (N - 1)], max(years))))
    colnames(out) <- paste(yr0, yrlast, sep = "-" )
    # this will remove 1 year groupings from 5 or 10 year aggregations
    out       <- out[, yr0 != yrlast]
    return(out)
  }  # end function definition
)


