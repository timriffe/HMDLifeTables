#' @title \code{ltper_mx_v7} fits Kannisto mortality to old ages according to version 7 of the methods protocol
#' 
#' @description This function optimizes to Kannisto likelihood to ages 70+ in the data, and it returns the entire 'age by year' \code{mx} matrix, properly imputed according to \code{extrap.ages.i} (the output from \code{ltper_getDx100()}). This function is internal, and called only by \code{ltper_AxN()}. It differs from \code{ltper_mx_v5()} in that ages 70-79 are included, but their information is exponentially weighted according to the 80+ population size: smaller populations use more information from the 70-79 range, while large populations asymptotically approach version 5 smoothing.
#' 
#' @param Dx 'age by year' numeric matrix of death counts, as passed from \code{ltper_AxN()}.
#' @param Exp 'age by year' numeric matrix of exposures, as passed from \code{ltper_AxN()}.
#' @param extrap.ages.i vector of age indices where smoothed mx values should start to be imputed.
#' 
#' @return an 'age by year' numeric matrix of mx values, to be used for all subsequent lifetable calculations.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#'
#' @export
#' 
#' @importFrom compiler cmpfun
# Author: triffe
###############################################################################
ltper_mx_v7 <- function(Dx, Exp, extrap.ages.i){
  
  # ---------------------------------------------------
  # define internal functions
  # ---------------------------------------------------
  # Kannisto MLE fitting functions as of MP version 5
  # Eq 53, pg 36 MPv5
  # a,b, are unconstrained in these functions, since we can set bounds with L-BFGS-B method.
  KannistoMu <- cmpfun(function(pars, x = c(-9.5:30.5)){
      a        <- pars["a"]
      b        <- pars["b"]
      (a * exp(b * x)) / (1 + a * exp(b * x))
    })
  # Kannisto MLE functions for weighted version, possibly for v6
  KannistoLik_w <- cmpfun(function(pars, .Dx, .Exp, .x. = c(-9.5:30.5)){
      mu       <- KannistoMu(pars, x = .x.)
      # take negative and minimize it (default optimizer behavior)
      Exp80p   <- sum(.Exp[11:41])
      wi       <- rep(1, length(.Exp))
      wi[1:10] <- exp(-c(9.5:.5) * Exp80p / 10000)
      # insert weights
      -sum(wi * ( .Dx * log(mu) - .Exp * mu ), na.rm = TRUE) 
    })
# Gradient function.
  KannistoGr_w <- cmpfun(function(pars, .Dx, .Exp, .x. = c(-9.5:30.5)){
      a        <- pars["a"]
      b        <- pars["b"]
      # make weights vector
      Exp80p   <- sum(.Exp[11:41])
      wi       <- rep(1, length(.Exp))
      wi[1:10] <- exp(-c(9.5:.5) * Exp80p / 10000)
      # simply multiply in the weights
      d.a      <- wi * (a * exp(b * .x.) * .Exp + ( -a * exp(b * .x.) - 1) * .Dx) /
        (a ^ 3 * exp(2 * b * .x.) + 2 * a ^ 2 * exp(b * .x.) + a)
      d.b      <- wi * (a * .x. * exp(b * .x.) * .Exp + ( -a * .x. * exp(b * .x.) - .x.) * .Dx) /
        (a ^ 2 * exp(2 * b * .x.) + 2 * a * exp(b * .x.) + 1)
      colSums(cbind(a = d.a, b = d.b), na.rm = TRUE) 
    })
  
  # starting value finder, since sometimes the L-BFGS-L can fail to converge if it goes off
  # in the wrong direction. unfortunately. Compiled due to double loop.
  startingValueGridSearch_w <- cmpfun(function(
      a.vals = seq(.001, .701, by = .1),
      b.vals = seq(.001, .701, by = .1), 
      .Exp, .Dx, i){
      
      MP.likM <- matrix(ncol = length(a.vals), nrow = length(b.vals))
      for (a. in 1:length(a.vals)){
        for (b. in 1:length(b.vals)){
          MP.likM[b., a.] <- KannistoLik_w(c(a = a.vals[a.], b = b.vals[b.]), 
            .Dx = .Dx, .Exp = .Exp, .x. = c(-9.5:30.5))                     
        } 
      }
      c(a = a.vals[col(MP.likM)[which.min(MP.likM)]], 
        b = b.vals[row(MP.likM)[which.min(MP.likM)]])
    })
  # -------------------------------------------------------------
  # begin actual fitting 
  
  years   <- colnames(Dx)
  par.est <- matrix(nrow = length(years), ncol = 2, dimnames = list(years, c("a","b")))
  
  # in a for-loop because we're iterating over 2 items..Dx and Exp
  by.j    <- c(.1, .05, .02, .005, .001, .0002, .00005, NA) # last entry, nr 8, is moot. 
  # I think all cases converge .001,
  # at most 5 or so need .0002 grids
  # i <-1
  for (i in 1:ncol(Dx)){
    if (!all(is.na(Dx[,i]))){
      # the first set of starting values typically suffices. It's over a wide range, but quite rough
      afromto.j <- c(.001, .701)
      bfromto.j <- c(.001, .701)
      starts.j  <- startingValueGridSearch_w(a.vals = seq(afromto.j[1], afromto.j[2], by = by.j[1]),
        b.vals = seq(bfromto.j[1], bfromto.j[2], by = by.j[1]), 
        .Dx = Dx[71:111, i], 
        .Exp = Exp[71:111, i])
      # jump into optimization loop- keep trying better starts until it converges j <- 1
      for (j in 1:7){
        # try optimizing with present starts
        ab.j <- try(optim(
            starts.j, 
            fn = function(...) KannistoLik_w(...) * 1e-6, # likelihood function rescaled in the optim call
            gr =  function(...) KannistoGr_w(...) * 1e-6, # gradient function rescaled in the optim call
            method = c("L-BFGS-B"),                     # BFGS with par constraints
            upper = c(5, 5), lower = c(0, 0),           # upper, lower bounds for a,b
            control = list(factr = 1e-10, lmm = 20),    # precision
            .Dx = Dx[71:111, i], .Exp = Exp[71:111, i], .x. = c(-9.5:30.5) # args to pass in
          )$par, 
          silent = TRUE
        ) 
        # if it didn't converge, and error will be thrown of class 'try-error'
        if (class(ab.j) != "numeric"){
          # in this case we select from a finer grid of starting values 
          # within narrower limits and try again
          afromto.j <- c(max(starts.j["a"] - 2 * by.j[j], 0), starts.j["a"] + 2 * by.j[j])
          bfromto.j <- c(max(starts.j["b"] - 2 * by.j[j], 0), starts.j["b"] + 2 * by.j[j])
          starts.j  <- startingValueGridSearch_w(a.vals = seq(afromto.j[1], afromto.j[2], by = by.j[j + 1]),
            b.vals = seq(bfromto.j[1], bfromto.j[2], by = by.j[j + 1]), 
            .Dx = Dx[71:111, i], 
            .Exp = Exp[71:111, i]) 
          next
        } else {
          break
        }       
      }
    } else {
      ab.j   <- c(NA, NA)
    }
    par.est[i, ] <- ab.j            
  }
  
  # fitted values over entire age range
  mx.est  <- apply(par.est, 1, function(parsest){
      KannistoMu(parsest, x = (.5:110.5 - 80))
    }
  )
  # additional redundant check for parameters:
  if (class(par.est[, 1]) != "numeric" | class(par.est[, 2]) != "numeric"){
    stop("looks like one of the parameters didn't converge properly, check: ", sex, ctry)
  }
  #  if (write.log){
#    Totb0 <- sum(zapsmall(par.est[,2]) == 0)
#    cat(paste0("* optimal b = 0 in ",Totb0," cases"),
#      ifelse(Totb0 > 0, 
#        paste0("* smoothed mx will be flat in the following years: ",
#          paste(years[zapsmall(par.est[,"b"]) == 0],collapse = ", ")),""),
#      file = log.file, append = TRUE, sep = "\n")
#  }
  # raw Mx ncol(
  mx <- Dx / Exp
  # impute according to 'extrap.ages.i', which is supplied as such
  for (i in 1:length(years)){
    if (!all(is.na(mx[,i]))){
      ind1        <- extrap.ages.i[i]:nrow(mx)
      ind2        <- (extrap.ages.i[i] - 5):extrap.ages.i[i]
      w           <- c(1, .8, .6, .4, .2, 0)
      mx[ind1, i] <- mx.est[ind1, i]
      mx[ind2, i] <- w * mx[ind2, i] + (1 - w) * mx.est[ind2, i]
    }
  }
  # return the finished mx
  mx
}
