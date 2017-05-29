#' @title \code{ltper_mx_v5} fits Kannisto mortality to old ages according to version 5 of the methods protocol
#' 
#' @description This function optimizes to Kannisto likelihood to ages 80+ in the data, and it returns the entire 'age by year' \code{mx} matrix, properly imputed according to \code{extrap.ages.i} (the output from \code{ltper_getDx100()}). This function is internal, and called only by \code{ltper_AxN()}.
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


ltper_mx_v5 <- function(Dx, Exp, extrap.ages.i){
  
  # ---------------------------------------------------
  # define internal functions
  # ---------------------------------------------------
  # Kannisto MLE fitting functions as of MP version 5
  # Eq 53, pg 36 MPv5
  # a,b, are unconstrained in these functions, since we can set bounds with L-BFGS-B method.
  KannistoMu <- function(pars, x = .5:30.5){
      a        <- pars["a"]
      b        <- pars["b"]
      (a * exp(b * x)) / (1 + a * exp(b * x))
  }
 
  # likelihood, parameters constrained to positive using abs() only
  # This is actually partial log likelihood, which matters when getting the information
  # from the Hessian
  KannistoLik <- function(pars, .Dx, .Exp, .x. = .5:30.5){
      mu       <- KannistoMu(pars, x = .x.)
      # take negative and minimize it (default optimizer behavior)
      -sum(.Dx * log(mu) - .Exp * mu, na.rm = TRUE)
  }

  # vector of partials with respect to a and b, respectively
  KannistoGr <- function(pars, .Dx, .Exp, .x. = .5:30.5){
      a   <- pars["a"]
      b   <- pars["b"]
      d.a <- (a * exp(b * .x.) * .Exp + (-a * exp(b * .x.) - 1) * .Dx) /
              (a ^ 3 * exp(2 * b * .x.) + 2 * a ^ 2 * exp(b * .x.) + a)
      d.b <- (a * .x. * exp(b * .x.) * .Exp + (-a * .x. * exp(b * .x.) - .x.) * .Dx) /
              (a ^ 2 * exp(2 * b * .x.) + 2 * a * exp(b * .x.) + 1)
      colSums(cbind(a = d.a, b = d.b), na.rm = TRUE) 
  }

  # starting value finder, since sometimes the L-BFGS-L can fail to converge if it goes off
  # in the wrong direction. unfortunately. 
  startingValueGridSearch <- function(
      a.vals = seq(.001, .701, by = .1),
      b.vals = seq(.001, .701, by = .1), 
      .Exp, .Dx, .x. = .5:30.5, i){
    #a.vals = seq(.001, .701, by = .1)
    #b.vals = seq(.001, .701, by = .1)   
    # vectors of length = len(a)xlen(b) -- vectorized grid
    ab.grid <- expand.grid(a.vals, b.vals)
    a.vec <- ab.grid[, 1]
    b.vec <- ab.grid[, 2]
      
    valsLik <- sapply(seq(along=a.vec), 
                      function(.i){ KannistoLik(c(a = a.vec[.i], b = b.vec[.i]), 
                                  .Dx = .Dx, .Exp = .Exp, .x. = .x.)} )
      # MP.likM <- matrix(ncol = length(a.vals), nrow = length(b.vals))
      # for (a. in 1:length(a.vals)){
      #   for (b. in 1:length(b.vals)){
      #     MP.likM[b., a.] <- KannistoLik(c(a = a.vals[a.], b = b.vals[b.]), 
      #       .Dx = .Dx, .Exp = .Exp, .x. = .x.)                     
      #   } 
      # }
    arg.min <- which.min(valsLik)
      # c(a = a.vals[col(MP.likM)[which.min(MP.likM)]], 
      #   b = b.vals[row(MP.likM)[which.min(MP.likM)]])
    
    return( c(a = a.vec[arg.min], b = b.vec[arg.min]) )
  }
  
  # -------------------------------------------------------------
  # begin actual fitting 
  ages    <- as.integer(rownames(Dx))
  ages.i  <- ages[81:length(ages)] + .5 - 80
  years   <- colnames(Dx)
  par.est <- matrix(nrow = length(years), ncol = 2, dimnames = list(years, c("a","b")))
  se.est <- matrix(nrow = length(years), ncol = 2, dimnames = list(years, c("a","b")))
  
  # in a for-loop because we're iterating over 2 items..Dx and Exp
  by.j    <- c(.1, .05, .02, .005, .001, .0002, .00005, .00001) # last entry, nr 8, is moot. 
  # I think all cases converge .001,
 
    for (i in 1:ncol(Dx)){ # iterate over years
      if(!all(is.na(Dx[,i]))){  # added for BEL
        # the first set of starting values typically suffices. It's over a wide range, but quite rough
        afromto.j <- c(.001, .701)
        bfromto.j <- c(.001, .701)
        # intitial starting values
        starts.j  <- startingValueGridSearch(a.vals = seq(afromto.j[1], afromto.j[2], by = by.j[1]),
          b.vals = seq(bfromto.j[1], bfromto.j[2], by = by.j[1]), 
          .Dx = Dx[81:nrow(Dx), i], 
          .Exp = Exp[81:nrow(Exp), i],
          .x. = ages.i )
        
        # jump into optimization loop- keep trying better starts until it converges
        
        for (j in 1:7){ # j <- 1
          # try optimizing with present starts
          result <- try(optim(
              starts.j, 
              fn = function(...) KannistoLik(...) , # likelihood function rescaled in the optim call
              gr =  function(...) KannistoGr(...) , # gradient function rescaled in the optim call
              method = c("L-BFGS-B"),                     # BFGS with par constraints
              upper = c(5, 5), lower = c(0, 0),           # upper, lower bounds for a,b
              control = list(lmm = 20, fnscale=1.0e6, trace=0), hessian=TRUE,   # precision
              .Dx = Dx[81:nrow(Dx), i], .Exp = Exp[81:nrow(Exp), i], .x. = ages.i # args to pass in
            ), 
            silent = FALSE
          );
          
          # if it didn't converge, an error will be thrown of class 'try-error'
          if (class(result) == "try-error" || result$convergence > 0){
            ##cat(result$message, "\n\n")  # what went wrong?
            ## try again, but with a different fnscale, 
            result <- try(optim(
              starts.j, 
              fn = function(...) KannistoLik(...) , # likelihood function rescaled in the optim call
              gr =  function(...) KannistoGr(...) , # gradient function rescaled in the optim call
              method = c("L-BFGS-B"),                     # BFGS with par constraints
              upper = c(5, 5), lower = c(0, 0),           # upper, lower bounds for a,b
              control = list(lmm = 20, fnscale=1.0e2, trace=2), hessian=TRUE,   # precision
              .Dx = Dx[81:nrow(Dx), i], .Exp = Exp[81:nrow(Exp), i], .x. = ages.i # args to pass in
            ), 
            silent = FALSE
            );
          }
            
          
          # if it still didn't  converge, refine the grid search and retry
          if (class(result) == "try-error" || result$convergence > 0){
            ##cat(result$message, "\n\n")  # what went wrong?
            # in this case we select from a finer grid of starting values 
            # within narrower limits and try again ('smart' grid search)
            afromto.j <- c(max(starts.j["a"] - 2 * by.j[j], 0), starts.j["a"] + 2 * by.j[j])
            bfromto.j <- c(max(starts.j["b"] - 2 * by.j[j], 0), starts.j["b"] + 2 * by.j[j])
            starts.j  <- try(startingValueGridSearch(a.vals = seq(afromto.j[1], afromto.j[2], by = by.j[j + 1]),
                b.vals = seq(bfromto.j[1], bfromto.j[2], by = by.j[j + 1]), 
                .Dx = Dx[81:nrow(Dx), i], 
                .Exp = Exp[81:nrow(Exp), i], .x. = ages.i))
            next
          } else {
            ## return parameter estimates and their SE
            ab.j <- result$par
            # stats.stackexchange.com:   "in R given an output from optim  hessian how to calculate parameter confidence intervals"
            fisher_info <- solve( result$hessian )
            se.j  <-  sqrt(diag(fisher_info))
            break
          }       
        } # j
       
        # it only gets this far in extreme testing, allow for sake of simulations
        if (class(result) == "try-error"){
          ab.j <- starts.j
          se.j <- c(NA,NA)
        }
        
      } else {
        ab.j <- c(NA, NA) +1  # class should be numeric
        se.j <- c(NA, NA)
      }
      
      par.est[i, ] <- ab.j    
      se.est[i,]   <- se.j
      cat(paste(i, ab.j[1], ab.j[2], "\n"))
    } #i
  
  # fitted values over entire age range
  mx.est  <- apply(par.est, 1, function(parsest){
      KannistoMu(parsest, x = c(0:(nrow(Dx)-1) - 79.5))
    }
  )
  
  # additional redundant check for parameters:
  if (class(par.est[,1]) != "numeric" | class(par.est[,2]) != "numeric"){
    stop("looks like one of the parameters didn't converge properly, check: ", sex, ctry)
  }

  # raw Mx
  mx <- Dx / Exp
  # impute according to 'extrap.ages.i', which is supplied as such
  for (i in seq_along(years)){
    if (!all(is.na(mx[,i]))){ # added for BEL
      ind1        <- extrap.ages.i[i]:nrow(mx)
      mx[ind1, i] <- mx.est[ind1, i]
    }
  }
  # return the finished mx
  #list(mx=mx, par=par.est,  se=se.est)
  return(mx)
}


