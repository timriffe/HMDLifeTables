#' @title \code{Abridge} a lifetable from single ages to 5-year abridged age groups
#'
#' @description
#' \code{Abridge} is an auxiliary function called by any of the four top level lifetable functions (\code{LTper_AxN()}, ...). This function can be called on its own if you dispose of three equally dimensioned matrices of \code{lx, Tx, ex} values (age by year) in single ages (years may or may not be grouped, it doesn't matter). \code{OPENAGE} is only used for age labelling of the output, which is a full-fledged lifetable as a \code{data.frame} in long (stacked) format, unrounded and unformatted. Procedure were roughly described on page 44 of MPv5.
#'
#' @details the \code{ax} method was not documented in MPv5, but was gleaned from the previous matlab version. Closeout procedures are basically the same as with single age lifetables. 
#'
#' @param lx a matrix (age by year) of lx values, as passed in by e.g. \code{LTper_AxN()}
#' @param Tx a matrix (age by year) of Tx values, as passed in by e.g. \code{LTper_AxN()}
#' @param ex a matrix (age by year) of ex values, as passed in by e.g. \code{LTper_AxN()}
#' @param OPENAGE the lifetable open age, as passed in by e.g. \code{LTper_AxN()} where the default value is 110
#' 
#' @return a \code{data.frame} of the full lifetable in long (stacked) format, unrounded and unformatted (except the \code{Age} column, which is centered on the age-group hyphens).
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
# Author: triffe
###############################################################################
Abridge <- function(lx, Tx, ex, OPENAGE){
  abr.ages <- c(0, 1, seq(5, OPENAGE, by = 5))
  n   <- diff(abr.ages)
  n   <- matrix(n, ncol = ncol(lx), nrow = nrow(lx))
  
  lx  <- lx[as.integer(rownames(lx)) %in% abr.ages, ]
  Tx  <- Tx[as.integer(rownames(Tx)) %in% abr.ages, ]
  ex  <- ex[as.integer(rownames(ex)) %in% abr.ages, ]
 
  # px, dx, qx, Lx, mx all 1 row shorter, with special close-out
  px  <- lx[2:nrow(lx), ] / lx[1:(nrow(lx) - 1), ]
  px[is.infinite(px)] <- NA
  dx  <- lx[1:(nrow(lx) - 1), ] - lx[2:nrow(lx), ]
  qx  <- 1 - px
  Lx  <- Tx[1:(nrow(Tx) - 1), ] - Tx[2:nrow(Tx), ]
  mx  <- dx / Lx
  mx[is.infinite(mx)] <- NA
  ax  <- (Lx - n * lx[2:nrow(lx), ]) / dx
  # TR: 5-Aug-2020. Hypothetically could get 0 deaths in age 1-4.
  ind <-dx == 0 & !is.na(qx)
  ax[ind] <- n[ind] / 2
 
  # prepare open age group entries:
  dxOP <- lx[nrow(lx), ]
  LxOP <- Tx[nrow(Tx), ]
  mxOP <- dxOP / LxOP
  mxOP[is.infinite(mxOP) | is.nan(mxOP)] <- NA
  qxOP <- ifelse(qx[nrow(qx), ] != 1 & !is.na(qx[nrow(qx), ]), 1, NA)
  axOP <- ex[nrow(ex), ]
  
  # stick together:
  dx <- rbind(dx, dxOP)
  Lx <- rbind(Lx, LxOP)
  mx <- rbind(mx, mxOP)
  qx <- rbind(qx, qxOP)
  ax <- rbind(ax, axOP)
 
  # prepare output:
  # ages formatted for txt printing later
  ages <- c("    0  ", 
            paste(sprintf("%3s", c(1, seq(5, OPENAGE - 5, by = 5))),
                  sprintf("%-3s", c(4, seq(9, OPENAGE - 1, by = 5))), 
                  sep = "-"),
            paste0(OPENAGE,"+   ")
            )
  yrs <- colnames(Tx)
  # stick together
  output <- data.frame(Year = rep(yrs, each = length(ages)),
    Age = rep(ages, length(unique(yrs))),
    mx = as.vector(mx),
    qx = as.vector(qx),
    ax = as.vector(ax),
    lx = as.vector(lx),
    dx = as.vector(dx),
    Lx = as.vector(Lx),
    Tx = as.vector(Tx),
    ex = as.vector(ex),
    stringsAsFactors = FALSE
  )
} # end function
