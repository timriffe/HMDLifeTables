#' @title \code{AKm02a0} estimates a0 using the Andreev-Kingkade rule of thumb.
#'
#' @description \code{AKm02a0} is an auxiliary function used by version 6 of the four HMD lifetable functions, \code{ltper_AxN()}, \code{ltcoh_AxN()}, \code{ltperBoth_AxN()}, \code{ltcohBoth_AxN()}. This is an approximation provided in an appendix by the authors. \code{AKm02a0_direct()} provides an analytical solution, which is overly precise and more laborious to explain in the MP. We prefer a separate table for m0 because it is only an approximation.
#'
#' @param m0 a value or vector of values of m0, the death probability for age 0 infants.
#' @param sex either "m" or "f"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export


AKm02a0 <- function(m0, sex = "m"){
  sex <- rep(sex, length(q0))
  ifelse(sex == "m", 
    ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
      ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
    # f
    ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
      ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}

#' @title \code{AKm02a0_direct} estimates a0 using an analytic translation of the Andreev-Kinkade rule of thumb.
#'
#' @description
#' \code{AKm02a0_direct} is an auxiliary function that could be used by version 6 of the Methods Protocol. Instead we opt for a simple approximate solution. This function calls \code{AKm02q0()} (also deprecated) to help get the work done, since the HMD needed to adapt the Andreev-Kingkade formulas to work with the period lifetable flow.
#'
#' @param m0 a value or vector of values of m0, the death rate for age 0 infants.
#' @param sex either "m" or "f"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export

AKm02a0_direct <- function(m0,sex="m"){
  sex <- rep(sex,length(m0))
  ifelse(sex == "m",
    ifelse(m0 < 0.02306737, 0.1493 - 2.0367 * AKm02q0(m0, 0.1493, -2.0367),
      ifelse(m0 < 0.0830706, 0.0244 + 3.4994 * AKm02q0(m0, 0.0244, 3.4994), .2991)),
    ifelse(m0 < 0.01725977, 0.1490 - 2.0867 * AKm02q0(m0, 0.1490, -2.0867),
      ifelse(m0 < 0.06919348, 0.0438 + 4.1075 * AKm02q0(m0, 0.0438, 4.1075), 0.3141))
  )
}




# Initial testing:
# backstop against automatic build
#runthis <- FALSE
#if (runthis){
## formula based on MPIDR Working Paper WP 2011-016
## Average Age at death in infancy: reconsidering
## the Coale-Demeny formulas at current levels of low mortality
#
## the formula given in that paper was revised, and the revision
## is at this writing not yet public. New forumula sent in email
## April 13th, 2013 from Dmitri Jadnov (Dima)
## -------------------------------------------------------
## new table (also in q0)
## Males
## lower q0     |   upper q0   | formula
## 0            | 0.0226       | 0.1493 - 2.0367 q0
## 0.0226       | 0.0785       | 0.0244 + 3.4994 q0
## 0.0785       |  +           | .2991
## ________________________________________
## Females      |              |
## 0            | 0.0170       | 0.1490 - 2.0867 q0
## 0.0170       | 0.0658       | 0.0438 + 4.1075 q0
## 0.0658       | +            | 0.3141
## -------------------------------------------------------
## 
## equation to get mx from qx and ax:
## mx = qx / ((ax-1)*qx + 1)
##
##
#qxax2mx <- function(ax,qx){
#  qx / ((ax-1)*qx + 1)
#}
#
#q02a0 <- function(q0, sex = "m"){
#  sex <- rep(sex, length(q0))
#  ifelse(sex == "m", 
#    ifelse(q0 < .0226, {0.1493 - 2.0367 * q0},
#      ifelse(q0 < 0.0785, {0.0244 + 3.4994 * q0},.2991)),
#    ifelse(q0 < 0.0170, {0.1490 - 2.0867 * q0},
#      ifelse(q0 < 0.0658, {0.0438 + 4.1075 * q0}, 0.3141))
#  )
#}
## using automagically generated equations
#m02a0Complex <- function(m0, sex = "m"){
#  sex <- rep(sex, length(m0))
#  ifelse(sex == "m", 
#    ifelse(m0 < qxax2mx(0.1493 - 2.0367 * 0.0226, 0.0226), {0.1493 - 2.0367 * (sqrt(887049049*m0^2+170140000*m0+100000000)-8507*m0-10000)/(40734 * m0)},
#      ifelse(m0 < qxax2mx(0.0244 + 3.4994 * 0.0785, 0.0785), {0.0244 + 3.4994 * -(sqrt(-81536279 * m0^2+12195000*m0+6250000)-2439*m0-2500)/(17497*m0)},.2991)),
#    ifelse(m0 < qxax2mx(0.1490 - 2.0867 * 0.0170, 0.0170), {0.1490 - 2.0867 * (5*sqrt(9071001 * m0^2 + 1702000*m0+1000000)-4255*m0-5000)/(20867*m0)},
#      ifelse(m0 < qxax2mx(0.0438 + 4.1075 * 0.0658, 0.0658), {0.0438 + 4.1075 * -(sqrt(-387892039 * m0^2 + 47810000*m0+25000000)-4781*m0-5000)/(41075*m0)}, 0.3141))
#  )
#}
## imputing rounded answers, changing the slope to match that calculated below...
#m02a0Simplified <- function(m0, sex = "m"){
#  sex <- rep(sex, length(m0))
#  ifelse(sex == "m", 
#    ifelse(m0 < 0.02306749, {0.1493 - 1.999268 * m0},
#      ifelse(m0 < 0.08307058, {0.0244 + 3.263742 * m0},.2991)),
#    ifelse(m0 < 0.01726011, {0.1490 - 2.057638 * m0},
#      ifelse(m0 < 0.06891019, {0.0438 + 3.887712 * m0}, 0.3141))
#  )
#}
#
## get new slopes...
## males
## first breakpoint
#m1 <- .001
#a1 <- 2.0367 * (sqrt(887049049*m1^2+170140000*m1+100000000)-8507*m1-10000)/(40734 * m1)
#m2 <- .02
#a2 <- 2.0367 * (sqrt(887049049*m2^2+170140000*m2+100000000)-8507*m2-10000)/(40734 * m2)
## slope, rise over run
#(a2-a1)/(m2-m1)
## second breakpf.m2a<- function(m,g,h){
#
#f.finda0 <-function(m0, sex="m"){
#
#
#m1 <- .03
#a1 <- 3.4994 * -(sqrt(-81536279 * m1^2+12195000*m1+6250000)-2439*m1-2500)/(17497*m1)
#m2 <- .07
#a2 <- 3.4994 * -(sqrt(-81536279 * m2^2+12195000*m2+6250000)-2439*m2-2500)/(17497*m2)
## slope
#(a2-a1)/(m2-m1)
## -------------------------
## females
## first breakpoint
#m1 <- .001
#a1 <- 2.0867 * (5*sqrt(9071001 * m1^2 + 1702000*m1+1000000)-4255*m1-5000)/(20867*m1)
#m2 <- .015
#a2 <- 2.0867 * (5*sqrt(9071001 * m2^2 + 1702000*m2+1000000)-4255*m2-5000)/(20867*m2)
## slope
#(a2-a1)/(m2-m1)
## second breakpoint
#m1 <- .02
#a1 <- 4.1075 * -(sqrt(-387892039 * m1^2 + 47810000*m1+25000000)-4781*m1-5000)/(41075*m1)
#m2 <- .06
#a2 <- 4.1075 * -(sqrt(-387892039 * m2^2 + 47810000*m2+25000000)-4781*m2-5000)/(41075*m2)
## slope
#(a2-a1)/(m2-m1)
##############################################################
## test values
#q0 <- seq(0,.1,by = .0001)
#m0m <- qxax2mx(q02a0(q0,"m"),q0)
#m0f <- qxax2mx(q02a0(q0,"f"),q0)
##plot(m0m,type='l')
##lines(m0f,col="red")
## the complex one matches
#plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
#
#plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
#lines(q0,m02a0Complex(m0m,"m"), col = "red",lty=2,lwd=2)
#
#lines(q0,q02a0(q0,"f"), col = "black")
#lines(q0,m02a0Complex(m0f,"f"), col = "blue",lty=2,lwd=2)
## the simplified one is rather off...WTF?
## ---------------------------------------------------
#plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
#lines(q0,m02a0Complex(m0m,"m"), col = "red",lty=2)
#lines(q0,m02a0Simplified(m0m,"m"), col = "red",lty=4,lwd=3)
#
#lines(q0,q02a0(q0,"f"), col = "black",lty=2)
#lines(q0,m02a0Complex(m0f,"f"), col = "blue",lty=2,lwd=2)
#lines(q0,m02a0Simplified(m0f,"f"), col = "blue",lty=4,lwd=3)
#
## ----------------------------------------------------
## second attempt
#m02a0Complex2 <- function(m0, sex = "m"){
#  sex <- rep(sex, length(m0))
#  ifelse(sex == "m", 
#    ifelse(m0 < qxax2mx(q02a0(0.0226), 0.0226), {0.1493 - 2.0367 * (sqrt(887049049*m0^2+170140000*m0+100000000)-8507*m0-10000)/(40734 * m0)},
#      ifelse(m0 < qxax2mx(q02a0(0.0785), 0.0785), {0.0244 + 3.4994 * -(sqrt(-81536279 * m0^2+12195000*m0+6250000)-2439*m0-2500)/(17497*m0)},.2991)),
#    ifelse(m0 < qxax2mx(q02a0(0.0170), 0.0170), {0.1490 - 2.0867 * (5*sqrt(9071001 * m0^2 + 1702000*m0+1000000)-4255*m0-5000)/(20867*m0)},
#      ifelse(m0 < qxax2mx(q02a0(0.0658), 0.0658), {0.0438 + 4.1075 * -(sqrt(-387892039 * m0^2 + 47810000*m0+25000000)-4781*m0-5000)/(41075*m0)}, 0.3141))
#  )
#}
#m02a0Simplified2 <- function(m0, sex = "m"){
#  sex <- rep(sex, length(m0))
#  ifelse(sex == "m", 
#    ifelse(m0 < 0.02306737, {0.1493 - 1.993632 * m0},
#      ifelse(m0 < 0.0830706, {0.0244 + 3.2601 * m0},.2991)),
#    ifelse(m0 < 0.01725977, {0.1490 -  2.053429 * m0},
#      ifelse(m0 < 0.06919348, {0.0438 + 3.880473 * m0}, 0.3141))
#  )
#}
#scal1 <- qxax2mx(q02a0(0.0226), 0.0226)/0.0226
#scal2 <- qxax2mx(q02a0(0.0785), 0.0785)/0.0785
#
## get new slopes for the above
## males
## first breakpoint
#m1 <- .001
#a1 <- 2.0367 * (sqrt(887049049*m1^2+170140000*m1+100000000)-8507*m1-10000)/(40734 * m1)
#m2 <- 0.02306737
#a2 <- 2.0367 * (sqrt(887049049*m2^2+170140000*m2+100000000)-8507*m2-10000)/(40734 * m2)
## slope
#(a2-a1)/(m2-m1)
#0.02306737
## second breakpoint
#m1 <- 0.02306737
#a1 <- 3.4994 * -(sqrt(-81536279 * m1^2+12195000*m1+6250000)-2439*m1-2500)/(17497*m1)
#m2 <- 0.0830706
#a2 <- 3.4994 * -(sqrt(-81536279 * m2^2+12195000*m2+6250000)-2439*m2-2500)/(17497*m2)
## slope
#(a2-a1)/(m2-m1)
#a1 - m1 * (a2-a1)/(m2-m1)
#
## -------------------------
## females
## first breakpoint
#m1 <- .001
#a1 <- 2.0867 * (5*sqrt(9071001 * m1^2 + 1702000*m1+1000000)-4255*m1-5000)/(20867*m1)
#m2 <- 0.01725977
#a2 <- 2.0867 * (5*sqrt(9071001 * m2^2 + 1702000*m2+1000000)-4255*m2-5000)/(20867*m2)
## slope
#(a2-a1)/(m2-m1)
## second breakpoint
#m1 <- 0.01725977
#a1 <- 4.1075 * -(sqrt(-387892039 * m1^2 + 47810000*m1+25000000)-4781*m1-5000)/(41075*m1)
#m2 <- 0.06919348
#a2 <- 4.1075 * -(sqrt(-387892039 * m2^2 + 47810000*m2+25000000)-4781*m2-5000)/(41075*m2)
## slope
#(a2-a1)/(m2-m1)
#plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
#lines(q0,m02a0Complex2(m0m,"m"), col = "red",lty=2)
#lines(q0,m02a0Simplified2(m0m,"m"), col = "red",lty=4,lwd=3)
#
#lines(q0,q02a0(q0,"f"), col = "black",lty=2)
#lines(q0,m02a0Complex2(m0f,"f"), col = "blue",lty=2,lwd=2)
#lines(q0,m02a0Simplified2(m0f,"f"), col = "blue",lty=4,lwd=3)
#
#
#
## Carl's solution/
##f.m2a<- function(m,g,h){
##  t2 = (g * m);
##  t3 = (m ^ 2);
##  t8 = (g ^ 2);
##  t13 = sqrt((t3 - 2 * m - 2 * g * t3 + 1 + 2 * t2 +
##        t8 * t3 + 4 * t3 * h));
##  t16 = 0.1e1 / m * (m - 0.1e1 + t2 + t13) / 0.2e1;
##  return(t16)
##}
#f.m2q<- function(m,g,h){
#  t1 = (m * g);
#  t2 = (m ^ 2);
#  t6 = (g ^ 2);
#  t12 = sqrt((t2 - 2 * g * t2 + 2 * m + t6 * t2 - 2 * t1 + 1 - 4 * h * t2));
#  t19 = -0.1e1 / h / m * (-m + t1 - 0.1e1 + t12) / 0.2e1;
#  return(t19)
#}
#
#
#
##f.finda0 <-function(m0, sex="m"){
##  if(sex=="m"){
##    q.changepoints<- c( 0, 0.0226, 0.0785 );
##    g.changepoints<- c( 0.1493, 0.0244, .2991 );
##    h.changepoints<- c(-2.0367, 3.4994, 0 );
##  } else {
##    q.changepoints<- c( 0, 0170, 0.0658 );
##    g.changepoints<- c( 0.1490, 0.0438, .3141 );
##    h.changepoints<- c(-2.0867, 4.1075, 0 );
##  }
##  
##  a.changepoints <- g.changepoints + h.changepoints * q.changepoints
##  m.changepoints <- q.changepoints / (1 - (1-a.changepoints)*q.changepoints)
##  
##  
##  isel <- findInterval(m0, m.changepoints)
##  
##  a0 <-switch(isel,
##    f.m2a(m0,g.changepoints[1],h.changepoints[1]),
##    f.m2a(m0,g.changepoints[2],h.changepoints[2]),
##    f.m2a(m0,g.changepoints[3],h.changepoints[3])
##  )
## 
##}
#
#f.finda02 <-function(m0, sex="m"){
#  if(sex=="m"){
#    q.changepoints<- c( 0, 0.0226, 0.0785 );
#    g.changepoints<- c( 0.1493, 0.0244, .2991 );
#    h.changepoints<- c(-2.0367, 3.4994, 0 );
#  } else {
#    q.changepoints<- c( 0, .0170, 0.0658 );
#    g.changepoints<- c( 0.1490, 0.0438, .3141 );
#    h.changepoints<- c(-2.0867, 4.1075, 0 );
#  }
#  
#  a.changepoints <- g.changepoints + h.changepoints * q.changepoints
#  m.changepoints <- q.changepoints / (1 - (1-a.changepoints)*q.changepoints)
#  
#  isel <- findInterval(m0, m.changepoints)
#  
#  ifelse(isel == 1,g.changepoints[1]+h.changepoints[1]*f.m2q(m0,g.changepoints[1],h.changepoints[1]),
#    ifelse(isel == 2, g.changepoints[2]+h.changepoints[2]*f.m2q(m0,g.changepoints[2],h.changepoints[2]),
#      g.changepoints[3]+h.changepoints[3]*f.m2q(m0,g.changepoints[3],h.changepoints[3])))
#}
#
#
#
#
#m02a0Complex2
#
#m0 <- m0m
#plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
#lines(q0,m02a0(m0m,"m"), col = "blue",lty=4,lwd=3)
#lines(q0,f.finda02(m0m,sex="m"), col = "red",lty=2)
#
#lines(q0,m02a0Complex(m0m,"m"), col = "red",lty=2)
#lines(q0,m02a0Simplified(m0m,"m"), col = "red",lty=4,lwd=3)
#lines(q0,m02a0(m0m,"m"), col = "blue",lty=4,lwd=3)
#
#
#lines(q0,q02a0(q0,"f"), col = "black",lty=2)
#lines(q0,m02a0Complex(m0f,"f"), col = "blue",lty=2,lwd=2)
#lines(q0,m02a0Simplified(m0f,"f"), col = "blue",lty=4,lwd=3)
#lines(q0,m02a0Simplified(m0f,"f"), col = "blue",lty=4,lwd=3)
#
#
#f.m2q<- function(m,g,h){
#  t1 = (m * g);
#  t2 = (m ^ 2);
#  t6 = (g ^ 2);
#  t12 = sqrt((t2 - 2 * g * t2 + 2 * m + t6 * t2 - 2 * t1 + 1 - 4 * h * t2));
#  t19 = -0.1e1 / h / m * (-m + t1 - 0.1e1 + t12) / 0.2e1;
#  return(t19)
#}
#m02a0 <- function(m0,sex="m"){
#  sex <- rep(sex,length(m0))
#  ifelse(sex == "m",
#      ifelse(m0 < 0.02306737, 0.1493 - 2.0367 * f.m2q(m0, 0.1493, -2.0367),
#      ifelse(m0 < 0.0830706, 0.0244 + 3.4994 * f.m2q(m0, 0.0244, 3.4994), .2991)),
#    ifelse(m0 < 0.01725977, 0.1490 - 2.0867 * f.m2q(m0, 0.1490, -2.0867),
#    ifelse(m0 < 0.06919348, 0.0438 + 4.1075 * f.m2q(m0, 0.0438, 4.1075), 0.3141))
#    )
#}
##plot(q0,q02a0(q0,"m"),type = 'l',ylim = c(0.05,.35))
##lines(q0,m02a0(m0m,"m"), col = "blue",lty=4,lwd=3)
##lines(q0,f.finda02(m0m,sex="m"), col = "red",lty=2)
##
##plot(q0,q02a0(q0,"f"),type = 'l',ylim = c(0.05,.35))
##lines(q0,m02a0(m0f,"f"), col = "blue",lty=4,lwd=3)
##lines(q0,f.finda02(m0f,sex="f"), col = "red",lty=2)
#
## code not activated just in case package tries to get built...
#}

# Males
# lower q0     |   upper q0   | formula
# 0            | 0.0226       | 0.1493 - 2.0367 q0
# 0.0226       | 0.0785       | 0.0244 + 3.4994 q0
# 0.0785       |  +           | .2991
# ________________________________________
# Females      |              |
# 0            | 0.0170       | 0.1490 - 2.0867 q0
# 0.0170       | 0.0658       | 0.0438 + 4.1075 q0
# 0.0658       | +            | 0.3141
# -------------------------------------------------------

# -------------------------------------------------
# do fake test against iterative:
# assume q01 = m0, then iterate:
#m0 <- .02
#q0 <- m0
#for (i in 1:10){
#  a0 <- AKq02a0(q0)
#  q0 <- m0 / (1 + (1 - a0) * m0)    
#}
## I consider this verified:
#AKm02q0(m0, 0.1493, -2.0367)-q0
#AKm02a0(m0)-a0
