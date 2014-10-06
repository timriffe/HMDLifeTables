\name{AKm02a0_direct}
\alias{AKm02a0_direct}
\title{\code{AKm02a0_direct} estimates a0 using an analytic translation of the Andreev-Kinkade rule of thumb.}
\usage{
  AKm02a0_direct(m0, sex = "m")
}
\arguments{
  \item{m0}{a value or vector of values of m0, the death
  rate for age 0 infants.}

  \item{sex}{either "m" or "f"}
}
\value{
  a0, the estimated average age at death of those dying in
  the first year of life, either a single value or a vector
  of a_0 values.
}
\description{
  \code{AKm02a0_direct} is an auxiliary function that could
  be used by version 6 of the Methods Protocol. Instead we
  opt for a simple approximate solution. This function
  calls \code{AKm02q0()} (also deprecated) to help get the
  work done, since the HMD needed to adapt the
  Andreev-Kingkade formulas to work with the period
  lifetable flow.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

