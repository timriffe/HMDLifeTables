\name{AKm02a0}
\alias{AKm02a0}
\title{\code{AKm02a0} estimates a0 using the Andreev-Kinkade rule of thumb.}
\usage{
  AKm02a0(m0, sex = "m")
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
  \code{AKm02a0} is an auxiliary function used by version 6
  of the four HMD lifetable functions, \code{ltper_AxN()},
  \code{ltcoh_AxN()}, \code{ltperBoth_AxN()},
  \code{ltcohBoth_AxN()}. This function calls
  \code{AKm02q0()} to help get the work done, since the HMD
  needed to adapt the Andreev-Kingkade formulas to work
  with the period lifetable flow.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

