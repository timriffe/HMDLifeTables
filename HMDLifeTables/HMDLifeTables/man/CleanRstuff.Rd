\name{CleanRstuff}
\alias{CleanRstuff}
\title{\code{CleanRstuff} DELETE the folders Rbin, RSTATS and RCHECKS permanently. Only used when doing testing. Function not to be included in shared package.}
\usage{
  CleanRstuff(WORKING = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{NULL}, since this is a consequential function.}
}
\value{
  function used for its side effects. USE EXTREME CAUTION.
}
\description{
  Only used when doing testing. Function not to be included
  in shared package. Before running the full gamut of
  lifetable functions it is often desirable to have totally
  clean folders, since some functions use existing objects
  where possible. Likewise, the RCHECKS folder ought to be
  up to date, so this removes detritus. The lifetable
  functions overwrite output, so this function is not
  strictly necessary. The most important folder to have
  clear before a clean run (sometimes) would be Rbin. These
  folder names are fixed (Rbin, RSTATS, RCHECKS). Do this
  type of thing manually if you want anything else.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

