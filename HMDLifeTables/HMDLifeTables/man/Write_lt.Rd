\name{Write_lt}
\alias{Write_lt}
\title{\code{Write_lt} a function to write formatted .txt files for both period and cohort lifetables.}
\usage{
  Write_lt(WORKING = getwd(), STATSFOLDER = "RSTATS",
    MPVERSION = 5, XXX = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{STATSFOLDER}{the folder name where output is to be
  written to (not a full path). Default \code{"RSTATS"}.}

  \item{MPVERSION}{5 or 6. Default 5. Here this only
  affects file headers.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  last path part.}
}
\value{
  function called for its side effect of creating the
  lifetable txt output files, e.g. \code{mltper_1x1.txt}
  and other time/sex configurations. No value returned.
}
\description{
  This function requires that all pertinent
  \code{ltper_AxN.Rdata} and \code{ltcoh_AxN()} data
  objects be present in the folder \code{WORKING/Rbin/} for
  males, females and both-sex tables. Objects are selected
  by grepping, which works for the present HMD file naming
  scheme. \code{.Rdata} files are read in, rounded,
  formatted and written back out into
  \code{WORKING/RSTATS/} by default, although this can be
  changed with the \code{STATSFOLDER} argument. This
  function must be called after ltper_AxN() (for all time
  intervals and sexes), which must have been run with the
  argument \code{save.bin = TRUE}.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

