\name{Write_Deaths_lexis}
\alias{Write_Deaths_lexis}
\title{\code{Write_Deaths_lexis} a function to prepare and save the file \code{deaths_lexis.txt}.}
\usage{
  Write_Deaths_lexis(WORKING = getwd(), males = NULL,
    females = NULL, STATSFOLDER = "RSTATS", OPENAGE = 110,
    MPVERSION = 5, XXX = NULL, LDBPATH = LDBPATH)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{STATSFOLDER}{the folder name where output is to be
  written to (not a full path). Default \code{"RSTATS"}.}

  \item{OPENAGE}{the desired open age. Default value is
  110.}

  \item{MPVERSION}{5 or 6. Default 5. Here this only
  affects file headers.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  last path part.}

  \item{LDBPATH}{in case the LexisDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "LexisDB")}}
}
\value{
  function called for its side effect of creating the files
  \code{Population.txt} or \code{Population5.txt}. No value
  returned.
}
\description{
  This function is essentially a stand-alone that prepares
  \code{deaths_lexis.txt} straight from LexisDB files. At
  this time is makes no indispensable calls to other
  functions, although it could be made to use
  \code{getPeriodComponents} as other functions do, but
  this is not a high priority.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

