\name{Write_Births}
\alias{Write_Births}
\title{\code{Write_Births} a function to prepare, format and write \code{Births.txt}.}
\usage{
  Write_Births(WORKING = getwd(), STATSFOLDER = "RSTATS",
    LDBPATH = NULL, MPVERSION = 5, XXX = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{STATSFOLDER}{the folder name where output is to be
  written to (not a full path). Default \code{"RSTATS"}.}

  \item{LDBPATH}{in case the LexisDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "LexisDB")}}

  \item{PVERSION}{5 or 6. Default 5. Here this only affects
  file headers.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  last path part.}
}
\value{
  function called for its side effect of creating the file
  \code{Births.txt}. No value returned.
}
\description{
  This function pulls birth count data from the Lexis
  Database (age = 0, lexis = 1) and sets up the standard
  output with columns \code{Year}, \code{Female},
  \code{Male}, and \code{Total}. This function makes no
  consequential calls to other functions, and is not called
  by any other function, nor does it require another
  function to have been run previously.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

