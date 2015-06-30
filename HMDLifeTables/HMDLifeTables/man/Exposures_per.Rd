\name{Exposures_per}
\alias{Exposures_per}
\title{\code{Exposures_per} a function to prepare either v5 or v6 period exposures.}
\usage{
  Exposures_per(WORKING = getwd(), pop1 = NULL,
    pop2 = NULL, dl = NULL, du = NULL, perComp = NULL,
    sex = "m", OPENAGE = 110, save.bin = TRUE,
    MPVERSION = 6, XXX = NULL, LDBPATH = NULL,
    IDBPATH = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{pop1}{optional. An 'age by year' matrix of Jan 1st
  population counts. Taken from \code{perComp} or derived
  if necessary. #'}

  \item{pop2}{optional. An 'age by year' matrix of Dec 31st
  population counts territorially adjusted back from the
  january 1st year t+1 back to year t. Taken from
  \code{perComp} or derived if necessary.}

  \item{dl}{optional. An 'age by year' matrix of lower
  triangle death counts. Taken from \code{perComp} or
  derived if necessary.}

  \item{du}{optional. An 'age by year' matrix of upper
  triangle death counts. Taken from \code{perComp} or
  derived if necessary.}

  \item{perComp}{optional. Output from
  \code{getPeriodComponents()}.}

  \item{sex}{\code{"m"} or \code{"f"}.}

  \item{OPENAGE}{the desired open age. Default value is
  110}

  \item{save.bin}{logical. Default \code{TRUE}. Should the
  output be saved as e.g. \code{Rbin/cohExposuresRaw.Rdata}
  as well.}

  \item{MPVERSION}{5 or 6. Version 5 exposures assume
  uniformity, v6 allows for non-uniformity across birthdays
  in the death distribution by taking information from
  monthly birth distributions.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  last path part.}

  \item{LDBPATH}{in case the LexisDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "LexisDB")}}

  \item{IDBPATH}{in case the InputDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "InputDB")}}
}
\value{
  an 'age by year' matrix of exposures. if in test mode, a
  list with matrices of both deaths and exposures is
  returned.
}
\description{
  This function is meant to be internal. It is called by
  \code{Exposures_Deaths_Mx_AxN()}, \code{ltper_AxN()} and
  \code{ltperBoth_AxN()}. If necessary, it calls
  \code{getPeriodComponents()}. Exposures from the
  appropriate MP version are returned in an 'age by year'
  matrix already cut down to \code{OPENAGE}.
}
\details{
  This function can take data arguments in-memory, from the
  \code{/Rbin/} folder, or else derive them by calling
  \code{getPeriodComponents()}. This function is separate
  because using monthly birth data for v6 exposures is a
  fair amount of code. Since a few functions need these
  exposures it's best to modularize and call where needed.
  This function uses a package dataset
  \code{monthDurations.rda} in order to save having to
  recompute it on each run- this is a matrix with the
  number of days in each month every year across several
  centuries, reaching way beyond the potential of the HMD
  to as to avoid future bugs. It gets cut to size as
  needed. This function also calls \code{AC2AP()}.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

