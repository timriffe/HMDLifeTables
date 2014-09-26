\name{Exposures_Deaths_Mx_AxN}
\alias{Exposures_Deaths_Mx_AxN}
\title{\code{Exposures_Deaths_Mx_AxN} a function to prepare unrounded Exposures, Deaths and Mx data.frames}
\usage{
  Exposures_Deaths_Mx_AxN(WORKING = getwd(), males = NULL,
    females = NULL, abridged = FALSE, N = 1, OPENAGE = 110,
    MPVERSION = 5, save.bin = TRUE, XXX = NULL,
    LDBPATH = NULL, IDBPATH = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation.  Default
  \code{getwd()}.}

  \item{males}{the cohort components \code{data.frame}
  object, as created by \code{getCohortComponents()} for
  males. If left as \code{NULL},
  \code{getCohortComponents()} is called.}

  \item{females}{the cohort components \code{data.frame}
  object, as created by getCohortComponents() for females.
  If left as \code{NULL}, \code{getCohortComponents()} is
  called.}

  \item{abridged}{logical. Default \code{FALSE}.}

  \item{N}{integer number of years to aggregate 1, 5 or 10.
  Default 1 (no aggregatation).}

  \item{OPENAGE}{the desired open age. Default value is
  110.}

  \item{MPVERSION}{5 or 6. Default 5. Here this only
  affects exposure calculation. For both versions, exposure
  calculations are performed by \code{Exposures_per()}.}

  \item{save.bin}{logical. Default \code{TRUE}. Should the
  output be saved as e.g. \code{Rbin/cExposures_1x1.Rdata}
  as well. Appropriate name is derived systematically. In
  this case both objects are saved separately.}

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
  a list of three \code{data.frame}s of the full Deaths, Mx
  and Exposures output, in long format. These have 5
  columns: \code{Year}, \code{Age}, \code{Female},
  \code{Male} and \code{Total}. Age is a formatted
  character vector, while the other columns are numeric.
}
\description{
  All countries produce Deaths_1x1, Deaths_1x5.. 1x10, 5x1,
  5x5 and 5x10, and likewise for Exposures and Mx
  estimates. Since Mx requires the former two in order to
  be calculated, and all items have a similar output
  construction, they are all produced in one function.
  Output is returned as a list of (long)
  \code{data.frame}s, and optionally saved as appropriately
  named .Rdata objects in the /Rbin folder.
}
\details{
  This function is a top-level function, and can be called
  independently. If arguments \code{males}, \code{females}
  are missing, they are created by calls to
  \code{getPeriodComponents()}
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

