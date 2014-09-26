\name{cExposures_and_cMx_AxN}
\alias{cExposures_and_cMx_AxN}
\title{\code{cExposures_and_cMx_AxN} a function to prepare unrounded cExposures and cMx data.frames}
\usage{
  cExposures_and_cMx_AxN(WORKING = getwd(), males = NULL,
    females = NULL, abridged = FALSE, N = 1, OPENAGE = 110,
    min.years = 30, MPVERSION = 5, save.bin = TRUE,
    XXX = NULL, LDBPATH = NULL, IDBPATH = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
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
  Default 1 (no aggregatation)}

  \item{OPENAGE}{the desired open age. Default value is
  110}

  \item{min.years}{the minimum number of data years that we
  need in order to publish anything for a cohort. Default
  30.}

  \item{MPVERSION}{5 or 6. Default 5. Here this only
  affects the selection of which cohorts to show. For v5,
  matlab programs had an out-of-protocol but reasonable
  requirement of \code{min.years} continuous years of at
  most age 100. v6 will take the \code{min.years}
  consecutive years literally unless changed in the MP.}

  \item{save.bin}{logical. Default \code{TRUE}. Should the
  output be saved as e.g. \code{Rbin/cExposures_1x1.Rdata}
  as well. Appropriate name is derived systematically. In
  this case both objects are saved separately.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  past path part.}

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
  a list of two \code{data.frame}s of the full cMx and
  cExposures output, in long format. These have 5 columns:
  \code{Year}, \code{Age}, \code{Female}, \code{Male} and
  \code{Total}. Age is a formatted character vector, while
  the other columns are numeric.
}
\description{
  IFF cohort measures are published for country \code{XXX},
  this function will prepare the \code{data.frame}s for cMx
  and cExposures, saving them to \code{/Rbin/} in order to
  be used as such or else to be formatted for publishing
  later by WriteExposuresDeathsMx().
}
\details{
  This function has no native check for whether or not
  results should be run. That check needs to happen before
  this function runs. In case there aren't a sufficient
  number of cohorts with \code{min.years} worth of data, no
  warning is given at this time. If you run this function
  on a country with less than \code{min.years} worth of
  period data, it will throw an uninformative native R
  error, as no check is done.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

