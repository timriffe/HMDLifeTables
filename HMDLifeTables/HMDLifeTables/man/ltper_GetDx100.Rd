\name{ltper_GetDx100}
\alias{ltper_GetDx100}
\title{\code{ltper_GetDx100} determines the age at which Kannisto smoothed mx values should be imputed over the raw Mx for each year of data.}
\usage{
  ltper_GetDx100(WORKING = getwd(), perComp, OPENAGE = 110,
    N = 1, XXX = NULL, LDBPATH = NULL, IDBPATH = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{OPENAGE}{the desired open age. Default value is
  110}

  \item{N}{year interval: 1, 5 or 10. Other intervals would
  also work in theory.}

  \item{XXX}{the HMD country abbreviation. If left
  \code{NULL}, this is extracted from \code{WORKING} as the
  last path part.}

  \item{LDBPATH}{LDBPATH in case the LexisDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "LexisDB")}.}

  \item{IDBPATH}{in case the InputDB is not in
  \code{WORKING} (local testing), the full path to the
  LexisDB folder. If left as \code{NULL} it is assumed to
  be \code{file.path(WORKING, "InputDB")}.}
}
\value{
  a vector with the age index (age + 1) for each year where
  smoothed mx values should be imputed.
}
\description{
  Smoothed mx values are imputed starting at the lowest age
  where male or female deaths drop to at most 100 with a
  minimum age of 80 and a maximum age of 95. The same age
  is used for both males and females.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

