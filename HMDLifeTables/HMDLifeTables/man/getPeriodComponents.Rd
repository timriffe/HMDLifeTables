\name{getPeriodComponents}
\alias{getPeriodComponents}
\title{\code{getPeriodComponents} a function to read in and prepare basic data components for period calculations}
\usage{
  getPeriodComponents(WORKING = getwd(), sex = "f",
    OPENAGE = 110, save.bin = TRUE, XXX = NULL,
    LDBPATH = NULL, IDBPATH = NULL)
}
\arguments{
  \item{WORKING}{path to working directory, which typically
  ends with the HMD country abbreviation. Default
  \code{getwd()}.}

  \item{sex}{either \code{"m"} or \code{"f"}}

  \item{save.bin}{logical. Default \code{TRUE}. Should the
  output be saved as e.g.
  \code{Rbin/f_periodComponents.Rdata.Rdata} /
  \code{Rbin/m_periodComponents.Rdata.Rdata} as well?
  Appropriate name is derived automatically.}

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
  a long \code{data.frame} with 5 columns:
  \code{"Year"},\code{"Age"}, \code{"dl"}, \code{"du"} and
  \code{"pop1"} (Jan 1st), \code{"pop2"} (Dec 31st).
}
\description{
  This function modularizes the reading in, staggering, and
  post-territorial adjusting of the 4 basic (v5) data
  components of period calculations: pop1, pop2, dl and du,
  and sticks them in a long-format \code{data.frame}
  amenable to reshaping in downstream functions.
}
\details{
  This function is called by \code{ltper_Axn()},
  \code{ltperBoth_AxN()} and
  \code{Exposures_Deaths_Mx_AxN()} when needed, though
  these functions can also take the output of
  \code{getPeriodComponents()} as input for in-memory data
  object passing (faster). The function
  \code{RunHMDCountry()} uses it this way. Another
  alternative is to save off the output of this function
  and read it in using \code{load()}, which is another
  built-in feature of the above functions. This function
  calls \code{perTadj()} when pertinent.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

