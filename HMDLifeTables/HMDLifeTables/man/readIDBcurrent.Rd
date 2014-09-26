\name{readIDBcurrent}
\alias{readIDBcurrent}
\title{\code{readIDBcurrent()} read in any IDB file to a data.frame}
\usage{
  readIDBcurrent(XXX = "SWE", item = "death", old = FALSE,
    full.path = NULL)
}
\arguments{
  \item{XXX}{character string; HMD country short code}

  \item{item}{character string; one of
  \code{c("death","pop","birth","tadj","monthly")}}

  \item{old}{logical should the file come from XXX_Old?}

  \item{full.path}{optional. *If you're not at UCB, or if
  you're no working from the main HMD working directory
  \code{"~/HMDWORK"}, then give the full path to the file,
  including \code{".txt"}.}
}
\value{
  data.frame for the given input file.
}
\description{
  Character columns are read in as character. Missing
  values are assumed to be indicated with \code{"."}, due
  to convention. We assume there are column headers. File
  must be comma separated. All of things will hold for
  standard InputDB files.
}
\examples{
#bo <- readIDBcurrent("SWE","birth")
#bo <- resortBirths(bo)
}

