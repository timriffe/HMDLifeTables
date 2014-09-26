\name{cohTadj}
\alias{cohTadj}
\title{\code{cohTadj} a function for doing lifetable territorial adjustments for cohort data}
\usage{
  cohTadj(LDBobj, tadj, sex)
}
\arguments{
  \item{LDBobj}{Lexis database object, with proper column
  names, to be passed in by \code{getCohortComponents()}.}

  \item{tadj}{a territorial adjustment object, as read
  straight from the Input Database.}

  \item{sex}{\code{"m"} or \code{"f"}}
}
\value{
  returns LDBobj back again, having 'un'-adjusted deaths
  and pop counts for pertinent years.
}
\description{
  cohort territorial adjustments 'un'adjust death counts
  from the year t+1 upper triangle as well as bring January
  1 population estimates from year t + 1 to December 31 so
  as to to have exposure and deaths both refer to year t.
  It's weird that this step needs to happen in the
  lifetable programs in my (TR) opinion- this should
  eventually be moved to the LexisDB programs, which would
  ideally return a more comprehensive data object. This
  function is called only by \code{getCohortComponents()}.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

