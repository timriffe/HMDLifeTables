\docType{package}
\name{RLifeTables}
\alias{RLifeTables}
\alias{RLifeTables-package}
\title{The full suite of lifetable functions used by the HMD}
\description{
  RLifeTables provides the entirety of post LexisDB
  functionality required by the HMD for country updates and
  testing. There are four kinds of functions in this
  package: 1) top-level functions, used for creating
  in-memory R \code{data.frame} output for the various time
  and sex configurations of period and cohort lifetables.
  These include \code{ltper_AxN()}, \code{ltcoh_AxN()},
  their both-sex variants, and the omnibus
  \code{RunHMDCountry()}. 2) auxiliary functions that do
  minor things. These have been written in the name of
  modularity and include functions \code{Abrdige()},
  \code{YearAgg()}, code{CDa0()} and others. These can be
  useful, but need not be handled at the top level for a
  mere country update. Top level functions call these as
  necessary. 3) data preparation functions such,
  \code{getCohortComponents()} and
  \code{getPeriodComponents()}, which feed data into
  downstream functions in a standard and useful format
  (territorial adjustments in-place) and 4) output
  formatting functions for preparing the .txt output files
  in a STATS folder (here RSTATS, unless otherwise
  specified by the user).
}
\details{
  Generating lifetables. For a country update, try: (if
  you're at UCB)
}
\examples{
\dontrun{
   RunHMDCountry(WORKING = "/hdir/0/hmd/HMDWORK/ISL") # Does Iceland
}
Assuming the standard HMD folder structure, LexisDB and InputDB file formats, the
entire body of HMD output will be generated. Intermediate files (unrounded) are also
produced as R binary (\\code{.Rdata}) files in the \\code{WORKING/Rbin/} folder. These
can be useful either as an alternative final product or for testing. \\code{RunHMDCountry()}
fully erases the specified \\code{STATSFOLDER} and \\code{Rbin} folder to ensure that
it is always working with the most recent files. In this way, each run is a full and
fresh run.

While this is the only necessary function for typical country updates, all other functions
are directly available to users, such that one could dynamically change arguments to produce
different output. Imagine for instance that you wanted cohort lifetables for a wider range of
years than the default; you could specify:
}
\references{
  Wilmoth, J.R. and Andreev, K. and Jdanov, D. and Glei, DA
  (2007). Methods protocol for the human mortality
  database. University of California, Berkeley, and Max
  Planck Institute for Demographic Research, Rostock,
  [version 31/05/2007] \url{http://mortality.org}.
}

