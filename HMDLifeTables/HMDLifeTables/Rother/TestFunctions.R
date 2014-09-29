 
# Author: triffe
###############################################################################

# Functions for testing output in RSTATS vs STATS

# test contents of STATS versus RSTATS
test_Contents <- function(WORKING = getwd()){
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  all.files     <- union(STATS.files, RSTATS.files)
  all.files.m   <- matrix(FALSE, nrow = length(all.files), ncol = 2, dimnames = list(all.files, c("STATS","RSTATS")))
  all.files.m[STATS.files, "STATS"]   <- TRUE
  all.files.m[RSTATS.files, "RSTATS"] <- TRUE
  sums <- rowSums(all.files.m)
  if (any(sums == 1)){
    cat("\n\nComparing contents of /STATS/ vs /RSTATS/ \n All files match except:\n")
    return(all.files.m[sums == 1, ])
  } else {
    return(TRUE)
  }
}
# test_Contents("/data/commons/hmd/HMDWORK/BEL")

# tests ages 0-80 of the first 5 columns of any 'ltper' output
test_ltper <- function(WORKING = getwd(), ignore.pattern = NULL){
  # since we know Kannisto smoothing will generally differ, we
  # test only mx, qx, ax up to age 80, as these must be equal.
  # *we allow, but report, differences of 1 unit in the last decimal place
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  STATS.ltper   <- grep(grep(STATS.files, pattern = "ltper", value = TRUE), pattern = ".txt", value = TRUE)
  RSTATS.ltper  <- grep(grep(RSTATS.files, pattern = "ltper", value = TRUE), pattern = ".txt", value = TRUE)
  
  files.to.test <- intersect(STATS.ltper, RSTATS.ltper)
  if (!is.null(ignore.pattern)){
    for (i in 1:length(ignore.pattern)){
      files.to.test <- grep(files.to.test, pattern = ignore.pattern[i], value = TRUE, invert = TRUE)
    }
  }
  
  # read in both sets of ltper files
  All.data <- lapply(files.to.test, function(fname, WORKING){
      ages.keep <- as.character(0:80)
      cols.keep <- c("Year", "Age", "mx", "qx", "ax")
      Rfile <- read.table(file.path(WORKING, "RSTATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      Mfile <- read.table(file.path(WORKING, "STATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      list(Rfile = Rfile[Rfile$Age %in% ages.keep, cols.keep], Mfile = Mfile[Mfile$Age %in% ages.keep, cols.keep])
    }, WORKING = WORKING)
  names(All.data) <- files.to.test

  equality.test <- do.call(rbind,lapply(All.data, function(x){
        c(Year = identical(x[["Rfile"]][["Year"]], x[["Mfile"]][["Year"]]),
          Age = identical(x[["Rfile"]][["Age"]], x[["Mfile"]][["Age"]]),
          mx = identical(x[["Rfile"]][["mx"]], x[["Mfile"]][["mx"]]),
          qx = identical(x[["Rfile"]][["qx"]], x[["Mfile"]][["qx"]]),
          ax = identical(x[["Rfile"]][["ax"]], x[["Mfile"]][["ax"]]))
      }
    )
  )
  if (!all(equality.test)){
    cat("\n\nComparing all files with 'ltper' in their names, mx, qx, ax columns up to age 80\nOnly following files showed any discrepancy (extreme detail comparison):\n")
    return(equality.test[rowSums(equality.test) != 5, colSums(equality.test) != nrow(equality.test)])
  } else {
    return(TRUE)
  }
}

# tests ages 0-80 of the first 5 columns of any 'ltcoh' output
test_ltcoh <- function(WORKING = getwd(), ignore.pattern = NULL){
  # since we know Kannisto smoothing will generally differ, we
  # test only mx, qx, ax up to age 80, as these must be equal.
  # *we allow, but report, differences of 1 unit in the last decimal place
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  STATS.coh   <- grep(grep(STATS.files, pattern = "ltcoh", value = TRUE), pattern = ".txt", value = TRUE)
  RSTATS.coh  <- grep(grep(RSTATS.files, pattern = "ltcoh", value = TRUE), pattern = ".txt", value = TRUE)
  
  files.to.test <- intersect(STATS.coh, RSTATS.coh)
  if (!is.null(ignore.pattern)){
    for (i in 1:length(ignore.pattern)){
      files.to.test <- grep(files.to.test, pattern = ignore.pattern[i], value = TRUE, invert = TRUE)
    }
  }
  
  # read in both sets of ltper files
  All.data <- lapply(files.to.test, function(fname, WORKING){
     
      Rfile <- read.table(file.path(WORKING, "RSTATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      Mfile <- read.table(file.path(WORKING, "STATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      list(Rfile = Rfile, Mfile = Mfile)
    }, WORKING = WORKING)
  names(All.data) <- files.to.test
 
  equality.test <- do.call(rbind,lapply(All.data, function(x){
        c(Year = identical(x[["Rfile"]][["Year"]], x[["Mfile"]][["Year"]]),
          Age = identical(x[["Rfile"]][["Age"]], x[["Mfile"]][["Age"]]),
          mx = identical(x[["Rfile"]][["mx"]], x[["Mfile"]][["mx"]]),
          qx = identical(x[["Rfile"]][["qx"]], x[["Mfile"]][["qx"]]),
          ax = identical(x[["Rfile"]][["ax"]], x[["Mfile"]][["ax"]]),
          lx = identical(x[["Rfile"]][["lx"]], x[["Mfile"]][["lx"]]),
          dx = identical(x[["Rfile"]][["dx"]], x[["Mfile"]][["dx"]]),
          Lx = identical(x[["Rfile"]][["Lx"]], x[["Mfile"]][["Lx"]]),
          Tx = identical(x[["Rfile"]][["Tx"]], x[["Mfile"]][["Tx"]]),
          ex = identical(x[["Rfile"]][["ex"]], x[["Mfile"]][["ex"]])
        )
      }
    )
  )
  if (!all(equality.test)){
    cat("\n\nComparing all files with 'ltcoh' in their names, all columns, all ages\nOnly following files-columns showed any discrepancy (extreme detail comparison):\n")
    return(equality.test[rowSums(equality.test) != ncol(equality.test), colSums(equality.test) != nrow(equality.test)])
  } else {
    return(TRUE)
  }
}

# function returns number of differences by column for a particular file (pattern)
# file.pattern = "Population.txt"
test_file_nr_diff  <- function(WORKING = getwd(), file.pattern = "ltper"){
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  STATS.f       <- grep(grep(STATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  RSTATS.f      <- grep(grep(RSTATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  
  files.to.test <- intersect(STATS.f, RSTATS.f)
  
  # read in both sets of ltper files
  All.data      <- lapply(files.to.test, function(fname, WORKING){
      Rfile <- read.table(file.path(WORKING, "RSTATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      Mfile <- read.table(file.path(WORKING, "STATS", fname), 
        header = TRUE,
        skip = 2, 
        stringsAsFactors = FALSE)
      list(Rfile = Rfile, Mfile = Mfile)
    }, WORKING = WORKING)
  names(All.data) <- files.to.test
  equality.test   <- do.call(rbind, lapply(All.data, function(x){
        round(colSums(x[["Rfile"]] != x[["Mfile"]]) / nrow(x[["Rfile"]]), digits = 4)
      }
    )
  )
  if (!all(equality.test == 0)){
    cat("\n\nComparing all files with specified name pattern\nProportion of rows per column with a discrepancy of any magnitude\n")
    return(equality.test)
  } else {
    return(TRUE)
  }
}

# opens a tkdiff device to take a closer look
tkdiff_file <- function(WORKING = getwd(), file = "mltper_1x1"){
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  STATS.f       <- grep(grep(STATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  RSTATS.f      <- grep(grep(RSTATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  
  files.to.test <- intersect(STATS.f, RSTATS.f)
  
  system(paste0("tkdiff ", file.path(WORKING, "RSTATS", file), " ", file.path(WORKING, "STATS", file)))
}

# opens a tkdiff device to take a closer look
tkdiff_file <- function(WORKING = getwd(), file = "mltper_1x1"){
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  STATS.f       <- grep(grep(STATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  RSTATS.f      <- grep(grep(RSTATS.files, pattern = file.pattern, value = TRUE), pattern = ".txt", value = TRUE)
  
  files.to.test <- intersect(STATS.f, RSTATS.f)
  cat("\nOpening tkdiff, can take a few secs...\n")
  system(paste0("tkdiff ", file.path(WORKING, "RSTATS", file), " ", file.path(WORKING, "STATS", file)))
}

# lists files eligible for comparison
list_files <- function(WORKING = getwd()){
  STATS.files   <- list.files(file.path(WORKING, "STATS"))
  RSTATS.files  <- list.files(file.path(WORKING, "RSTATS"))
  intersect(STATS.files, RSTATS.files)
}