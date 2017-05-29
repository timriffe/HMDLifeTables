
 readTadj <- function(Tadj.path, log.file = stdout() ){
 
 # test file path
 
if(! file.exists(Tadj.path) ){
  stop(paste(Tadj.path,": invalid tadj.txt path, no such file"))
}

header.expect.tadj  <- c("PopName", "Year", "Age", "Area1", "Area2", "Sex", "RefCode", 
                         "Access", "Type", "Value", "NoteCode1", "NoteCode2", "NoteCode3", "LDB")
                         
first.line <- readLines( con = Tadj.path, n = 1)

#    Allow for the file to exist but be empty
if (length(first.line) == 0){
  return(Tadj = NULL)
}

#  match first line with expected Headers.  Error out if not a perfect match, "Missing or erroneous headers"
#  order of columns can differ in a CSV file and still be valid
headers.read <- strsplit( gsub(" ", "", first.line), "," )[[1]]
headers.missing <- setdiff( header.expect.tadj, headers.read )

LDBisMissing <- ifelse(headers.missing == "LDB", TRUE, FALSE)
if ! LDBisMissing || length(headers.missing == 0){
   stop(paste(Tadj.path, ": missing columns:", paste(headers.missing, collapse=", ") )
}
 
## reread the file, since it has proper headers, it can be a CSV file 
 # read in entire file using headers, with expected column type designations.  The columns can be in any order. LDB
 # is optional, but put it in
Tadj  <- read.table(Tadj.path, 
                     sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE,
                     colClasses="character")
                     
if (LDBisMissing){
    cat("\nThere was no LDB column in ", Tadj.path, 
            ".\nPreviously, for some reason, the tadj file did not have an LDB column, but it is supposed to now.  
              \nPlease add this to the original csv version of the file, as it will not be written out\n\n",
            file = log.file, append = TRUE) 
    stop("LDBisMissing")
}
                     
                     
 # integrity tests
 ## duplicate entries
 testit <- Tadj[, c("Year","Age","Sex","Type")]
 duplicated.rows <- duplicated(testit)
 if(any(duplicated.rows)){
    cat("\nThere duplicate entries in ", Tadj.path, file = log.file, append = TRUE)
    print("Duplicated Tadj entries:")
    print(Tadj[ duplicated.rows, ] )
    stop("Duplicated rows in Tadj")
}

TadjVx <- Tadj[ Tadj$Type == "Vx", ]
YearSex<- split(TadjVx$Age, c(TadjVx$Year, TadjVx$Sex)
AgesRequired <- 0:130  # InputDB Docs state that this is required

fMissingAges <- function(x){
  missing.ages <- AgesRequired[ is.na( match(AgesRequired, x) ) ]
  if (length(missing.ages) == 0) {
    missing.ages <- NULL  # use NULL for autodrop of list entries
  }
  return(missing.ages)  
}
  
TadjVx.missing.ages <- sapply( as.numeric(TadjVx$Ages),  fMissingAges)
iFindNulls <- sapply(TadjVx.missing.ages, is.na)

if (any( ! iFindNulls) ){
   print("Missing Age entries for Year.Sex from tadj.txt")
   print( TadjVx[ !iFindNulls, ] )
   stop( "Tadj: missing some ages")
}

return(Tadj)
}

 
 if (tadjTF){
    Tadj        <- read.table(file.path(InputDB.path, grep(pattern = "tadj", x = files.we.want, value = TRUE)), 
                     sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  }

  
   header.expect.tadj      <- c("PopName", "Year", "Age", "Area1", "Area2", "Sex", "RefCode", 
      "Access", "Type", "Value", "NoteCode1", "NoteCode2", "NoteCode3", "LDB")

         # Territorial adjustment
    if (tadjTF){
      if (ncol(Tadj) == length(header.expect.tadj)){
        if (!all(colnames(Tadj) == header.expect.tadj)){
          if (all(toupper(colnames(Tadj)) == toupper(header.expect.tadj))){
            colnames(Tadj) <- header.expect.tadj
            cat("\nThere was a case error in the column names of ", 
              potential.file.names["tadj"], "\nnames reassigned correctly, but you should change this!\n\n", 
              file = log.file, append = TRUE)
          } else {
            stop("\nEither the spelling or the order of the column names is off in ", 
              potential.file.names["tadj"], ".\nCorrect this before continuing. The correct columns should be:\n\n", 
              paste(header.expect.tadj, collapse = ", "),"\n")
          }
        } 
      } else {
       # special case for tadj, since technically it is not documented to have an LDB column:
        if (all(toupper(colnames(Tadj)) == toupper(header.expect.tadj[-length(header.expect.tadj)]))){
          Tadj$LDB <- 1
          colnames(Tadj) <- header.expect.tadj
          cat("\nThere was no LDB column in ", potential.file.names["tadj"], 
            ".\nPreviously, for some reason, the tadj file did not have an LDB column, but it is supposed to now.\nThe column has been added with an assumed value of 1 in all rows.\nPlease add this to the original csv version of the file, as it will not be written out\n\n",
            file = log.file, append = TRUE)
        } else {
          stop("\nThe columns of ", 
            potential.file.names["tadj"], ".\nneed to be fixed\nCorrect this before continuing. The correct columns should be:\n\n", 
            paste(header.expect.tadj, collapse = ", "),"\n")
        }
      } 
    }
