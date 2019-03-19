#' Read a Nexus alignment file
#'
#' This function reads a nexus alignment file and converts it
#'      to an easy-to-read data frame object.
#' @param file The name of the file to be read. The file name 
#'      should be provided within quotes, either as an absolute 
#'      file path, or a relative file path given the current 
#'      working directory, getwd().
#' @return Returns a data frame with taxon names as row names 
#'      and a column for each site in the alignment.
#' @keywords read nexus DNA alignment
#' @export
#' @examples
#' read.nexus.alignment()

read.nexus.alignment <- function (file)
{
  # Simplified NEXUS data parser.
  #
  # Version: 03/15/2019 12:23 PM EST
  #
  # By:      Johan Nylander, nylander @ scs.fsu.edu
  # Edited for DNAplyR by: Kathryn M. Everson, kathryn.everson @ uky.edu
  #
  # WARNING: This is parser reads a restricted nexus format,
  #          see README for details.
  #
  # Argument (x) is a nexus formatted data file.
  #
  #------------------------------------------------------------------
  
  "find.ntax" <- function (x)
  {
    for (i in 1:NROW(x)) {
      if(any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
        ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                               "\\3", x[i], perl = TRUE, ignore.case = TRUE))
        break
      }
    }
    ntax
  }
  
  "find.nchar" <- function (x)
  {
    for (i in 1:NROW(x)) {
      if(any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
        nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
        break
      }
    }
    nchar
  }
  
  "find.matrix.line" <- function (x)
  {
    for (i in 1:NROW(x)) {
      if(any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
        matrix.line <- as.numeric(i)
        break
      }
    }
    matrix.line
  }
  
  "trim.whitespace" <- function (x)
  { 
    gsub("\\s+", "", x)
  }
  
  "trim.semicolon" <- function (x)
  {
    gsub(";", "", x)
  }
  
  #
  # Check that we have read permission on the file:
  # 
  file <- path.expand(file) 
  if(file.access(file, mode = 4) != 0) stop(paste("R does not have permission to read", file))
  
  #
  # Start reading the file
  #
  X <- scan(file = file, what = character(), sep = "\n",
            quiet = TRUE, comment.char = "[", strip.white = TRUE)
  ntax <- find.ntax(X)
  nchar <- find.nchar(X)
  matrix.line <- find.matrix.line(X)
  start.reading <- matrix.line + 1
  Obj <- list()
  length(Obj) <- ntax
  i <- 1
  pos <- 0
  tot.nchar <- 0
  tot.ntax <- 0
  
  for (j in start.reading:NROW(X)) {
    Xj <- trim.semicolon(X[j])
    if(Xj == "") {
      break
    }
    if(any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE))) {
      break
    }
    ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
    if (length(ts) > 2) {
      stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
    }
    if (length(ts) !=2) {
      stop("nexus parser failed to read the sequences (ts!=2)")
    }
    Seq <- trim.whitespace(ts[2])
    Name <- trim.whitespace(ts[1])
    nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
    if (any(l <- grep(nAME, names(Obj)))) {
      tsp <- strsplit(Seq, NULL)[[1]]
      for (k in 1:length(tsp)) {
        p <- k + pos
        Obj[[l]][p] <- tsp[k]
        chars.done <- k
      }
    }
    else {
      names(Obj)[i] <- Name
      tsp <- strsplit(Seq, NULL)[[1]]
      for (k in 1:length(tsp)) {
        p <- k + pos
        Obj[[i]][p] <- tsp[k]
        chars.done <- k
      }
    }
    tot.ntax <- tot.ntax + 1
    if (tot.ntax == ntax) {
      i <- 1
      tot.ntax <- 0 
      tot.nchar <- tot.nchar + chars.done
      if (tot.nchar == nchar*ntax) {
        print("ntot was more than nchar*ntax")
        break
      }
      pos <- tot.nchar
    }
    else {
      i <- i + 1
    }
  }
  if (tot.ntax != 0) {
    cat("ntax:",ntax,"differ from actual number of taxa in file?\n")
    stop("nexus parser did not read names correctly (tot.ntax!=0)")
  }
  for (i in 1:length(Obj)) {
    if (length(Obj[[i]]) != nchar) {
      cat(names(Obj[i]),"has",length(Obj[[i]]),"characters\n")
      stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
    }
  }
  Obj <- lapply(Obj, tolower)
  ObjMat<-matrix(unlist(Obj), ncol = length(Obj[[1]]), byrow = TRUE)
  rownames(ObjMat)<-names(Obj)
  ObjMat<-as.data.frame(ObjMat)
  ObjMat
}