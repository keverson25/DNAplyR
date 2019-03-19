#' Read a fasta alignment file
#'
#' This function reads a fasta alignment file and converts it
#'      to an easy-to-read data frame object. Each individual 
#'      must have the same sequence length. The file cannot
#'      contain any comments.
#' @param file The name of the file to be read. The file name 
#'      should be provided within quotes, either as an absolute 
#'      file path, or a relative file path given the current 
#'      working directory, getwd().
#' @return Returns a data frame with taxon names as row names 
#'      and a column for each site in the alignment.
#' @keywords read fasta DNA alignment
#' @export
#' @examples
#' read.fasta.alignment()

read.fasta.alignment <- function(file) {
  lines <- readLines(file)
  ind <- which(substr(lines, 1L, 1L) == ">")
  nseq <- length(ind)
  if(nseq == 0){
    stop("no line starting with a > character found")
  }
  #
  # Localize sequence data:
  #
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  #
  # Read sequences:
  #
  sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))
  #
  # Read sequence names:
  #
  nomseq <- lapply(seq_len(nseq), function(i){
    firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })
  ObjMat <- matrix(nrow = length(nomseq), ncol = length(strsplit(sequences[[1]], split="")[[1]]))
  rownames(ObjMat) <- nomseq
  for(i in 1:nrow(ObjMat)){
    ObjMat[i,]<-strsplit(sequences[[i]], split="")[[1]] 
  }
  ObjMat<-as.data.frame(ObjMat)
  ObjMat
}