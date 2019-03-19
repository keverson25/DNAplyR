#' Read a clustal (.aln) alignment file
#'
#' This function reads a clustal (.aln) alignment file and 
#'      converts it to an easy-to-read data frame object.
#' @param file The name of the file to be read. The file name 
#'      should be provided within quotes, either as an absolute 
#'      file path, or a relative file path given the current 
#'      working directory, getwd().
#' @return Returns a data frame with taxon names as row names 
#'      and a column for each site in the alignment.
#' @keywords read clustal aln DNA alignment
#' @export
#' @examples
#' read.clustal.alignment()

read.clustal.alignment <- function (file)
{
  # 
  # Check that we have read permission on the file:
  # 
  file <- path.expand(file) 
  if(file.access(file, mode = 4) != 0) stop(paste("R does not have permission to read", file))
  
  #
  # Start reading the file
  #
  X <- scan(file = "Primates.aln", skip=1, what = character(),
            quiet = TRUE, comment.char = "[", strip.white = TRUE)
  #
  # Remove the lines that display conserved sites
  #
  X <- X[!grepl("\\.", X)]
  X <- X[!grepl("\\:", X)]
  X <- X[!grepl("\\*", X)]
  #
  # Prep the data
  #
  X<-matrix(X, ncol=2, byrow=T)
  taxaSort <- unique(X[,1])
  X<-X[order(X[,1]),]
  taxa<-unique(X[,1])
  #
  # Generate the data frame
  #
  ObjMat<-matrix(unlist(strsplit(X[,2],split="")), nrow=length(taxa), ncol=length(strsplit(X[grep(taxa[i], X)+1], byrow=T)
  rownames(ObjMat)<-taxa
  ObjMat<-ObjMat[order(match(rownames(ObjMat),taxaSort)),]
}