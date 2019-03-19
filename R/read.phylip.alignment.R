#' Read a Phylip alignment file
#'
#' This function reads a phylip alignment file and converts it
#'      to an easy-to-read data frame object.
#' @param file The name of the file to be read. The file name 
#'      should be provided within quotes, either as an absolute 
#'      file path, or a relative file path given the current 
#'      working directory, getwd().
#' @return Returns a data frame with taxon names as row names 
#'      and a column for each site in the alignment.
#' @keywords read phylip DNA alignment
#' @export
#' @examples
#' read.phylip.alignment()

read.phylip.alignment <- function(file) {
  spp.bps <- scan(file = file, n = 2, quiet = TRUE)
  Obj <- scan(file = file, what = "", skip = 1, quiet = TRUE)
  ObjMat <- matrix(nrow = spp.bps[[1]], ncol = spp.bps[[2]])
  rownames(ObjMat) <- Obj[seq(from=1, to=length(Obj), by=2)]
  for(i in seq(from=2, to=length(Obj), by=2)){
    ObjMat[(i/2),1:spp.bps[[2]]]<-strsplit(Obj[i], split="")[[1]] 
  }
  ObjMat<-as.data.frame(ObjMat)
  ObjMat
}