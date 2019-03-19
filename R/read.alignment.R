read.alignment <- function(file, format)
{
  #
  # Check that we have read permission on the file:
  # 
  file <- path.expand(file) 
  if(file.access(file, mode = 4) != 0) stop(paste("R does not have permission to read", file))
  #
  # call the correct function to actually read the file
  #
  alignment <- switch( tolower(format),
                 fasta = .Call("read_fasta", file, PACKAGE = "seqinr"), 
                 mase = .Call("read_mase", file, PACKAGE = "seqinr"),
                 phylip = .Call("read_phylip_align", file, PACKAGE = "seqinr"),
                 msf = .Call("read_msf_align", file, PACKAGE = "seqinr"),
                 clustal = .Call("read_clustal_align", file, PACKAGE = "seqinr"),
                 stop("Wrong format name: Available formats are fasta, mase, phylip, msf, clustal")
  )
}
  