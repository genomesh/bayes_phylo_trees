library('seqinr')
library('ape')

create_nexus_from_fasta <- function(file_name) {
  data <- read.fasta(file_name)
  write.nexus.data(data, 'gene.nex')
}

file_name = commandArgs(trailingOnly=TRUE)

main <- function() {
  if (length(file_name)!=1) {
    stop("One argument must be supplied (input file).fasta", call.=FALSE)
  } else {
    create_nexus_from_fasta(file_name)
  }
}

main()