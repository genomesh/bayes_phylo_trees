library(ape)

loadTrees <- function(fileName){
  trees <- read.nexus(file=fileName, tree.names = NULL, force.multi = FALSE)
  startBurn <- ceiling(0.25 * length(trees))
  return(trees[startBurn:length(trees)])
}

partition <- function(trees){
  partitions <- list(1)
  getPartition <- function(tree,id) {
    for (p in 1:length(partitions)) {
      comparison_id <- partitions[[p]][1]
      if (all.equal(tree, trees[[comparison_id]],use.edge.length = F)) { #.phylo ?
        partitions[[p]] <- c(partitions[[p]],id)
        return(partitions)
      }
    }
    partitions[length(partitions)+1] <- id
    return(partitions)
  }
  for (i in 2:length(trees)) { #replace 5 with numTrees
    tree <- trees[[i]]
    partitions <- getPartition(tree, i)
  }
  partitionSizes <- length(partitions[[1]])
  if (length(partitions) != 1) {
    for (i in 2:length(partitions)) {
      partitionSizes <- c(partitionSizes, length(partitions[[i]]))
    }
  }
  orderedPartitionIndices <- order(partitionSizes,decreasing=T)
  largestPartitionIndex <- orderedPartitionIndices[1]
  
  return(list(
    partitions=partitions,
    orderedPartitionIndices=orderedPartitionIndices,
    partitionSizes = sort(partitionSizes,decreasing=T),
    largestPartitionIndex=largestPartitionIndex,
    largestPartition=partitions[[largestPartitionIndex]],
    majority_tree=trees[[partitions[[largestPartitionIndex]][1]]],
    proportion_in_majority = length(partitions[[largestPartitionIndex]]) / length(trees)
  ))
}

# ------------------- Large Gene Data Set Analysis ---------------- #

data_to_csv <- function(gene_id) {
  mcmc_trees <- loadTrees('gene.nex.t')
  mcmc_parts <- partition(mcmc_trees)
  mle_tree <- read.tree('gene.fasta.treefile')
  gene_data <- data.frame(
    gene_id = gene_id,
    prop_in_majority = mcmc_parts$proportion_in_majority,
    min_length = min(mle_tree$edge.length),
    mle_in_majority = all.equal.phylo(mle_tree, mcmc_parts$majority_tree, use.edge.length = F)
  )
  write.table(gene_data, file = '~/dev/stori/housekeeping_analysis.csv',quote=F, row.names=F, col.names=F, append=T,sep=',')
  #write.csv(gene_data, file = '~/dev/stori/housekeeping_analysis.csv',quote=F, row.names=F)
}

gene_id = commandArgs(trailingOnly=TRUE)

main <- function() {
  if (length(gene_id)!=1) {
    stop("One argument must be supplied (input file).fasta", call.=FALSE)
  } else {
    data_to_csv(gene_id)
  }
}

main()
