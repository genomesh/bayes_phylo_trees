
rm(list=ls())

setwd('~/dev/stori')
#install.packages("ggplot2")
#library('phangorn')
library('ggplot2')
library('seqinr')
library('ape')
library('stringr')

get_nexus_from_fasta <- function(file_in, file_out) {
  data <- read.fasta(file_in)
  write.nexus.data(data, file=file_out)
}

# get_nexus_from_fasta('data/eas-ccl-wg_10014_2_1-dna-trimmed.fasta', 'housekeeping.nex')

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

# -------------- Alkaloid model comparisons -------------- #

analyse_alkaloids <- function() {
  
  alkaloids_gtr_trees <- loadTrees("output/alkaloids_gtr/alkaloids.nex.run1.t")
  alkaloids_gtr_rates_trees <- loadTrees("output/alkaloids_gtr_rates/alkaloids.nex.run1.t")
  alkaloids_jc_trees <- loadTrees("output/alkaloids_jc/alkaloids.nex.run1.t")
  alkaloids_jc_rates_trees <- loadTrees("output/alkaloids_jc_rates/alkaloids.nex.run1.t")
  
  alkaloids_gtr_trees <- loadTrees("output/alkaloids_gtr/alkaloids.nex.run1.t")
  alkaloids_gtr_rates_trees <- loadTrees("output/alkaloids_gtr_rates/alkaloids.nex.run1.t")
  alkaloids_jc_trees <- loadTrees("output/alkaloids_jc/alkaloids.nex.run1.t")
  alkaloids_jc_rates_trees <- loadTrees("output/alkaloids_jc_rates/alkaloids.nex.run1.t")
  
  alkaloids_gtr_parts <- partition(alkaloids_gtr_trees)
  # 1075 t3 id1, 211 id20, 150 id66, 39 id51
  alkaloids_gtr_rates_parts <- partition(alkaloids_gtr_rates_trees)
  # 502 t1 id1, 249 t2 id7, 213 t3 id4, 188 t4 id82
  alkaloids_jc_parts <- partition(alkaloids_jc_trees)
  # 1427 t3 id1, 36 id215, 1, 1
  alkaloids_jc_rates_parts <- partition(alkaloids_jc_rates_trees)
  # 533 t1 id2, 514 t3 id5, 301 id22, 27 id180
  
  all.equal.phylo(
    alkaloids_jc_rates_trees[[22]],
    alkaloids_gtr_rates_trees[[7]],
    use.edge.length = F
  )
  
  plot(alkaloids_jc_rates_trees[[22]], main='JC rates 3rd big', edge.width = 2)
  plot(alkaloids_jc_rates_trees[[2]], main='T1, rates 1st big', edge.width = 2)
  plot(alkaloids_gtr_rates_trees[[4]], main='T3, rates 3rd big', edge.width = 2)
  
  print(all.equal.phylo(alkaloids_jc_partition$majority_tree, alkaloids_trees[[4]],use.edge.length = F))
}

# --------------------- 5 Gene Variance Analysis --------------------- #

five_gene_analysis <- function() {
  alkaloid_genes <- c(
    "eas-ccl-wg_16390-dna-trimmed.fasta",
    "eas-ccl-wg_16472-dna-trimmed.fasta",
    "eas-ccl-wg_16584-dna-trimmed.fasta",
    "eas-ccl-wg_16606-dna-trimmed.fasta",
    "eas-ccl-wg_17201-dna-trimmed.fasta"
  )
  
  housekeeping_genes <- c(
    "eas-ccl-wg_16158-dna-trimmed.fasta",
    "eas-ccl-wg_16165-dna-trimmed.fasta",
    "eas-ccl-wg_16182-dna-trimmed.fasta",
    "eas-ccl-wg_16190-dna-trimmed.fasta",
    "eas-ccl-wg_16495-dna-trimmed.fasta"
  )
  
  #for (i in 1:5) {
  #  get_nexus_from_fasta(alkaloid_genes[i], paste('alkaloids/alkaloid_', i, '.nex', sep=''))
  #}
  
  fetch_info <- function(index, alkaloids = T){
    if (alkaloids) {
      setwd(paste('~/dev/stori/data/variance_test/alkaloids/alkaloid_', index, sep=''))
      mcmc_trees <- loadTrees(paste("alkaloid_",index,".nex.run1.t",sep=''))
      mcmc_parts <- partition(mcmc_trees)
      mle_tree <- read.tree(paste(alkaloid_genes[index],'.treefile', sep=''))
    } else {
      setwd(paste('~/dev/stori/data/variance_test/housekeeping/housekeeping_', index, sep=''))
      mcmc_trees <- loadTrees(paste("housekeeping_",index,".nex.run1.t",sep=''))
      mcmc_parts <- partition(mcmc_trees)
      mle_tree <- read.tree(paste(housekeeping_genes[index],'.treefile', sep=''))
    }
    list(
      mcmc_trees = mcmc_trees,
      mcmc_parts = mcmc_parts,
      mle_tree = mle_tree,
      proportion_in_majority = mcmc_parts$proportion_in_majority,
      min_length_in_mle = min(mle_tree$edge.length)
    )
  }
  
  alkaloid_data <- list(
    fetch_info(1),
    fetch_info(2),
    fetch_info(3),
    fetch_info(4),
    fetch_info(5)
  )
  
  housekeeping_data <- list(
    fetch_info(1, alkaloid = F),
    fetch_info(2, alkaloid = F),
    fetch_info(3, alkaloid = F),
    fetch_info(4, alkaloid = F),
    fetch_info(5, alkaloid = F)
  )
  
  #majority_tree$tip.label <- sub('-','_',majority_tree$tip.label)
  
  lengths <- alkaloid_data[[1]]$min_length_in_mle
  proportions <- alkaloid_data[[1]]$proportion_in_majority
  num_partitions <- length(alkaloid_data[[1]]$mcmc_parts$partitions)
  
  for (i in 2:5) {
    lengths <- c(lengths, alkaloid_data[[i]]$min_length_in_mle)
    proportions <- c(proportions, alkaloid_data[[i]]$proportion_in_majority)
    num_partitions <- c(num_partitions, length(alkaloid_data[[i]]$mcmc_parts$partitions))
  }
  
  alkaloid_df <- data.frame(lengths = lengths, proportions = proportions, num_partitions = num_partitions)
  
  lengths <- housekeeping_data[[1]]$min_length_in_mle
  proportions <- housekeeping_data[[1]]$proportion_in_majority
  num_partitions <- length(housekeeping_data[[1]]$mcmc_parts$partitions)
  
  for (i in 2:5) {
    lengths <- c(lengths, housekeeping_data[[i]]$min_length_in_mle)
    proportions <- c(proportions, housekeeping_data[[i]]$proportion_in_majority)
    num_partitions <- c(num_partitions, length(housekeeping_data[[i]]$mcmc_parts$partitions))
  }
  
  housekeeping_df <- data.frame(lengths = lengths, proportions = proportions, num_partitions = num_partitions)
  
  plot(alkaloid_df_2[c('lengths', 'proportions')], main='small dataset', ylim=c(0,1))
  points(housekeeping_df_2[c('lengths', 'proportions')], pch = 2)
  
  #plot(housekeeping_df[c('lengths', 'num_partitions')], main='alkaloids')
  
  for (i in 1:5) {
    plot(alkaloid_data[[i]]$mle_tree, main=paste('Alkaloid',i,'MLE Tree'), edge.width = 2)
  }
}

# ------------------- Large Gene Data Set Analysis ---------------- #

#rm(list=ls())

setwd('~/dev/stori')
alkaloids_data <- read.csv('alkaloid_analysis.csv')
housekeeping_data <- read.csv('housekeeping_analysis.csv')

alkaloids_data$dataset <- 'alkaloids'
housekeeping_data$dataset <- 'housekeeping'
combined_data <- rbind(alkaloids_data, housekeeping_data)

plot(alkaloids_data[c('min_length', 'prop_in_majority')], main = 'alkaloids', ylim=c(0,1))
plot(housekeeping_data[c('min_length', 'prop_in_majority')], main = 'housekeeping', ylim=c(0,1))

plot(alkaloids_data[c('min_length', 'prop_in_majority')], main = 'both', ylim=c(0,1))
points(housekeeping_data[c('min_length', 'prop_in_majority')], pch = 2)

table(alkaloids_data['mle_in_majority'])
table(housekeeping_data['mle_in_majority'])

plot(alkaloids_data[alkaloids_data$mle_in_majority,][c('min_length', 'prop_in_majority')], main = 'alkaloids where mle in majority', ylim=c(0,1))
plot(alkaloids_data[!alkaloids_data$mle_in_majority,][c('min_length', 'prop_in_majority')], main = 'alkaloids where mle not in majority', ylim=c(0,1))

plot(housekeeping_data[housekeeping_data$mle_in_majority,][c('min_length', 'prop_in_majority')], main = 'housekeeping where mle in majority', ylim=c(0,1))
plot(housekeeping_data[!housekeeping_data$mle_in_majority,][c('min_length', 'prop_in_majority')], main = 'housekeeping where mle not in majority', ylim=c(0,1))

plot(alkaloids_data[alkaloids_data$mle_in_majority,][c('min_length', 'prop_in_majority')], main = 'alkaloids by whether mle in majority', ylim=c(0,1))
points(alkaloids_data[!alkaloids_data$mle_in_majority,][c('min_length', 'prop_in_majority')], pch=2)

ggplot(combined_data, aes(x=min_length, y=prop_in_majority, colour=dataset)) + 
  geom_point(alpha=0.95) +
  ylim(0,1)

# add legend
# conclusions?
# better way to represent data
# analysis of MLE agreeing with majority tree
# analyse alkaloid gene_262 <- 56% in majority but not MLE! min branch length 0.006 (big)
setwd('~/dev/stori/data/new_alkaloids/gene_262')
mcmc_trees <- loadTrees('gene.nex.t')
mcmc_parts <- partition(mcmc_trees)
mle_tree <- read.tree('gene.fasta.treefile')
all.equal.phylo(mle_tree, mcmc_parts$majority_tree, use.edge.length = F)
all.equal.phylo(mle_tree, mcmc_trees[[2]], use.edge.length = F)
plot(mle_tree)
plot(mcmc_parts$majority_tree)


mcmc_trees_2 <- loadTrees("~/dev/stori/data/alkaloid_262_test/gene.nex.run1.t")
mcmc_parts_2 <- partition(mcmc_trees_2)
all.equal.phylo(mle_tree, mcmc_trees_2[[7]], use.edge.length = F)

mle_2 <- read.tree('~/dev/stori/data/alkaloid_262_test/eas-ccl-wg_16214-dna-trimmed.fasta.treefile')

all.equal.phylo(mle_tree, mle_2, use.edge.length = F)

# RESULT: reproducable mcmc giving majority topology which excludes mle
# 56.2% in majority topology
# 38.4% in mle topology


renamed_mle <- mle_tree
renamed_majority <- mcmc_parts$majority_tree


relabel <- function(cur) {
  labels <- mle_tree$tip.label
  new_label <- which(labels == cur)
  return(new_label)
}

for (i in 1:11) {
  cur <- renamed_majority$tip.label[i]
  renamed_majority$tip.label[i] <- relabel(cur)
}

plot(renamed_mle, type = 'tidy')
plot(renamed_majority, type = 'tidy')

plot(housekeeping_data[housekeeping_data$mle_in_majority,][c('min_length', 'prop_in_majority')],
     ylim=c(0,1),
     xlab='Minimum branch length in MLE',
     ylab='Proportion in MLE Topology'
)
# 632, 441, other 53
my_list <- list(majority_topol = 632, mle_topol = 441, other = 53)
plot(my_list)
barplot(c(632, 441, 53) / 1126, names.arg = c('Majority', 'MLE', 'Other'), main = 'Proportion of MCMC samples in different topologies')

x <- seq(-20,20,length.out = 1000)
y <- dnorm(x, mean=-10, sd=1) + 3.5 * dnorm(x, mean=6, sd=4)
plot(x,y, type='l', ylab='Posterior Density', xlab='Branch Length')
abline(v = 0)
