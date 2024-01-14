rm(list=ls())

setwd('~/dev/stori')

#install.packages('gtools')

library('seqinr')
library('ape')
library('gtools')

custom_trees <- read.nexus(file="output/custom_jc/custom.nex.run1.t", tree.names = NULL, force.multi = FALSE)
custom_tree <- custom_trees[[1]]
custom_data <- read.nexus.data(file="data/custom.nex")

final_tree <- custom_trees[[101]]

plot(custom_trees[[101]], main='custom tree', edge.width = 2)

#primates_data <- read.nexus.data(file="data/primates.nex")

get_likelihood <- function(tree, gene_data, model=jc_model){
  taxa_names <- names(gene_data)
  num_loci <- length(gene_data[[taxa_names[1]]])
  likelihood <- 1
  num_nodes <- tree$Nnode
  num_taxa <- num_nodes + 2
  nucs <- c('a', 'g', 'c', 't')
  node_nuc_perms <- permutations(4, num_nodes, v=c('a','g','c','t'), repeats.allowed = T)
  for (loci in 1:num_loci) {
    loci_nucs <- gene_data[[taxa_names[1]]][loci]
    for (i in 2:num_taxa) {
      loci_nucs[i] <- gene_data[[taxa_names[i]]][loci]
    }
    sum <- 0
    for (perm in 1:(4^num_nodes)) {
      loci_nucs[num_taxa + 1:num_nodes] <- node_nuc_perms[perm,]
      #print(loci_nucs)
      perm_likelihood <- model(tree$edge, tree$edge.length, loci_nucs)
      #print(perm_likelihood)
      sum <- sum + perm_likelihood
    }
    likelihood <- likelihood * sum
  }
  likelihood
}

jc_model <- function(edges, edge.lengths, nucs, mu = 1) {
  # mu is overall substitution rate.
  prod <- 1
  num_edges <- length(edges[,1])
  for (e in 1:num_edges) {
    if (nucs[edges[e,1]] == nucs[edges[e,2]]) {
      prod <- prod * (0.25 + 0.75 * exp(-mu*edge.lengths[e]))
    } else {
      prod <- prod * (0.25 - 0.25 * exp(-mu*edge.lengths[e]))
    }
  }
  prod
}



log(get_likelihood(custom_tree, custom_data))

propose_tree <- function(cur_tree, cur_likelihood, gene_data) {
  num_edges <- length(cur_tree$edge.length)
  edge_to_change <- sample(1:num_edges, 1)
  current_length <- cur_tree$edge.length[edge_to_change]
  unif <- runif(1, min = -1, max = 1)
  lambda <- 1/4
  proposed_length <- exp(log(current_length) + lambda * unif)
  proposed_tree <- cur_tree
  proposed_tree$edge.length[edge_to_change] <- proposed_length
  proposed_likelihood <- get_likelihood(proposed_tree, gene_data)
  prob_accept <- min(1, proposed_likelihood / cur_likelihood)
  list(tree = proposed_tree,
       prob_accept = prob_accept,
       likelihood = proposed_likelihood
  )
}

run_mcmc <- function(start_tree, gene_data, num_iterations = 10) {
  cur_tree <- start_tree
  all_trees <- cur_tree
  cur_likelihood <- get_likelihood(start_tree, gene_data)
  all_likelihoods <- cur_likelihood
  for (n in 1:num_iterations) {
    proposal <- propose_tree(cur_tree, cur_likelihood, gene_data)
    unif <- runif(1, min = 0, max = 1)
    if (unif < proposal$prob_accept) {
      cur_tree <- proposal$tree
      cur_likelihood <- proposal$likelihood
    }
    all_trees <- c(all_trees, cur_tree)
    all_likelihoods <- c(all_likelihoods, cur_likelihood)
  }
  list(trees = all_trees, likelihoods = all_likelihoods)
}

# ----------------- MCMC Analysis ----------------- #

mcmc_run_1 <- run_mcmc(custom_tree, custom_data, num_iterations = 50000)

mcmc_run_2 <- run_mcmc(final_tree, custom_data, num_iterations = 10000)

which.max(mcmc_run_2$likelihoods)
max(mcmc_run_2$likelihoods)

plot(log(mcmc_run_2$likelihoods), main='start=final_tree, lambda = 1/4, log likelihoods', type = 'l')


plot(log(mcmc_run_1$likelihoods), main='start=final_tree, lambda = 1/4, log likelihoods', type = 'l')


for (i in c(1,100,1550,3000,10000)) {
  plot(mcmc_run_2$trees[[i]], main = paste('mcmc tree', i) , edge.width = 2)
}
get_likelihood(mcmc_run_1$trees[[35432]], custom_data)

# ----------------- Branch Length Analysis ----------------- #

analyse_branches <- function(tree, gene_data, edge_to_change = 1) {
  trees <- tree
  likelihoods <- get_likelihood(tree, gene_data)
  print(likelihoods)
  for (a in seq(from = -1, to = 1, length.out = 50)) {
    lambda <- 1/4
    current_length <- tree$edge.length[edge_to_change]
    new_length <- exp(log(current_length) + lambda * a)
    new_tree <- tree
    new_tree$edge.length[edge_to_change] <- new_length
    trees <- c(trees, new_tree)
    likelihoods <- c(likelihoods, get_likelihood(new_tree, gene_data))
  }
  list(trees = trees, likelihoods = likelihoods)
}

custom_tree_analysis <- analyse_branches(mcmc_run_1$trees[[35432]], custom_data, edge_to_change = 7)

plot(log(custom_tree_analysis$likelihoods))
