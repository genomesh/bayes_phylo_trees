# Bayesian Inference of Phylogenetic Trees

This work was done during a research internship at the STOR-i CDT at Lancaster University, under the supervision of [Georgios Aliatimis](https://www.lancaster.ac.uk/dsi/about-us/members/george-aliatimis). Please see the attached poster or slides for a more detailed description of the project.

### Project Outline

Researching [phylogenetic trees](https://en.wikipedia.org/wiki/Phylogenetic_tree) is important for tasks such as predicting fast-evolving viruses like HIV, and understanding the origin of life through projects such as the Tree of Life. This was first done by Darwin, who used morphological features of birds to infer their evolutionary relationships. With the abundance of DNA sequence data now available, it is possible to infer phylogenetic trees with less uncertainty. A very common approach is to define a [model of DNA evolution](https://en.wikipedia.org/wiki/Models_of_DNA_evolution), then find tree phylogenetic tree that maximises this likelihood - the MLE tree (this can be achieved by using programs such as IQ-TREE [1]). However, what if we want to incorporate prior beliefs about the tree, and understand more about the distribution of different trees than just the mode? One approach is using [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) to sample from this distribution, and then analyse the samples.

### Overview of Code

Apologies for my unclear file names. In more recent projects I have improved this.

| File | Description |
| ---- | ----------- |
| run_analysis.sh | The bulk of the computation. Given a gene, does some preprocessing, finds the MLE using IQ-TREE [1] then runs MCMC using MrBayes [2]. |
| data_to_csv.R | Analyses the MCMC sample trees and MLE tree, and saves the results in a csv. We were particularly interested in the relationship between the proportion of samples with the same topology as the MLE tree, and the length of the smallest branch in the MLE tree. |
| script.sh | Performs run_analysis.sh for each gene. |
| data_to_csv.sh | Performs data_to_csv.R for each gene. |
| davidBayes.R | My own implementation of MCMC for phylogenetic trees. |
| mb_script.nex | Mr Bayes MCMC and model parameters. |
| script.R | Final analyis for all results. |

### References

1. IQ-TREE: B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534. https://doi.org/10.1093/molbev/msaa015

2. MrBayes: Ronquist, F., M. Teslenko, P. van der Mark, D.L. Ayres, A. Darling, S. Höhna, B. Larget, L. Liu, M.A. Suchard, and J.P. Huelsenbeck. 2012. MRBAYES 3.2: Efficient Bayesian phylogenetic inference and model selection across a large model space. Syst. Biol. 61:539-542. https://nbisweden.github.io/MrBayes/index.html