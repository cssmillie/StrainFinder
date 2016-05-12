# StrainFinder
StrainFinder

# Overview

Strain Finder takes as input a reference genome alignment and the number of strains to infer. It finds maximum likelihood estimates for the strain genotypes and the strain frequencies in each metagenomics sample.

The entry point to Strain Finder is the 'StrainFinder.py' script. This uses the EM algorithm to perform the optimization. Because the EM algorithm converges to a local optimum, but not necessarily a global optimum, you need to run Strain Finder with many initial conditions and select the estimate with the best likelihood score. In addition, because the number of strains is usually not known in advance, you need to run Strain Finder for N=2,3,...,k strains. You can find the optimal number of strains with model selection criteria, such as the AIC/BIC scores or the LRT.

# Inputs

1. Numpy array (--aln)

Initially, the input to Strain Finder is a cPickle numpy array with dimensions: (M, L, 4). You specify the name of this array with the '--aln' flag.

M = number of samples
L = number of alignment positions (alignment length)
4 = nucleotides (A,C,G,T)

This array is your alignment. Each entry (i, j, k) represents how often you observe nucleotide k at position j in sample i.

2. EM file (--em)

If you have previously run Strain Finder and saved the results as an EM file, you can input this file directly into Strain Finder using the '--em' flag. You can refine estimates that are in your EM file or search more initial conditions.

3. Simulated data (--sim)

You can simulate alignments with the options in the 'Simulation' group. Alignments can be simulated on a phylogenetic tree using the --phylo and -u options. You can also disrupt alignments by randomly perturbing sites with the --noise option.

# Search strategies

Strain Finder first generates random guesses for the strain genotypes and the strain frequencies. These are the initial conditions. It then uses the EM algorithm to optimize these guesses until they converge onto a final estimate. Because the EM algorithm can take a long time to converge, it may be effective to divide your searches into 'shallow' searches and 'deep' searches. Shallow searches are meant to quickly explore many initial conditions, while deep searches are meant to refine these estimates.

For example, you might use shallow searches to quickly explore 1000 starting conditions, each time running 3 iterations of the EM algorithm. You can then refine the estimate with the best likelihood using a deep search with 100 iterations of the EM algorithm.

At any given time, Strain Finder will only hold a fixed number of your total estimates. You can control this number with the '--n_keep' flag. For example, '--n_keep 3' tells Strain Finder to save the 3 estimates with the best log-likelihoods. If it finds a new estimate with a better likelihood, it will replace the estimate that has the lowest likelihood.

You can also run each search until it converges with the '--converge' flag. This runs a deep search on every initial condition.

