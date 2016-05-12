# Strain Finder

If you have any questions about Strain Finder, feel free to contact me (csmillie AT mit DOT edu).

Strain Finder 1.0 was written by Jonathan Friedman and can be found here:
https://bitbucket.org/yonatanf/strainfinder/

## Overview

Strain Finder takes as input a reference genome alignment and the number of strains to infer. It finds maximum likelihood estimates for the strain genotypes and the strain frequencies in each metagenomics sample.

The entry point to Strain Finder is the 'StrainFinder.py' script. This uses the EM algorithm to perform the optimization. Because the EM algorithm converges to a local optimum, but not necessarily a global optimum, you need to run Strain Finder with many initial conditions and select the estimate with the best likelihood score. In addition, because the number of strains is usually not known in advance, you need to run Strain Finder for N=2,3,...,k strains. You can find the optimal number of strains with model selection criteria, such as the AIC/BIC scores or the LRT.

## Inputs

• Numpy array (--aln)

Initially, the input to Strain Finder is a cPickle numpy array with dimensions: (M, L, 4). You specify the name of this array with the '--aln' flag.

M = number of samples
L = number of alignment positions (alignment length)
4 = nucleotides (A,C,G,T)

This array is your alignment. Each entry (i, j, k) represents how often you observe nucleotide k at position j in sample i.

• EM file (--em)

If you have previously run Strain Finder and saved the results as an EM file, you can input this file directly into Strain Finder using the '--em' flag. You can refine estimates that are in your EM file or search more initial conditions.

• Simulated data (--sim)

You can simulate alignments with the options in the 'Simulation' group. Alignments can be simulated on a phylogenetic tree using the --phylo and -u options. You can also disrupt alignments by randomly perturbing sites with the --noise option.

## Search strategies

Strain Finder first generates random guesses for the strain genotypes and the strain frequencies. These are the initial conditions. It then uses the EM algorithm to optimize these guesses until they converge onto a final estimate. Because the EM algorithm can take a long time to converge, it may be effective to divide your searches into 'shallow' searches and 'deep' searches. Shallow searches are meant to quickly explore many initial conditions, while deep searches are meant to refine these estimates.

For example, you might use shallow searches to quickly explore 1000 starting conditions, each time running 3 iterations of the EM algorithm. You can then refine the estimate with the best likelihood using a deep search with 100 iterations of the EM algorithm.

At any given time, Strain Finder will only hold a fixed number of your total estimates. You can control this number with the '--n\_keep' flag. For example, '--n\_keep 3' tells Strain Finder to save the 3 estimates with the best log-likelihoods. If it finds a new estimate with a better likelihood, it will replace the estimate that has the lowest likelihood.

You can also run each search until it converges with the '--converge' flag. This runs a deep search on every initial condition.

## Local convergence

After many iterations, the EM algorithm will hit a point where it stops making significant gains in likelihood. When this happens, it is better to search more initial conditions than to refine your current estimate. After each iteration, Strain Finder measures the log-likelihood increase. If this increase is less than '--dtol' for '--ntol' iterations, then the current estimate has converged.

## Global convergence

In addition to local convergence, you can also specify global convergence criteria (i.e. convergence between estimates). For example, suppose that after searching through 10,000 initial conditions, your best estimates have all converged on a common solution. The degree to which these estimates must converge can be specified with the '--min\_fdist', and '--min\_gdist' options.

## Parallelization

Strain Finder uses the '--log' file to support parallelization. By reading and writing to this file, multiple processes can communicate the results of their optimizations with each other. It is much faster to read this file than it is to load an existing EM object, only to discover it has already been optimized.

A few options can assist with parallelization. The '--min\_reps' and '--max\_reps' options specify the minimum and maximum numbers of initial conditions to explore. After an EM object has exceeded '--max\_reps', it will no longer perform additional searches. The '--max\_time' option will interrupt an existing search if it has hit the time limit and save the results before exiting. This is useful if your cluster has time limits on the jobs you submit.

## Output files
• EM file (--em_out)

An EM object is a binary cPickle file. This object holds: (1) the input alignment and simulated data (a Data object), and (2) the strain genotypes and the strain frequencies (Estimate objects).

To load the EM object:
import cPickle
em = cPickle.load(open(em_file, 'rb'))

To access the alignment data:
em.data # data object
em.data.x # numpy array alignment

To access the estimates:
em.estimates # list of estimate objects
em.select_best_estimates(1) # select best estimate object

To access the strain genotypes
em.estimates[0].p # strain genotype of first estimate

To access the strain frequencies
em.estimates[0].z # strain frequencies of first estimate

• Alignment (--aln_out)

If you simulated an alignment, you can save it as a cPickled numpy array using this option.

• Data object (--data_out)

If you simulated an alignment, you can save the alignment, along with the underlying strain genotypes and strain frequencies, using this option.

• OTU table (--otu\_out)
This writes the strain genotypes and strain frequencies as an OTU table. The strain genotypes are included in the OTU names.

## Other options
There are also some options for robust estimation (automatically ignore incompatible alignment sites) and to exhaustively search genotypes (instead of numerically optimize them).

In general, always use --merge\_out and --force\_update


## How to use it
Once you have generated a cPickled numpy alignment, you are ready to use Strain Finder. If you are only estimating strains for a single reference genome, you can do something like this:

python em.py -N 5 --aln test.aln --s\_reps 100 --s\_iter 3 --d\_reps 1 --d\_iter 100 --n\_keep 3 --dtol 1 --ntol 2 --em\_out test.em --log test.log --force\_update --merge\_output --otu\_out otu\_table.txt

This will estimate the genotypes and frequencies of 5 strains to best explain your alignment. The search strategy is to quickly search 100 initializations for 3 iterations of the EM algorithm each. Then, it selects the 3 best estimates and refines them with deep searches for 100 iterations each. * You can always add new estimates and refine existing estimates later, using the '--em' option *
