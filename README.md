# Strain Finder
If you have any questions about Strain Finder, feel free to contact me (csmillie AT mit DOT edu).

Strain Finder 1.0 was written by Jonathan Friedman and can be found here:
https://bitbucket.org/yonatanf/strainfinder/

## Overview
Strain Finder takes as input a reference genome alignment and the number of strains to infer. It finds maximum likelihood estimates for the strain genotypes and the strain frequencies across all metagenomes.

Strain Finder uses the EM algorithm to perform the optimization. Because EM only converges to a local optimum, but not necessarily a global optimum, you should run Strain Finder with many initial conditions and select the estimate with the best likelihood. Additionally, because the number of strains is not known in advance, you should run it for 2-N strains. You can select the optimal number of strains with model selection criteria, such as AIC, BIC, or the LRT.

## Quick start
The input to Strain Finder is a cPickled numpy alignment. To generate this file, map your metagenomes against a reference, then use a variant caller (such as mpileup) to count the SNPs at every position. After you have generated this file (see details below), you are ready to use Strain Finder. The easiest way to run Strain Finder:

python StrainFinder.py -N 5 --aln aln.cpickle --em em.cpickle --max\_reps 10 --n\_keep 3 --dtol 1 --ntol 2 --max\_iter 25 --max\_time 3600 --converge --em\_out em.cpickle --force\_update --merge\_output --otu\_out otu\_table.txt --log log.txt 

This command reads the alignment data from aln.cpickle (or em.cpickle if it exists). It estimates strains from --max\_reps 10 initial conditions, keeping the 3 best estimates. Each search terminates after the local convergence criteria (specified by --dtol, --ntol, and --max\_iter) have been met, or if the --max\_time limit of 3600 seconds has been reached. New searches started with the --converge command will pick up where the last search left off. It saves the results in em.cPickle and writes the strain profiles to otu\_table.txt. For parallelization, you can submit many identical jobs and they will optimize different estimates, communicating via log.txt. Usually, I send a lot of these jobs to the cluster overnight, wake up, and take the day off knowing that Strain Finder worked so hard overnight.

## Inputs
• Numpy array (--aln)

To start, the input to Strain Finder is a cPickle numpy array with dimensions (M, L, 4), with:

M = number of samples
L = number of alignment positions (alignment length)
4 = nucleotides (A,C,G,T)

This array is your alignment. Each entry (i,j,k) represents how often you observe nucleotide k at position j in sample i of your alignment.

• EM file (--em)

Alternatively, if you have run Strain Finder and saved the results as an EM object, you can input this file directly using the '--em' option. This lets you refine estimates in your EM object or explore more initial conditions.

• Simulated data (--sim)

You can also simulate alignments with the options in the 'Simulation' group. Strain genotypes can be simulated on a phylogenetic tree with the --phylo and -u options. You can also add noise to your alignment with the --noise option.

## Search strategies
Strain Finder starts with an initial guess for the strain genotypes. This guess is informed by the majority SNPs in each sample. Strain Finder then uses the EM algorithm to iteratively optimize the strain frequencies and genotypes until they converge onto a final estimate. The --max\_reps option specify the minimum and maximum numbers of initial conditions to explore. After an EM object has exceeded --max\_reps, it will no longer perform additional searches. 

At any given time, Strain Finder will only hold a fixed number of estimates. You can control this number with the --n\_keep option. For example, --n\_keep 3 tells Strain Finder to save the 3 estimates with the best log-likelihoods. If it finds a better estimate, it will automatically replace the estimate with the lowest log-likelihood.

## Local convergence
At a certain point, Strain Finder will stop making significant gains in likelihood. When this happens, it is better to search more initial conditions than to refine your current estimate. After each iteration, Strain Finder measures the log-likelihood increase. If this increase is less than --dtol for --ntol iterations, then the current estimate has converged.

## Global convergence
You can also specify global convergence criteria (i.e. convergence between estimates). For example, suppose that after searching 100 initial conditions, your best estimates have all converged on a common solution. The degree to which these estimates must converge can be specified with the --min\_fdist and --min\_gdist options.

## Parallelization
Strain Finder uses the --log file to support parallelization. By reading and writing to this file, multiple processes can communicate the results of their optimizations with each other. It is much faster to read this file than it is to load an EM object, only to discover it has already been optimized.

The --max\_time option interrupts a search if it has hit the time limit and saves the results before exiting. This is useful if your cluster imposes time limits on submitted jobs.

## Output files
• EM file (--em_out)

An EM object is a binary cPickled file. This object holds: (1) the input alignment, (2) simulated data (a Data object), (3) the strain genotypes, and the strain frequencies (Estimate objects).

To load the EM object:
from StrainFinder import \*
import cPickle
em = cPickle.load(open(em_file, 'rb'))

To access the alignment data:
em.data # data object
em.data.x # alignment data, dim = (M x L x 4)

To access the estimates:
em.estimates # list of estimate objects
em.select_best_estimates(1) # estimate with best likelihood

To access the strain genotypes
em.estimates[0].p # genotypes of first estimate, dim = (N x L x 4)
em.estimates[0].p.get_genotypes() # fasta format

To access the strain frequencies
em.estimates[0].z # frequencies of first estimate, dim = (M x N)

• Alignment (--aln_out)

If you simulated an alignment, you can save it as a cPickled numpy array using this option.

• Data object (--data_out)

If you simulated an alignment, you can save the alignment, along with the underlying strain genotypes and strain frequencies, using this option.

• OTU table (--otu\_out)
This writes the strain genotypes and strain frequencies as an OTU table. The strain genotypes are included in the OTU names.

## Extras
Strain Finder also has options for robust estimation (automatically ignore incompatible alignment sites) and to exhaustively search strain genotype space (instead of numerical optimization).

## Notes
Do not use the shallow and deep search options. Instead, use the search strategy outlined above (--max_reps and local convergence)

The global convergence criteria do not work with insufficient data. Future plans to mask sites with low coverage when calculating the genotype distances among estimates.

In general, always use --merge\_out and --force\_update
