import argparse, cPickle
import numpy as np


# Input arguments
# ---------------

parser = argparse.ArgumentParser()
parser.add_argument('--aln', help='Numpy alignments (.cPickle)')
parser.add_argument('--map', help='Map of genomes to contigs (tab-delimited)')
parser.add_argument('--samples', help='List of samples')
parser.add_argument('--tlen', help='Number of base pairs to trim from beg/end of each contig', type=int, default=0)
parser.add_argument('--faln', help='Minimum fraction of aligned sites (per sample)', type=float, default=.5)
parser.add_argument('--mcov', help='Minimum mean coverage (per sample)', type=float, default=10)
parser.add_argument('--dcov', help='Remove sites with coverage > [dcov] standard deviations from the mean', type=float, default=1.5)
args = parser.parse_args()


# Read data
# ---------

# Numpy alignments
# x[contig] = numpy alignment
# y[genome] = concatenated numpy alignments
x = cPickle.load(open(args.aln))
y = {}

# Sample list
# M = [sample1, sample2, ...]
M = np.array([line.rstrip() for line in open(args.samples)])

# Contig map
# cmap[genome] = [contig1, contig2, ...]
cmap = {}
for line in open(args.map):
    line = line.rstrip().split()
    genome = line[0]
    contig = line[1]
    if genome not in cmap:
        cmap[genome] = []
    cmap[genome].append(contig)


# Concatenate contigs
# -------------------

for genome in cmap:
    contigs = cmap[genome]
    
    # Initialize array
    m = len(M)
    n = sum([(np.shape(x[contig])[1] - 2*args.tlen) for contig in contigs])
    k = 4
    y[genome] = np.zeros([m,n,k])
    
    # Add alignment data
    beg = 0
    end = 0
    for contig in contigs:
        end += (np.shape(x[contig])[1] - 2*args.tlen)
        if args.tlen == 0:
            y[genome][:, beg:end, :] = x[contig]
        else:
            y[genome][:, beg:end, :] = x[contig][:, args.tlen: -1*args.tlen, :]
        beg = end


# Alignment filtering
# -------------------

def coverage(x):
    # calculate coverage for each sample
    # returns MxN numpy array
    return x.sum(axis=2)
    
def z_coverage(x):
    # calculate standardized coverage for each sample
    # returns MxN numpy array
    cov = coverage(x)
    zcov = (cov - cov.mean(axis=1)[:, np.newaxis]) / cov.std(axis=1)[:, np.newaxis]
    return zcov


for genome in y:
    
    # Get alignment data
    x = y[genome]
    i = np.array(range(x.shape[0]))
    j = np.array(range(x.shape[1]))
    print '\n%s [%d samples x %d sites]' %(genome, len(i), len(j))
    
    # Filter samples by fraction of aligned sites
    if x.shape[0] > 0 and x.shape[1] > 0:
        cov = coverage(x)
        pos = ((cov == 0).sum(axis=1) > args.faln*x.shape[1])
        x = x[pos,:,:]
        i = i[pos]
        print '\tFiltering samples by fraction aligned [%d x %d]' %(len(i), len(j))
    
    # Select SNP positions
    if x.shape[0] > 0 and x.shape[1] > 0:
        pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
        x = x[:,pos,:]
        j = j[pos]
        print '\tSelecting polymorphic sites [%d x %d]' %(len(i), len(j))
    
    # Filter samples by mean coverage
    if x.shape[0] > 0 and x.shape[1] > 0:
        cov = coverage(x)
        pos = (cov.mean(axis=1) >= args.mcov)
        x = x[pos, :, :]
        i = i[pos]
        print '\tFiltering samples by mean coverage [%d x %d]' %(len(i), len(j))
    
    # Filter sites by atypical coverage
    if x.shape[0] > 0 and x.shape[1] > 0:
        zcov = z_coverage(x)
        x[abs(zcov) > args.dcov,:] = 0
        print '\tZeroing %d sites with atypical coverage' %((abs(zcov) > args.dcov).sum())
    
    # Select SNP positions
    if x.shape[0] > 0 and x.shape[1] > 0:
        pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
        x = x[:,pos,:]
        j = j[pos]
        print '\tSelecting polymorphic sites [%d x %d]' %(len(i), len(j))
    
    # Test empty alignment
    if x.shape[0] == 0 or x.shape[1] == 0:
        print'\tSkipping genome [%d samples x %d sites]' %(len(i), len(j))
        continue
    
    # Write alignment
    print '\tWriting files: %s.np.cPickle, %s.samples.txt, %s.sites.txt' %(genome, genome, genome)
    cPickle.dump(x, open('%s.np.cPickle' %(genome), 'w'))
    
    # Write samples
    out = open('%s.samples.txt' %(genome), 'w')
    for index in i:
        out.write('%s\n' %(M[index]))
    out.close()
    
    # Write alignment sites
    out = open('%s.sites.txt' %(genome), 'w')
    for index in j:
        out.write('%s\n' %(index))
    out.close()
