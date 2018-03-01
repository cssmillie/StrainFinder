import argparse, cPickle
import numpy as np


# Input arguments
# ---------------

parser = argparse.ArgumentParser()
parser.add_argument('--data', help='Numpy alignments (.cPickle)')
parser.add_argument('--contigs', help='Map of genomes to contigs (sep="\t")')
parser.add_argument('--samples', help='Sample list (sep="\n")')
parser.add_argument('--trim_bp', help='#bp to remove from beg/end of each gene', type=int, default=0)
parser.add_argument('--u_cov', help='Minimum mean(coverage) per sample', type=float, default=10)
parser.add_argument('--z_cov', help='Maximum z-score(coverage) per alignment site', type=float, default=1.5)
args = parser.parse_args()


# Read data
# ---------

# Numpy alignments
# x[contig] = numpy alignment
# y[genome] = concatenated numpy alignments
x = cPickle.load(open(args.data))
y = {}

# Sample list
# M = [sample1, sample2, ...]
M = np.array([line.rstrip() for line in open(args.samples)])

# Contig map
# cmap[genome] = [contig1, contig2, ...]
cmap = {}
for line in open(args.contigs):
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
    n = sum([(np.shape(data[contig])[1] - 2*args.trim_bp) for contig in contigs])
    k = 4
    y[genome] = np.zeros([m,n,k])
    
    # Add alignment data
    beg = 0
    end = 0
    for contig in contigs:
        end += (np.shape(data[contig])[1] - 2*args.trim_bp)
        y[genome][:, beg:end, :] = data[contig][:, args.trim_bp: -1*args.trim_bp, :]
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
    i = range(x.shape[0])
    j = range(x.shape[1])
    
    # Select polymorphic sites
    if x.shape[0] > 0 and x.shape[1] > 0:
        pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
        x = x[:,pos,:]
        j = j[pos]
    
    # Filter samples by coverage
    if x.shape[0] > 0 and x.shape[1] > 0:
        cov = coverage(x)
        pos = (cov.mean(axis=1) >= args.u_cov)
        x = x[pos, :, :]
        i = i[pos]
    
    # Filter sites by coverage
    if x.shape[0] > 0 and x.shape[1] > 0:
        zcov = z_coverage(x)
        x[abs(cov) > 1.5,:] = 0
    
    # Re-select polymorphic sites
    if x.shape[0] > 0 and x.shape[1] > 0:
        pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
        x = x[:,pos,:]
        j = j[pos]
    
    # Test empty alignment
    if x.shape[0] == 0 or x.shape[1] == 0:
        x = i = j = np.array([])
    
    # Write alignment
    cPickle.dump(x, open('%s.np.cPickle' %(genome), 'w'))
    
    # Write samples
    out = open('%s.samples.txt', 'w')
    for index in i:
        out.write('%s\n' %(index))
    out.close()
    
    # Write alignment sites
    out = open('%s.sites.txt', 'w')
    for index in j:
        out.write('%s\n' %(index))
    out.close()
