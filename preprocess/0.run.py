import argparse, re

# Read input arguments
parser = argparse.ArgumentParser()

g1 = parser.add_argument_group('Required arguments')
g1.add_argument('--fastqs', help='List of FASTQ files', required=True)
g1.add_argument('--ref', help='Reference database (FASTA)', required=True)
g1.add_argument('--map', help='Map of genomes to contigs (tab-delimited)', required=True)
g2 = parser.add_argument_group('BWA options')
g2.add_argument('--pct', help='Percent identity', type=float, default=90)
g2.add_argument('--len', help='Minimum length of mapped reads', type=float, default=40)
g3 = parser.add_argument_group('kpileup options')
g3.add_argument('--bqual', help='Minimum base quality', type=float, default=20)
g3.add_argument('--mqual', help='Minimum mapping quality', type=float, default=0)
g3.add_argument('--depth', help='Minimum read depth', type=float, default=10)
g4 = parser.add_argument_group('Filter options')
g4.add_argument('--tlen', help='Number of base pairs to trim from beg/end of alignment', type=int, default=0)
g4.add_argument('--faln', help='Minimum fraction of aligned sites (per sample)', type=float, default=.5)
g4.add_argument('--mcov', help='Minimum mean coverage (per sample)', type=float, default=10)
g4.add_argument('--dcov', help='Remove sites with coverage > [dcov] standard deviations from the mean', type=float, default=1.5)
g4.add_argument('--npos', help='Randomly subsample [npos] alignment sites', type=int, default=None)
args = parser.parse_args()

# Read input data
prefixes = [re.sub('.fastq', '', line.rstrip()) for line in open(args.fastqs)]

# 1) Run BWA on each FASTQ
for prefix in prefixes:
    print 'bwa mem -a %s %s.fastq > %s.sam' %(args.ref, prefix, prefix)

# 2) Filter SAM files
for prefix in prefixes:
    print 'python 1.filter_sam.py %s.sam %s %s > %s.filter.sam' %(prefix, args.pct, args.len, prefix)

# 3) Convert to BAM
for prefix in prefixes:
    print 'samtools view -bS -F 4 -o %s.bam %s.filter.sam' %(prefix, prefix)
    print 'samtools sort %s.bam -o %s.sorted' %(prefix, prefix)
    print 'samtools index %s.sorted.bam' %(prefix)

# 4) Make gene and sample files
print 'python 2.make_gene_file.py --fst %s --out gene_file.txt' %(args.ref)
out = open('samples.txt', 'w')
for prefix in prefixes:
    out.write('%s\n' %(prefix))
out.close()

# 5) Run kpileup
for prefix in prefixes:
    print 'perl 3.kpileup.pl %s %s.sorted.bam gene_file.txt %s %s %s > %s.kp.txt' %(prefix, prefix, args.bqual, args.mqual, args.depth, prefix)

# 6) Convert to numpy
print 'python 4.kp2np.py --samples samples.txt --gene_file gene_file.txt --out all_alignments.cPickle'

# 7) Concatenate and filter numpy alignments
print 'python 5.filter_np.py --aln all_alignments.cPickle --map %s --samples samples.txt --tlen %s --faln %s --mcov %s --dcov %s --npos %s > filter_np.log' %(args.map, args.tlen, args.faln, args.mcov, args.dcov, args.npos)
