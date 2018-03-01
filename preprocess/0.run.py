import argparse, re

# Read input arguments
parser = argparse.ArgumentParser()

g1 = parser.add_argument_group('Input files')
g1.add_argument('--fastqs', help='List of FASTQ files')
g1.add_argument('--ref_db', help='Reference database (FASTQ)')
g2 = parser.add_argument_group('BWA options')
g2.add_argument('--pct_id', help='Percent identity', type=float, default=90)
g2.add_argument('--minlen', help='Minimum length of mapped reads', type=float, default=40)
g3 = parser.add_argument_group('kpileup options')
g3.add_argument('--bqual', help='Minimum base quality', type=float, default=20)
g3.add_argument('--mqual', help='Minimum mapping quality', type=float, default=0)
g3.add_argument('--depth', help='Minimum read depth', type=float, default=10)

args = parser.parse_args()

# Read input data
prefixes = [re.sub('.fastq', '', line.rstrip()) for line in open(args.fastqs)]
ref_db = args.db

# 1) Run BWA on each FASTQ
for prefix in prefixes:
    print 'bwa mem -a %s %s.fastq > %s.sam' %(ref_db, prefix, prefix)

# 2) Filter SAM files
for prefix in prefixes:
    print 'python 1.filter_sam.py %s.sam %s %s > %s.filter.sam' %(prefix, args.pct_id, args.minlen, prefix)

# 3) Convert to BAM
for prefix in prefixes:
    print 'samtools view -bS -F 4 -o %s.bam %s.filter.sam' %(prefix, prefix)
    print 'samtools sort %s.bam -o %s.sorted' %(prefix, prefix)
    print 'samtools index %s.sorted.bam' %(prefix)

# 4) Make gene and sample files
print 'python 2.make_gene_file.py --fst %s --out gene_file.txt' %(ref_db)
out = open('sample_file.txt', 'w')
for prefix in prefixes:
    out.write('%s\n', prefix)
out.close()

# 5) Run kpileup
for prefix in prefixes:
    print 'perl 3.kpileup.pl %s %s.sorted.bam gene_file.txt %s %s %s > %s.kp.txt' %(prefix, prefix, args.bqual, args.mqual, args.depth, prefix)

# 6) Convert to numpy
print 'python 4.kp2np.py --sfile samples.txt --gfile gene_file.txt --out alignments.cPickle'

