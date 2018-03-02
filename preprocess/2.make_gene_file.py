import argparse

# Make kpileup "gene file" from FASTA reference
# Output columns:
# 1 = Contig ID
# 2 = Gene ID
# 3 = Beg position (1-indexed)
# 4 = End position (1-indexed)
# 5 = Strand (+ or -)
# 6 = Sequence

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='Input FASTA file')
parser.add_argument('--out', help='Output gene file')
args = parser.parse_args()

def iter_fst(fn):
    sid = ''
    seq = ''
    for line in open(fn):
        line = line.rstrip()
        if line.startswith('>'):
            if seq != '':
                yield [sid, seq]
            sid = line
            seq = ''
        else:
            seq += line
    yield [sid, seq]

# Make gene file
out = open(args.out, 'w')
for [sid, seq] in iter_fst(args.fst):
    sid = sid.rstrip().split()
    contig_id = sid[0][1:]
    gene_id = sid[0][1:]
    strand = '+'
    beg = 1
    end = len(seq)
    out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(contig_id, gene_id, beg, end, strand, seq))
out.close()

