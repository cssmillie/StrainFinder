use strict;

my $usage = qq{
    $0 <sampleId> <bam file> <gene file> <minimum BQ> <minimum MQ> <minimum D> 
    
    <sampleId>: sample ID
    <bam file>: an aligned bam file
    <gene file>: a tab-delimited file of six columns in this order: contigId, geneId, begin, end, strand, DNA transcript seq. (Note: begin<end)
    <minimum BQ>: the minimum base quality score of a sequenced base
    <minimum MQ>: the minimum MQ mapping score of the aligned reads
    <minimum D>: the minimum depth, an integer.
   
    This script is for parsing the BAM file and look for reads overlapping with the target genes and report the pileup.
    Reads must pass the minimum MQ threshold.
    For each position, an allele must have a depth >= <minimum D>.
    The frequency of each allele is calculated after the filterings.
    
    Note: the quality score must be Sanger Phred score (ascii 33).
    
    Dependencies:
    Use Samtools

    CIGAR parser: the embedded cigar parser is not versatile, please review it to make sure that it is handling the cigar code appropriately.


    The script generates a tab-delimited table directly to STDOUT.
    Each row is a base of the queried genes.
    A base could be reported repeatedly up to four times as four rows when polymorphisms are observed: A,T, G or C.
    The columns are:
    sample ID
    contig ID
    This Base Position
    Gene ID
    Ref allele of this Base
    Condon Position (if coding gene)
    Observed Consensus Allele of this Base (in this BAM)
    Observed Allele
    Coverage Depth of the Observed Allele
    Allele Frequency of the Observed Allele
    
};

die $usage unless scalar(@ARGV) == 6;
sub parseGene {
    # parse the input file
    my $file = shift;
    my %data = ();
    open(FILE, "$file") or die "Can't read $file:$!\n";
    while (<FILE>) {
        s/\n//;
        my @tabs = split (/\t/);
        my $gene = $tabs[1];
        %{$data{$gene}} = (
                'contig' => $tabs[0],
                'begin' => $tabs[2],
                'end' => $tabs[3],
                'strand' => $tabs[4],
                'seq' => $tabs[5]
                           
        );
 
    }
    close FILE;
    return %data;
    
}

sub parseBases {
    # go through each amphora gene, and add the nucleotide positions to %nuc
    my ($genesR) = @_;
    my %nuc = (); # indexed by contigs
    for my $g (keys %$genesR) {
        #print STDERR "Storing positions of $g:\n";
        my $begin = $genesR->{$g}{'begin'};
        my $end = $genesR->{$g}{'end'};
        my $c = $genesR->{$g}{'contig'};
        my $strand = $genesR->{$g}{'strand'};
        my @temp = split //, $genesR->{$g}{'seq'};
        for (my $i = $begin; $i <= $end; $i++){
            my $codonPos = ($i-$begin+1)%3;
            if ($strand eq '-' && $codonPos != 2){
                $codonPos = ($codonPos==0)?1:0;
            }
            $codonPos = ($codonPos == 0)?3:$codonPos;
            $nuc{$c}{$i} = $g."\t".$temp[$i-$begin]."\t".$codonPos; #value: gene name + the reference nucleotide (A, T, G, or C) + codon position (1,2, or 3)
        }
    }
    return %nuc;
}

sub decodeCigar {
    # rewrite the cigar code. for example, if cigar = 10M1I2M, then returns MMMMMMMMMMIMM
    my $cig = shift;
    my @list = split (/\d+/, $cig);
    shift(@list);
    my @list2 = split (/[MIDNSHP]/, $cig);
    my $newString = "";
    for (my $i = 0; $i<scalar(@list); $i++ ){
        $newString .= $list[$i]x$list2[$i];
    }
    #print STDERR $newString." new cigar code\n";
    return $newString; 
}

sub convertQual {
    my @list = split //, shift;
    
    my @scores = ();
    foreach my $b (@list){
        my $sanger_score = ord($b)-33;
        push @scores, $sanger_score;
        #print STDERR $b."\t".$sanger_score."\n";
    }
    return @scores;
}

my ($sampleId, $bamF, $geneF, $minBQ, $minMQ, $minD) = @ARGV;
my %genes = parseGene($geneF);
my %bases = parseBases(\%genes);
my %contigs = ();

for my $c (keys %bases) {
    print $c."===\n";
    $contigs{$c} = 1;
}


my %fTable = (); # allele frequency table. keys: contig_id, pos, the 4 nucleotides. values: counts

print STDERR "parse $bamF\n";
my $counter = 0;
open BAM, "samtools view $bamF | " or die "Can't execute samtools view:$!\n";
while (<BAM>) {
    my ($qname, $flag, $rname, $begin, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $info) = split /\t/;

    next if $rname eq "*";
    next if $mapq < $minMQ;
    next unless $contigs{$rname} == 1;

    my $end = $begin+length($seq)-1;
    # convert char scores to numerical scores
    my @qualScore = convertQual($qual); # base quality scores
    
    
    if ($bases{$rname}{$begin} || $bases{$rname}{$end}){ # if the query read overlaps with an amphora gene
        
        #print "match ".$_;
        $counter++;
        #last if $counter == 100;
        # align the cigar code properly with the base positions of the read
        my $s = decodeCigar($cigar); # this returns an extended cigar alignment (eg 5M2I3M = MMMMMIIMMM)
        
        # bug checking
        #if ($cigar =~ /([NSHP])/){
        #    print STDERR "Warning $cigar: $1 is not properly handled\n";
        #}
        #print STDERR $cigar,"\t", $s."\n";
        
        # this block should be written as a function
        my @b = split //, $seq; # the list of the bases on the read
        my @ci = split //, $s; # the cigar character string
        my @new = (); #converted to the actual alignment
            
	my $readI = 0;
        for (my $cigarI = 0; $cigarI <= scalar(@ci); $cigarI++){ # go through the cigar string
            my $base = "-"; # first initiate each base to N
            if ($ci[$cigarI] eq "D"){ # if the position is a deletion on the read
		$base = "-"; # then use "-" for the base. do not increment $readI
	    }
	    elsif ($ci[$cigarI] eq "H"){ # if the position is hard-clipped
		next; # ignore this position on the query read. do not increment $readI
	    }
	    elsif ($ci[$cigarI] eq "I" || $ci[$cigarI] eq "S"){ # if the position is an insertion on the read or soft-clipped, ignore the base. increment $readI.
		$readI += 1;
		next;
	    }
	    elsif ($qualScore[$readI] < $minBQ) { # if the cigar character is not D, H, S, or I, then inspect base quality
		# a low base quality position
		$base = "-"; # flag this position with an "-". increment $readI
		$readI += 1;
	    }
	    else { # everything looks ok, so we use the read base. increment readI.
		$base = $b[$readI];
		$readI += 1;
            }
            push @new, $base;
        }
                
        # map onto the reference genome
        @b = @new; #reconstructed read sequence including deletions
        for (my $i = $begin; $i<$begin+scalar(@b); $i++){ # this refers to the reference sequence genome
                    
            next unless defined($b[$i-$begin]); # get the nucleotide from the read
            my $nuc = $b[$i-$begin]; # read position
            if ($bases{$rname}{$i}){ # if this nucleotide position is in an amphora gene
                        
                if ($i==0){ #reset $i to some position for debuggin
                #if ($nuc eq "None"){
                    print STDERR "$begin:$i =================================\n";
                    print STDERR join "\t", $nuc,$i,"\n";
                    print STDERR $seq," ",length($seq),"\n"; # debugging purposes
                    print STDERR $s," ",length($s),"\n";
                    print STDERR join "", @b,"\n";
                    print STDERR join "\t",$cigar,$flag,$mapq."\n"; 
                }
                        
                $fTable{$rname}{$i}{$nuc}++ unless $nuc eq "-"; # F denotes a low base quality, so is excluded
                        
            }
        }
    }
}
close BAM;

# final tab-delimited report

print join "\t", "Sample", "Contig", "Position", "Gene", "Ref", "Codon", "Consensus","Allele","Counts", "Frequency"."\n"; 
foreach my $c (keys %fTable){
    foreach my $pos (sort {$a<=>$b} keys %{$fTable{$c}}) {
        my $total = 0;
        my $major = "";

        foreach my $nuc (sort keys %{$fTable{$c}{$pos}}){
            my $counts = $fTable{$c}{$pos}{$nuc};
            next if $counts<= $minD; # this nucleotide does not pass the minimum coverage depth threshold
            $total += $fTable{$c}{$pos}{$nuc};
            if ($major eq ""){
                $major = $nuc;
            }
            elsif ($fTable{$c}{$pos}{$major}<$fTable{$c}{$pos}{$nuc}){
                $major = $nuc;
            }
        }
        foreach my $nuc (sort keys %{$fTable{$c}{$pos}}){
            my $counts = $fTable{$c}{$pos}{$nuc};
            next if $counts <= $minD;
            my $percent = 100*$fTable{$c}{$pos}{$nuc}/$total;
            print join "\t", $sampleId, $c, $pos, $bases{$c}{$pos}, $major, $nuc, $fTable{$c}{$pos}{$nuc},$percent."\n";
        }
    }
}
