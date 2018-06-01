#!/usr/bin/perl -s
use strict;
use warnings;
use List::Util 'sum';

#############################################################################
# 
# This script extracts UMIs containing single base at mutation point
# (all reads with the same UMI must contain the same base),
# split the UMI's by the nucleotide at the mutation point, count them,
# calculate total number of UMI's and VAFs for the mutated base. 
#
# Usage:
# perl /hpc/pmc_vanboxtel/tools/scripts/smMIP_UMI_filtering.pl \
#	[mutation position] [input bam] [sample name] > [output]
#
# [mutation position]:
#    Comma separated with "chrom,position,refBase,mutBase" 
#    This file is expected to be SORTED
#
# [input bam]:
#	SORTED BAM file.
#	Before the alignment, the smMIP fastq has to be clipped so that the 
#   UMI's could be found back.
#   This BAM file also has to be SORTED.
#
#
# Version 2 (2018/05/22 - Rurika): 
#   Lowered the threashold of when the substitution was called as valid
#   from 100% (no mismatch allowed) to 70%
#
#
##############################################################################

 my $mip_pos_file = $ARGV[0];
 my $bam_file = $ARGV[1];
my $sample_name = $ARGV[2];

 my $wk_dir = '/hpc/pmc_vanboxtel/projects/Rurika_smMIP' ;

open( my $fh_pos, '<:encoding(UTF-8)', $mip_pos_file )
  or die "Could not open file '$mip_pos_file' $!";

print "chr\tpos\tref\talt\ttotalUMI\tA\tT\tG\tC\tVAF\n" ;

while ( my $row_pos = <$fh_pos> ) {

	if ( $row_pos =~ /[a-z]+/ ){
		next ;
	} else {
		chomp $row_pos ;
		my @pos = split /,/, $row_pos ;

		my $mip_chr = $pos[0] ;
		my $mip_pos = $pos[1] ;
		my $mip_ref = $pos[2] ;
		my $mip_mut = $pos[3] ;

			# coordinates to extract reads from BAM file
		my $bam_ext_start = $mip_pos - 200 ;
		my $bam_ext_end = $mip_pos + 200 ;

		my %UMI_counts = (
			A => "0",
			T => "0",
			C => "0",
			G => "0",
			);

			# extracting IDs of UMIs with more than or equal to 5 reads 
		`samtools view $bam_file $mip_chr\:$bam_ext_start\-$bam_ext_end | awk -v n=$mip_pos '{ OFS="\t"; gsub(/ +/, ""); print \$NF, substr(\$10, n+1-\$4, 1)}' | sort | uniq -c | awk '{ if ( \$1 > 4 ) print }' > $wk_dir\/tmp/UMI_ATCG_count_$sample_name\.txt` ;
		`awk '{ print \$2 }' $wk_dir\/tmp/UMI_ATCG_count_$sample_name\.txt | uniq > $wk_dir\/tmp/UMI_cov5_$sample_name\.txt` ;


		open( my $fh_UMIID, '<:encoding(UTF-8)', "$wk_dir\/tmp/UMI_cov5_$sample_name\.txt" )
		  or die "Could not open file '$wk_dir\/tmp/UMI_cov5_$sample_name\.txt' $!";

		    # counting number of each nucleotide at the mutation position for each UMI ID with >4xcov
		    # if more than 70% of the reads consists of the major nucleotide, keep the UMI
		while ( my $row_UMIID = <$fh_UMIID> ) {
			chomp $row_UMIID ; 

			my %read_counts = (
				A => "0",
				T => "0",
				C => "0",
				G => "0",
				NA => "0",
				);

			my @read_count_ATCG = `grep $row_UMIID $wk_dir\/tmp/UMI_ATCG_count_$sample_name\.txt` ;

			my $max_base = 0 ;

			for ( my $i = 0 ; $i < @read_count_ATCG ; $i++ ){

				chomp $read_count_ATCG[$i] ;
				my @read_count_ATCG_spl = split /\s+/, $read_count_ATCG[$i] ;
				my $read_count_num = $read_count_ATCG_spl[1] ;
				my $read_count_base = $read_count_ATCG_spl[3] ;

					# count only if a base is identified (sometimes it is empty
					# could be reads not covering the mutation position)
				if ( $read_count_base ){
					$read_counts{$read_count_base} = $read_count_num ;
				}
			}

			my $total_reads = sum values %read_counts;

			$max_base = max(%read_counts) ;
			if ( $total_reads ){
				my $max_base_freq = $read_counts{$max_base} / $total_reads ;
				if ( $max_base_freq >= 0.7 ){
					$UMI_counts{$max_base}++ ;
				}
			}
		}

		my @UMI_counts_ATCG = ( $UMI_counts{'A'}, $UMI_counts{'T'}, $UMI_counts{'G'}, $UMI_counts{'C'});
		my $total_UMIs = sum values %UMI_counts;
		my $VAF = "NA" ;

		if ( $total_UMIs ){
			$VAF = $UMI_counts{$mip_mut} / $total_UMIs ;
		} 		

		print join ( "\t", @pos ), "\t$total_UMIs\t", join ( "\t", @UMI_counts_ATCG ), "\t$VAF\n";
	}
}





sub max {
    my (%data) = @_;
 
    my $max;
    while (my ($key, $value) = each %data) {
        if (not defined $max) {
            $max = $key;
            next;
        }
        if ($data{$max} < $value) {
            $max = $key;
        }
    }
    return $max;
}
