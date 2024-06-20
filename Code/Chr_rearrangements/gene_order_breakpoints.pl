#!/usr/bin/perl
use strict;
use warnings;

# Command line arguments
my ($ref_genome, $bin_size) = @ARGV;

# Check for correct input
die "Usage: $0 <reference genome> <bin size>\n" unless $ref_genome && $bin_size;

# Input and output file names
my $input_filename = 'syntenicBlock_coordinates.csv';
my $output_filename = "gene_order_breakpoints_${ref_genome}.csv";

# Open CSV file for reading
open my $fh, '<', $input_filename or die "Could not open '$input_filename': $!\n";

# Data storage for breakpoints and gene order ranges
my %bins;

# Skip header line
my $header_line = <$fh>;

# Process each line in the input CSV file
while (my $line = <$fh>) {
    chomp $line;
    # Split the line into fields based on comma separator
    my @fields = split /,/, $line;

    # Extract relevant fields
    my ($genome1, $genome2, $chr1, $chr2, $blkID, $startBp1, $endBp1, $startOrd1, $endOrd1) = @fields;

    # Process only the relevant reference genome
    if ($genome1 eq $ref_genome) {
        # Calculate start and end bins
        my $start_bin = int(($startOrd1 - 1) / $bin_size);
        my $end_bin = int(($endOrd1 - 1) / $bin_size);
        
        # Update bin information for the chromosome
        for (my $bin = $start_bin; $bin <= $end_bin; $bin++) {
            $bins{$chr1}{$bin}{count}++;
            # Update minimum and maximum gene order ranges
            $bins{$chr1}{$bin}{min} = $startOrd1 if (!exists $bins{$chr1}{$bin}{min} || $startOrd1 < $bins{$chr1}{$bin}{min});
            $bins{$chr1}{$bin}{max} = $endOrd1 if (!exists $bins{$chr1}{$bin}{max} || $endOrd1 > $bins{$chr1}{$bin}{max});
        }
    }
}

# Close input file handle
close $fh;

# Open output file for writing
open my $out, '>', $output_filename or die "Could not open '$output_filename': $!\n";

# Write header line to the output file
print $out "Chromosome,Bin,Count,GeneOrderRange\n";

# Write bin information to the output file
foreach my $chr (sort keys %bins) {
    foreach my $bin (sort { $a <=> $b } keys %{ $bins{$chr} }) {
        my $range = $bins{$chr}{$bin}{min} . ".." . $bins{$chr}{$bin}{max};
        my $count = $bins{$chr}{$bin}{count};
        print $out "$chr,$bin,$count,$range\n";
    }
}

# Close output file handle
close $out;

# Inform the user that the process is complete
print "Gene order based breakpoint data written to '$output_filename'\n";
