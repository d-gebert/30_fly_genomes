#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $bin_len = 100_000;
# Global variables
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <genome.fa.out> <genome.fa.fai>\n";
unless ($ARGV[0]&&$ARGV[1]) {
   die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $rm_out_file = $ARGV[0];
my $fa_fai_file = $ARGV[1];

my $rep_data = get_repeatmask_data($rm_out_file);
my $chr_lens = seq_lens_from_fai($fa_fai_file);

my $out = open_outfile("$rm_out_file.bed");
foreach my $chr (sort {$a cmp $b} keys %{$rep_data}) {
    foreach my $rep (@{$rep_data->{$chr}}) {
        # Output repeats
        foreach my $i (4,5,6,9) {
            print($out "$rep->[$i]\t");
        }
        my $str = $rep->[8] eq '+' ? '+' : '-';
        print($out "1\t$str\n");
    }
}
close($out);

my %reps_per_bin = ();
# Go through each chromosome
foreach my $chr (sort {$a cmp $b} keys %{$rep_data}) {
    # Go through each repeat
    foreach my $rep (@{$rep_data->{$chr}}) {
        #my $bin = int($rep->[5]/100_000)*100_000;
        #my $rep_len =
        foreach my $pos ($rep->[5]..$rep->[6]) {
            my $bin = int($pos/$bin_len)*$bin_len;
            $reps_per_bin{$chr}{$bin} += 1;
        }
        #$reps_per_bin{$chr}{$bin} += 1;
    }
}

my $out2 = open_outfile("$rm_out_file.bins.bed");
foreach my $chr (sort {$a cmp $b} keys %reps_per_bin) {
    #foreach my $bin (sort {$a <=> $b} keys %{$reps_per_bin{$chr}}) {
    for (my $bin = 0; $bin < $chr_lens->{$chr}; $bin += $bin_len) {
        my $end = $bin + $bin_len;
        if ($end > $chr_lens->{$chr}) {
            $end = $chr_lens->{$chr};
        }
        if ($reps_per_bin{$chr}{$bin}) {
            $reps_per_bin{$chr}{$bin} = ($reps_per_bin{$chr}{$bin} / ($end-$bin+1)) * 100;
        }
        else {
            $reps_per_bin{$chr}{$bin} = 0;
        }

        #print($out2 "$chr\t$bin\t$end\t$reps_per_bin{$chr}{$bin}\n");
        printf($out2 "%s\t%d\t%d\t%.2f\n", $chr,$bin,$end,$reps_per_bin{$chr}{$bin});
    }
}
close($out2);

exit;

################################# subroutines #################################

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file) = @_;
	# Storage variable
	my %rep_data = ();
	# Get file data
	my @repeatmask_data = get_file_data_array($repeatmask_file);
	# Parse repeatmasker file
	foreach my $line (@repeatmask_data) {
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			my $line_s = $line;
			$line_s =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line_s);
			my $chr = $d[4];
			my $cla = $d[10];
			if ($cla =~ /Satellite/ || $cla =~ /Simple/ || $cla =~ /Tandem/ || $cla =~ /Low/) { next; }
			push(@{$rep_data{$chr}},[@d,$line]);
		}
	}
	return \%rep_data;
}

sub seq_lens_from_fai {
	# Take repeatmasker file name
	my($fai_file) = @_;
	# Storage variable
	my %seq_lens = ();
	# Get file data
	my @fai_data = get_file_data_array($fai_file);
	# Parse repeatmasker file
	foreach my $line (@fai_data) {
        # Get line data
        my @d = split(/\t/,$line);
        my $seq_id = $d[0];
        my $sq_len = $d[1];
        $seq_lens{$seq_id} = $sq_len;
	}
	return \%seq_lens;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}