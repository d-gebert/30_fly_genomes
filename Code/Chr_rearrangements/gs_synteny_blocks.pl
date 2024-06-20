#!/usr/bin/perl
use strict;
use warnings;
use Bio::TreeIO;

# Global constants
my $bin_len = 0.05;
# Global variables
# Options
$|=1; #Autoflush

# Muller elements
my @mullers = ('A','B','C','D','E','F');

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <genespace_results_dir> <muller_lengths.txt> <tree.newick>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $gs_dir = $ARGV[0];
my $l_file = $ARGV[1];
my $newick = $ARGV[2];

# Load newick tree
my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $newick);
my $tree = $treeio->next_tree;

# Muller elements lengths
my %mul_lens = ();
# Get file data
my @mul_lens_data = get_file_data_array($l_file);
# Chromosome list
my @chrs = ();
# Parse muller lengths file
for my $line (@mul_lens_data) {
    # Skip title line
    if ($line eq $mul_lens_data[0]) {
        # Chromosome names
        @chrs = split(/\t/,$line);
        shift(@chrs);
        next;
    }
    # Split line
    my @dat = split(/\t/,$line);
    my $spe = shift(@dat);
    # Go through each chromosome index
    for my $idx (0..$#chrs) {
        # Species name and genome file base
        my $chr = $chrs[$idx];
        my $len = $dat[$idx];
        $mul_lens{$spe}{$chr} = $len;
    }
}

# Get path to syntenic block coordinates file
my $syn_block_file = "$gs_dir/results/syntenicBlock_coordinates.csv";
# Get file data
my @syn_block_data = get_file_data_array($syn_block_file);
# Number of blocks per muller element hash
my %n_blocks_per_dist = ();
# Block sizes hash
my %block_lens = ();
my %block_lens_per_dist = ();
# Parse syntenic block coordinates data
for my $line (@syn_block_data) {
    # Skip title line
    next if $line eq $syn_block_data[0];
    # Split line
    my @dat = split(/,/,$line);
    # Block properties
    my $spe_x = $dat[0];
    my $spe_y = $dat[1];
    my $chr_x = $dat[2];
    my $chr_y = $dat[3];
    my $blkid = $dat[4];
    my $beg_x = $dat[5];
    my $end_x = $dat[6];
    my $orien = $dat[19];
    my $beg_y = $dat[20];
    my $end_y = $dat[21];
    my ($pair_id) = ($blkid =~ /(\w+):\s/);
    my ($blk_num) = ($blkid =~ /\w+:\s(\d+)/);
    # Block size
    my $len_x = $end_x - $beg_x + 1;
    # Skip if no interspecies pair
    next if !$pair_id;
    $pair_id = "$spe_x.$spe_y";
    # Get evolutionary distance between species
    my $evo_dist = calculate_combined_branch_length($tree, $spe_x, $spe_y);
    my $dist_bin = int($evo_dist/$bin_len)*$bin_len;
    $n_blocks_per_dist{$chr_x}{$dist_bin}{$pair_id} += 1/$mul_lens{$spe_x}{$chr_x}*1000000 if $mul_lens{$spe_x}{$chr_x};
    $n_blocks_per_dist{$chr_y}{$dist_bin}{$pair_id} += 1/$mul_lens{$spe_y}{$chr_y}*1000000 if $mul_lens{$spe_y}{$chr_y};
    # Save block size
    push(@{$block_lens{$chr_x}},$len_x);
    push(@{$block_lens_per_dist{$chr_x}{$dist_bin}},$len_x);
}

# Open output file
my $outfile1 = 'muller_block_sizes.txt';
my $out1 = open_outfile($outfile1);
# Get max index number n
my $n = 0;
for my $chr (@mullers) {
    $n = scalar(@{$block_lens{$chr}}) > $n ? scalar(@{$block_lens{$chr}}) : $n;
}
# Title line
for my $chr (@mullers) {
    print($out1 "$chr\t");
}
print($out1 "\n");
# Print each index as row
for my $i (0..$n-1) {
    for my $chr (@mullers) {
        $block_lens{$chr}[$i] = 'NA' if !$block_lens{$chr}[$i];
        print($out1 "$block_lens{$chr}[$i]\t");
    }
    print($out1 "\n");
}
close($out1);

# Scale for evolutionary distance
my @dist_range = map { $_ * 0.05 } (0..8);

# Open output file
my $outfile2 = 'muller_breaks_per_dist.txt';
my $out2 = open_outfile($outfile2);
# Title line
print($out2 "mull\tdist\tbrks\n");
for my $chr (@mullers) {
    for my $dist (@dist_range) {
        for my $pair (sort keys %{$n_blocks_per_dist{$chr}{$dist}}) {
            $n_blocks_per_dist{$chr}{$dist}{$pair} = 0 unless $n_blocks_per_dist{$chr}{$dist}{$pair};
            print($out2 "$chr\t$dist\t$n_blocks_per_dist{$chr}{$dist}{$pair}\n");
        }
    }
}
close($out2);

# Open output file
my $outfile3 = 'muller_block_sizes_per_dist.txt';
my $out3 = open_outfile($outfile3);
# Title line
print($out3 "mull\tdist\tbrks\n");
for my $chr (@mullers) {
    for my $dist (@dist_range) {
        for my $len (@{$block_lens_per_dist{$chr}{$dist}}) {
            print($out3 "$chr\t$dist\t$len\n");
        }
    }
}
close($out3);

# Means
for my $dist (@dist_range) {
    print("\t$dist");
}
print("\n");
for my $chr (@mullers) {
    print("$chr");
    for my $dist (@dist_range) {
        my @dist_list = ();
        for my $pair (sort keys %{$n_blocks_per_dist{$chr}{$dist}}) {
            $n_blocks_per_dist{$chr}{$dist}{$pair} = 0 unless $n_blocks_per_dist{$chr}{$dist}{$pair};
            push(@dist_list,$n_blocks_per_dist{$chr}{$dist}{$pair});
        }
        my $dist_mean = get_mean(\@dist_list);
        $dist_mean = round_float($dist_mean, 2);
        print("\t$dist_mean");
    }
    print("\n");
}
print("\n");
# Std divs
for my $dist (@dist_range) {
    print("\t$dist");
}
print("\n");
for my $chr (@mullers) {
    print("$chr");
    for my $dist (@dist_range) {
        my @dist_list = ();
        for my $pair (sort keys %{$n_blocks_per_dist{$chr}{$dist}}) {
            $n_blocks_per_dist{$chr}{$dist}{$pair} = 0 unless $n_blocks_per_dist{$chr}{$dist}{$pair};
            push(@dist_list,$n_blocks_per_dist{$chr}{$dist}{$pair});
        }
        my $dist_sdev = get_sample_standard_deviation(\@dist_list);
        $dist_sdev = round_float($dist_sdev, 2);
        print("\t$dist_sdev");
    }
    print("\n");
}
print("\n");

exit;

################################# subroutines #################################

# Function to calculate distance from a node to an ancestor
sub distance_to_ancestor {
    my ($node, $ancestor) = @_;
    my $distance = 0;

    # Traverse up the tree from node to ancestor
    while (defined $node && (!$ancestor || $node ne $ancestor)) {
        $distance += $node->branch_length if $node->branch_length;
        $node = $node->ancestor;
    }

    return $distance;
}

# Function to calculate combined branch length
sub calculate_combined_branch_length {
    my ($tree, $species1_name, $species2_name) = @_;

    # Find nodes for both species
    my $species1_node = $tree->find_node(-id => $species1_name);
    my $species2_node = $tree->find_node(-id => $species2_name);

    # Handle case where species are not found
    unless ($species1_node && $species2_node) {
        warn "One or both species not found in the tree.\n";
        return;
    }

    # Find least common ancestor (LCA)
    my $lca = $tree->get_lca(-nodes => [$species1_node, $species2_node]);

    # Calculate combined branch length
    my $distance1 = distance_to_ancestor($species1_node, $lca);
    my $distance2 = distance_to_ancestor($species2_node, $lca);

    return $distance1 + $distance2;
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

# Calculate mean of array values
sub get_mean {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	if ($N == 0) { return 0 }
	# Sum values
	my $sum = get_sum($array);
	# Calculate mean value
	my $mean = $sum/$N;
	# Return mean
	return $mean;
}

# Round decimal number to desired decimal places
sub round_float {
	# Take value and number of decimal places
	my($val, $n_dec) = @_;
	# Truncate to desired decimal places
	my $trunc_val = int($val*(10**$n_dec)+0.5)/(10**$n_dec);
	# Return truncated value
	return $trunc_val;
}

# Calculate sample variance of array values
sub get_sample_variance {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	# Get mean value
	my $mean = get_mean($array);
	# Sum squared differences
	my $diff_sum = 0;
	grep { $diff_sum += ($_-$mean)**2 } @{$array};
	# Calculate sample variance
	my $variance = $diff_sum/($N-1);
	# Return sample variance
	return $variance;
}

# Calculate sample standard deviation of array values
sub get_sample_standard_deviation {
	# Take array list
	my($array) = @_;
	# Get sample variation
	my $variance = get_sample_variance($array);
	# Calculate standard deviation from variance
	my $std_dev = sqrt($variance);
	# Return standard deviation
	return $std_dev;
}