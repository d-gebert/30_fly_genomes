use strict;
use warnings;

# Species (groups)
my @subg1 = ('Dmel','Dmau','Dsim','Dsec','Dtei','Dyak','Dere','Dtak');
my @subg2 = ('Dboc','Djam','Dkik','Dtri','Druf');
my @subg3 = ('Dmal_pal','Dmal_mal','Dbip','Dpar','Dpse_pse','Dana');
my @subg4 = ('Dper','Dpse','Dsub');
my @subg5 = ('Dequ','Dpau','Dwil','Dtro','Dins');
my @subg6 = ('Dvir','Dame','Dlit');
my @specs = (@subg1,@subg2,@subg3,@subg4,@subg5,@subg6);

# Set default threshold
my $max_breaks = 1;

# Input files
my $file_break_points = 'gene_breakpoints_Dmel.csv';
my $file_pangenes = 'Dmel_pangenes.txt';

# Get break point count valleys
my @break_valleys = get_break_valleys($file_break_points, $max_breaks);

# Get pangenes
my %pangenes = get_pangenes($file_pangenes);

# Valley index
my $i = 0;
# Go through each valley
for my $valley (@break_valleys) {
    # Get valley data
    my ($chr,$beg_bin,$end_bin,$size,$beg_og,$end_og) = @{$valley};
    # Get number of ref spec genes in valley
    my $count = 0;
    my $n_gen = $end_og-$beg_og+1;
    for my $go ($beg_og..$end_og) {
        $count++ if $pangenes{$go}{Dmel};
    }
    # Skip region with >1/3 NAs
    next if $count/$n_gen < 0.67;
    # Increment valley index
    $i++;
    # Region line
    print("REGION $i - chr:$chr bins:$beg_bin-$end_bin n_bins:$size gene_nums:$beg_og-$end_og\n");
    # Go through each species
    for my $sp (@specs) {
        print("Dmel#\t$sp#\tgene\tchr\tbeg\tend\n");
        for my $go ($beg_og..$end_og) {
            if ($pangenes{$go}{$sp}) {
                print("$go");
                for my $d (@{$pangenes{$go}{$sp}}) {
                    print("\t$d");
                }
                print("\n");
            }
            else {
                print("$go\tNA\tNA\tNA\tNA\tNA\n");
            }
        }
        print("\n");
    }
    print("\n");
}

exit;

################################# subroutines #################################

sub get_break_valleys {
    # Take in file name and threshold
    my ($file, $threshold) = @_;
    # Breakopoint valleys storage
    my @break_valleys = ();
    # Loop variables
    my $start_bin;
    my $end_bin;
    my $start_gene;
    my $end_gene;
    my $current_chr;
    # Open the CSV file
    open(my $fh, '<', $file) or die "Could not open file '$file' $!";
    # Read the header
    my $header = <$fh>;
    # Parse file
    while (my $line = <$fh>) {
        # Remove newline
        chomp $line;
        # Get line data
        my ($chr, $bin, $count, $range) = split(/,/, $line);
        my ($gene_x, $gene_y) = split(/\.\./, $range);
        # If the count is above threshold
        if ($count <= $threshold) {
            # If start variable is not defined
            # We have just entered a subthreshold count valley
            if (!defined $start_bin) {
                # Define start bin as current bin number
                $start_bin = $bin;
                $start_gene = $gene_x;
                # Define chromosome as current chromosome
                $current_chr = $chr;
            }
            # This will be updated as long as we remain in the valley
            # Define end bin as current bin number
            $end_bin = $bin;
            $end_gene = $gene_y;
        }
        # If the count is not above threshold
        else {
            # If start variable is defined
            # We have just passed a subthreshold count valley
            if (defined $start_bin) {
                # Get length of bin stretch, i.e, subthreshold count valley
                my $len = $end_bin-$start_bin+1;
                # Save break point valley that was just exited
                push(@break_valleys,[$current_chr,$start_bin,$end_bin,$len,$start_gene,$end_gene]);
                # Delete previous start and end bins
                undef $start_bin;
                undef $end_bin;
            }
            # If no start defined this line will be ignored
        }
    }
    # Print the last range if it exists
    if (defined $start_bin) {
        # Get length of bin stretch, i.e, subthreshold count valley
        my $len = $end_bin-$start_bin+1;
        # Save break point valley that was just exited
        push(@break_valleys,[$current_chr,$start_bin,$end_bin,$len,$start_gene,$end_gene]);
    }
    # Close file handle
    close($fh);
    # Return break valleys
    return @break_valleys;
}

sub get_pangenes {
    # Take in file name and threshold
    my ($file) = @_;
    # Breakopoint valleys storage
    my %pangenes = ();
    # Open the CSV file
    open(my $fh, '<', $file) or die "Could not open file '$file' $!";
    # Read the header
    my $header = <$fh>;
    # Parse file
    while (my $line = <$fh>) {
        # Remove newline
        chomp $line;
        # Get line data
        my @d = split(/\t/, $line);
        # Save pan gene info for each species
        my ($go,$sp,$gene,$chr,$beg,$end,$go_sp) = ($d[3],$d[5],$d[8],$d[9],$d[10],$d[11],$d[12]);
        $pangenes{$go}{$sp} = [$go_sp,$gene,$chr,$beg,$end];
    }
    # Return pangenes
    return %pangenes;
}