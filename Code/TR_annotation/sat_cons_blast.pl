#!/usr/bin/perl
use strict;
use warnings;

# Global constants
# Global variables
# Options
$|=1; #Autoflush

# Muller elements
my @mullers = ('A','B','C','D','E','F');
my %mullers = map { $_ => 1 } @mullers;
# Species
my @species = qw(
    D.m.pallens
    D.m.malerkotliana
    D.bipectinata
    D.parabipectinata
    D.p.pseudoananassae
    D.ananassae
);

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <CEN.cons.fa> <CEN.locs.txt>\n";
unless ($ARGV[0]&&$ARGV[1]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $cons_file = $ARGV[0];
my $locs_file = $ARGV[1];

# Satellite consensus families
my @fams = ('FAM1', 'FAM2a', 'FAM2b', 'FAM3');

# Get locs file data
my @locs_data = 
get_file_data_array($locs_file);
# Centromere locs hash
my %cen_locs = ();
# Current species
my $spec = '';
# Parse locs file data
for my $line (@locs_data) {
    # Species line
    if ($line =~ /^#/) {
        ($spec) = ($line =~ /#(.*)\|/);
    }
    elsif ($line !~ /^$/) {
        # Coordinates
        my($chr) = ($line =~ /^(\w+)\:/);
        my($beg) = ($line =~ /^\w+\:(\d+)-/);
        my($end) = ($line =~ /^\w+\:\d+-(\d+)/);
        # Save coordinates
        if ($beg && $end) {
            $cen_locs{$spec}{$chr} = [$beg,$end];
        }
    }
}

# Storage variable for satellite hits
my %sat_hits = ();
my %cen_hits = ();
# Go through each species
for my $spec (@species) {
    # Storage variable for satellite hit positions
    my %sat_hit_poss = ();
    my %cen_hit_poss = ();
    # Genome fasta file
    my $gnm_file = "$spec.hic_scaffolds_final.fa";
    # Blast output file
    my $bln_file = "$cons_file.$spec.bln";
    # Run blastn
    if (!-e $bln_file) {
        system("blastn -query $cons_file -subject $gnm_file -out $bln_file -dust no");
    }
    # Blast output file (outfmt 7)
    my $tbl_file = "$cons_file.$spec.bln.tbl";
    # Run blastn (outfmt 7)
    if (!-e $tbl_file) {
        system("blastn -query $cons_file -subject $gnm_file -out $tbl_file -outfmt 7 -dust no");
    }
    # Blast output file (outfmt 10)
    my $csv_file = "$cons_file.$spec.bln.csv";
    # Run blastn (outfmt 10)
    if (!-e $csv_file) {
        system("blastn -query $cons_file -subject $gnm_file -out $csv_file -outfmt 10 -dust no");
    }
    # Error: no blast output/run
    if (!-e $tbl_file) { die("\nCould not find blast output or run blast\n"); }
    # Extract blast output data
    my @bln_data = 
    get_file_data_array($tbl_file);
    # Parse blast output file
    for my $line (@bln_data) {
        # Skip info lines
        next if $line =~ /^#/;
        # Get hit data
        my @dat = split(/\t/,$line);
        # Get sat family and hit coordinates
        my($fam) = ($dat[0] =~ /(\w+)-/);
        my $chr = $dat[1];
        my $beg = $dat[8] < $dat[9] ? $dat[8] : $dat[9];
        my $end = $dat[8] > $dat[9] ? $dat[8] : $dat[9];
        # Save hit positions
        for my $pos ($beg..$end) {
            $sat_hit_poss{$spec}{$fam}{$chr}{$pos} = 1;
        }
        # Centromere coordinates
        next if !$cen_locs{$spec};
        next if !$cen_locs{$spec}{$chr};
        my $cen_beg = $cen_locs{$spec}{$chr}->[0];
        my $cen_end = $cen_locs{$spec}{$chr}->[1];
        # In cen bool
        my $in_cen = 0;
        if ($beg >= $cen_beg && $end <= $cen_end) {
            $in_cen = 1;
        }
        if ($in_cen) {
            # Save hit positions
            for my $pos ($beg..$end) {
                $cen_hit_poss{$spec}{$fam}{$chr}{$pos} = 1;
            }
        }
    }
    # Count put hits per chromosome
    for my $fam (@fams) {
        for my $chr (sort keys %{$sat_hit_poss{$spec}{$fam}}) {
            # Hit base pairs
            my $hit_bps = keys %{$sat_hit_poss{$spec}{$fam}{$chr}};
            # Unplaced scaffold
            $chr = 'U' if !$mullers{$chr};
            # Save number of hit base pairs for this chromosome
            $sat_hits{$spec}{$fam}{$chr} += $hit_bps;
        }
        for my $chr (sort keys %{$cen_hit_poss{$spec}{$fam}}) {
            # Hit base pairs
            my $hit_bps = keys %{$cen_hit_poss{$spec}{$fam}{$chr}};
            # Unplaced scaffold
            $chr = 'U' if !$mullers{$chr};
            # Save number of hit base pairs for this chromosome
            $cen_hits{$spec}{$fam}{$chr} += $hit_bps;
        }
    }
}

# Title line
for my $fam (@fams) {
    print("\t$fam");
}
print("\n");
# Go through each species
for my $spec (@species) {
    print("$spec");
    for my $fam (@fams) {
        my $fam_hits = 0;
        for my $chr (sort keys %{$sat_hits{$spec}{$fam}}) {
            my $hit_bps = $sat_hits{$spec}{$fam}{$chr};
            $fam_hits += $hit_bps;
        }
        print("\t$fam_hits");
    }
    print("\n");
}
print("\n");

# Title line
for my $fam (@fams) {
    print("\t$fam");
}
print("\n");
# Go through each species
for my $spec (@species) {
    print("$spec");
    for my $fam (@fams) {
        my $fam_hits = 0;
        for my $chr (sort keys %{$sat_hits{$spec}{$fam}}) {
            next if $chr ne 'U';
            my $hit_bps = $sat_hits{$spec}{$fam}{$chr};
            $fam_hits += $hit_bps;
        }
        print("\t$fam_hits");
    }
    print("\n");
}
print("\n");

# Title line
for my $fam (@fams) {
    print("\t$fam");
}
print("\n");
# Go through each species
for my $spec (@species) {
    print("$spec");
    for my $fam (@fams) {
        my $fam_hits = 0;
        for my $chr (sort keys %{$cen_hits{$spec}{$fam}}) {
            my $hit_bps = $cen_hits{$spec}{$fam}{$chr};
            $fam_hits += $hit_bps;
        }
        print("\t$fam_hits");
    }
    print("\n");
}
print("\n");

# Go through each sat fam
for my $fam (@fams) {
    print("$fam\tA\tB/C\tD/E\tF\tUn\n");
    # Go through each species
    for my $spec (@species) {
        for my $chr ('A','B','C','D','E','F', 'U') {
            $sat_hits{$spec}{$fam}{$chr} = 0 if !$sat_hits{$spec}{$fam}{$chr};
        }
        my $mulA  = $sat_hits{$spec}{$fam}{'A'};
        my $mulBC = $sat_hits{$spec}{$fam}{'B'}+$sat_hits{$spec}{$fam}{'C'};
        my $mulDE = $sat_hits{$spec}{$fam}{'D'}+$sat_hits{$spec}{$fam}{'E'};
        my $mulF  = $sat_hits{$spec}{$fam}{'F'};
        my $Unpld = $sat_hits{$spec}{$fam}{'U'};
        print("$spec\t$mulA\t$mulBC\t$mulDE\t$mulF\t$Unpld\n");
    }
}
print("\n");

exit;

################################# subroutines #################################

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