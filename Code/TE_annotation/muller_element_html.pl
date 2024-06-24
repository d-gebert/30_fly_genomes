#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my @rgb1 = (32,56,100);
my @rgb2 = (68,114,196);
my @rgb3 = (180,199,231);
my $chr_len_min_frct = 0.01;
my $v_space = 50;
my $x_space = 10;
my @file_bases = (
    'D.sechellia',
	'D.simulans',
	'D.mauritiana',
	'D.melanogaster',
	'D.teissieri.2733',
	'D.yakuba',
	'D.erecta',
	'D.takahashii',
	'D.jambulina',
	'D.bocqueti',
	'D.kikkawai',
	'D.triauraria',
	'D.rufa',
	'D.bipectinata',
	'D.parabipectinata',
	'D.m.malerkotliana',
	'D.m.pallens',
	'D.p.pseudoananassae',
	'D.ananassae',
	'D.persimilis',
	'D.pseudoobscura',
	'D.subobscura',
	'D.equinoxialis',
	'D.paulistorum.L12',
	'D.willistoni.00',
	'D.tropicalis',
	'D.insularis',
	'D.virilis',
	'D.americana',
	'D.littoralis'
);
my @specs = (
    'Dsec',
	'Dsim',
	'Dmau',
	'Dmel',
	'Dtei',
	'Dyak',
	'Dere',
	'Dtak',
	'Djam',
	'Dboc',
	'Dkik',
	'Dtri',
	'Druf',
	'Dbip',
	'Dpar',
	'Dmal_mal',
	'Dmal_pal',
	'Dpse_pse',
	'Dana',
	'Dper',
	'Dpse',
	'Dsub',
	'Dequ',
	'Dpau',
	'Dwil',
	'Dtro',
	'Dins',
	'Dvir',
	'Dame',
	'Dlit'
);
# Global variables
my $total_w = 2000;
my $total_h = 2000;
# Options
$|=1; #Autoflush

# Muller elements
my @mullers = ('A','B','C','D','E','F');
my %mullers = map { $_ => 1 } @mullers;

# Program name
print("\n--- $0 ---\n");

my $max_chr_len = 40_000_000;

# Open output file
my $outfile = "muller_rep_bins.html";
my $out = open_outfile($outfile);

# Print top part of svg html file
print($out "
<!DOCTYPE html>
<html>
<body>

<svg width=\"$total_w\" height=\"$total_h\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $total_w $total_h\" id=\"tooltip-svg\">
	<style>
		#tooltip {
			dominant-baseline: hanging;
		}
        .small { font: 12px sans-serif; }
        .medium { font: 14px sans-serif; }
        .large { font: 16px sans-serif; }
        .title { font: 18px sans-serif; }
	</style>

");

# Draw chromosomes and pirna cluster loci
my $chr_y_top = 0;
# Chromosome width factor
my $chr_w_f = 140;
# Chromosome x space
my $chr_x_space = 220;
# Go through each chromosome
my $j = 0;
foreach my $chr (@mullers) {
    my $txt_x = $x_space+($j*$chr_x_space);
    print($out "\t\t<text x=\"$txt_x\" y=\"15\" fill=\"black\" class=\"title\">$chr</text>\n");
    $j++;
}
for my $idx (0..$#specs) {
    $chr_y_top += 25;
    # Species name and genome file base
    my $spec = $specs[$idx];
    my $base = $file_bases[$idx];
    # Collect command line arguments
    my $genome_file = "genome_fas/$base.hic_scaffolds_final.fa";
    my $rm_bed_file = "$spec.HiTE.out.bins.bed";
    # Get chromosome lengths
    my $chr_lens = get_fasta_seq_lens($genome_file,1);
    # Get regional start and stop positions
    my $reg_beg = 1;
    my $reg_end = $max_chr_len;
    # Get genome rm bins
    my $rm_bin_data = get_tab_fields_id0($rm_bed_file);
    # Print species name
    my $txt_y = $chr_y_top+33;
    print($out "\t\t<text x=\"12600\" y=\"$txt_y\" fill=\"black\" class=\"large\">$spec</text>\n");
    # Initialize chromosome i
    my $i = 0;
    # Go through each chromosome
    foreach my $chr (@mullers) {
        # Get chromosome length
        my $chr_len = $chr_lens->{$chr};
        # Get chromosome coordinate system properties
        my $chr_y = $chr_y_top;
        my $chr_x = $x_space+($i*$chr_x_space);
        my $chr_w = ($chr_len-$reg_beg+1)/(($reg_end-$reg_beg+1)/$chr_w_f);
        # Repeat content bins
        my $chr_b_y = $chr_y + ($v_space/2.5);
        # Print chromosome rectangle
        print($out "\t\t<rect x=\"$chr_x\" y=\"$chr_b_y\" width=\"$chr_w\" height=\"15\" style=\"fill:rgb(250,250,250);stroke-width:1;stroke:rgb(250,250,250)\" />\n");
        # Print each bin
        foreach my $bin (@{$rm_bin_data->{$chr}}) {
            # Bin coordinates
            my $bin_beg = $bin->[1];
            my $bin_end = $bin->[2];
            my $bin_len = $bin_end-$bin_beg+1;
            # Get bin coordinate system properties
            my $bin_x = $chr_x+($bin_beg-$reg_beg+1)/(($reg_end-$reg_beg+1)/$chr_w_f);
            my $bin_w = ($bin_len-$reg_beg+1)/(($reg_end-$reg_beg+1)/$chr_w_f);
            # Repeat share
            my $rep_prc = $bin->[3];
            my $factor = (100-$rep_prc)/100;
            my @rgb = (235*$factor,235*$factor,235*$factor);
            # Print repeat bin rectangle
            print($out "\t\t<rect x=\"$bin_x\" y=\"$chr_b_y\" width=\"$bin_w\" height=\"15\" style=\"fill:rgb($rgb[0],$rgb[1],$rgb[2]);stroke-width:1;stroke:rgb($rgb[0],$rgb[1],$rgb[2])\"/>\n");
        }
        # Increment chromosome i
        $i++;
    }
}

# Print bottom part of svg html file
print($out "

    <g id=\"tooltip\" visibility=\"hidden\" >
		<rect x=\"2\" y=\"2\" width=\"100\" height=\"34\" fill=\"black\" opacity=\"0.4\" rx=\"2\" ry=\"2\"/>
		<rect width=\"100\" height=\"34\" fill=\"white\" rx=\"2\" ry=\"2\"/>
		<text x=\"4\" y=\"17\" class=\"medium\">Tooltip</text>
	</g>

	<script type=\"text/ecmascript\"><![CDATA[
		(function() {
			var svg = document.getElementById('tooltip-svg');
			var tooltip = svg.getElementById('tooltip');
			var tooltipText = tooltip.getElementsByTagName('text')[0];
			var tooltipRects = tooltip.getElementsByTagName('rect');
			var triggers = svg.getElementsByClassName('tooltip-trigger');
			for (var i = 0; i < triggers.length; i++) {
				triggers[i].addEventListener('mousemove', showTooltip);
				triggers[i].addEventListener('mouseout', hideTooltip);
			}
			function showTooltip(evt) {
				var CTM = svg.getScreenCTM();
				var x = (evt.clientX - CTM.e + 6) / CTM.a;
				var y = (evt.clientY - CTM.f + 20) / CTM.d;
				tooltip.setAttributeNS(null, \"transform\", \"translate(\" + x + \" \" + y + \")\");
				tooltip.setAttributeNS(null, \"visibility\", \"visible\");
				tooltipText.firstChild.data = evt.target.getAttributeNS(null, \"data-tooltip-text\");
				var length = tooltipText.getComputedTextLength();
				for (var i = 0; i < tooltipRects.length; i++) {
					tooltipRects[i].setAttributeNS(null, \"width\", length + 8);
				}
			}
			function hideTooltip(evt) {
				tooltip.setAttributeNS(null, \"visibility\", \"hidden\");
			}
		})()
    ]]></script>

</svg>

</body>
</html>
");

exit;

################################# subroutines #################################

sub get_tab_fields_id0 {
	# Take name of bin tree file
	my($infile, $id_i) = @_;
	# Set 0 as default index
	unless ($id_i) { $id_i = 0 }
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
        # Skip title line
        #if ($line eq $in_data[0]) { next; }
        # Get line data
        my @d = split(/\t/,$line);
        # Get id
        my $id = $d[$id_i];
        # Save data fields
        push(@{$data_fields{$id}},\@d);
    }
	# Return data fields
	return \%data_fields;
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

sub get_fasta_seq_lens {
	# Take fasta file name
	my($file,$short) = @_;
	# Variables
	my $name = '';
	my %seq_lens = ();
	# Open file
	my $in = open_infile($file);
	# Extract sequence and save with name
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		if ($line =~ /^>/) {
			($name) = ($line =~ />(.*)$/);
			if ($short) { $name =~ s/\s.*// }
		} else {
			$seq_lens{$name} += length($line);
		}
	}
	# Return sequence lengths
	return \%seq_lens;
}