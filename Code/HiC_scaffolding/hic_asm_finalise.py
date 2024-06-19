#!/usr/bin/env python3
import sys
import subprocess
import os.path
import argparse
import re
import statistics
from Bio import SeqIO


# Options
cores = 10  # Cores to be used

# Dependencies for this program
dependencies = ['juicer', 'blastn']
repeatmasker = '/usr/local/RepeatMasker-4.1.1-new/RepeatMasker'

# Global variables
# Dmel chromosome names to muller element names dict and arm definition
muller_dmel = {'chrX': 'A', 'chr2L': 'B',
               'chr2R': 'C', 'chr3L': 'D', 'chr3R': 'E', 'chr4': 'F'}
muldir_dmel = {'A': 'L', 'B': 'L', 'C': 'R', 'D': 'L', 'E': 'R', 'F': 'L'}


def main():
    # # # Dependencies and arguments checks # # #
    # Check dependencies
    test_dependencies(dependencies)
    # Handle command line arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-r", "--review_assembly", type=str,
                          help="review.assembly file", required=True, metavar='REV_ASM')
    required.add_argument("-l", "--liftover_agp", type=str,
                          help="liftover.agp file", required=True, metavar='LIFT_AGP')
    required.add_argument("-g", "--genome", type=str,
                          help="genome fasta file", required=True, metavar='GENOME')
    required.add_argument("-m", "--reference_genome", type=str,
                          help="reference fasta file (dm6.fa)", required=True, metavar='REF_GNM')
    args = parser.parse_args()
    # Get input files from command line arguments
    rev_asm_file = args.review_assembly
    lft_agp_file = args.liftover_agp
    genome_file = args.genome
    ref_gnm_file = args.reference_genome
    # Test if input files exist
    test_path_existence(rev_asm_file)
    test_path_existence(lft_agp_file)
    test_path_existence(genome_file)
    test_path_existence(ref_gnm_file)
    # All input tests passed

    # # # Create fasta file from review assembly # # #
    scaff_fas_file = juicer_post(rev_asm_file, lft_agp_file, genome_file)
    # Create genome assembly fasta index
    asm_fai_file = scaff_fas_file + '.fai'
    if not os.path.exists(asm_fai_file):
        subprocess.run(['samtools', 'faidx', scaff_fas_file], capture_output=True)
    # Get scaffold lengths
    scaff_lens = fai_seq_lens(asm_fai_file)
    # Get scaffold sequences
    scaff_seqs = fasta_sequences(scaff_fas_file)
    
    # # # Identify muller elements in scaffolds # # #
    muller_scaffs, muller_majstr = muller_scaffolds(scaff_fas_file, ref_gnm_file, muller_dmel)
    
    # # # Determine chromosome orientation by pericentromeric repeats # # #
    muller_oris = muller_orientations(scaff_fas_file, scaff_lens, muller_scaffs, muldir_dmel)
    
    # Name output prefix
    out_fas_file = re.sub(r'\.fa.*', '.hic_scaffolds_final.fa', genome_file)
    # Open output fasta file
    out_fas = open(out_fas_file, "w")
    # Go through each muller-allocated scaffold
    for scaff in muller_scaffs:
        # Special adjustment for Muller F (by alignment orientation)
        if muller_scaffs[scaff] == 'F':
            muller_oris[scaff] = muller_majstr[scaff]
        # Muller elements stdout
        print(scaff, '->', muller_scaffs[scaff], '<-', muller_oris[scaff])
        # Get the scaffold sequence
        scaff_seq = scaff_seqs[scaff]
        # Get reverse complement if scaffold reversed
        if muller_oris[scaff] == 'rev':
            scaff_seq = reverse_complement(scaff_seq)
        # Print fasta header
        print('>' + muller_scaffs[scaff], file=out_fas)
        # Print fasta sequence
        fas_lines = chop_dna(scaff_seq, 80)
        for fas_line in fas_lines:
            print(fas_line, file=out_fas)
    # Go through each non-muller-allocated scaffold
    for scaff in scaff_seqs:
        if scaff not in muller_scaffs:
            # Print fasta header
            print('>' + scaff, file=out_fas)
            # Print fasta sequence
            fas_lines = chop_dna(scaff_seqs[scaff], 80)
            for fas_line in fas_lines:
                print(fas_line, file=out_fas)


def test_path_existence(path: str):
    """Test whether given path exists, if not: exit"""
    if not os.path.exists(path):
        print("Error: File or directory '" + path +
              "' cannot be found! Exit.\n", file=sys.stderr)
        sys.exit()


def test_dependencies(depends: list):
    """Test whether dependencies are installed, if not: exit"""
    for depend in depends:
        result = subprocess.run(['which', depend], capture_output=True)
        if not result.returncode == 0:
            print("Error: Dependency '" + depend +
                  "' cannot be found! Exit.\n", file=sys.stderr)
            sys.exit()


def fai_seq_lens(file: str) -> dict[str, int]:
    """Take a fasta index file and return a sequence lengths dictionary"""
    fasta_lens: dict[str, int] = {}
    # Open fai file to read
    with open(file, "r") as fa:
        # Read the file line by line
        for line in fa:
            # Strip whitespace from line end
            line = line.rstrip()
            # Get tab field content
            tab_data = line.split("\t")
            # Get seq name and length
            ctg = tab_data[0]
            len = tab_data[1]
            # Save length for seq name
            if ctg not in fasta_lens:
                fasta_lens[ctg] = int(len)
    return fasta_lens


def fasta_sequences(file: str) -> dict[str, str]:
    """Take a fasta file and return a sequence lengths dictionary"""
    fasta_seqs = SeqIO.index(file, "fasta")
    fasta_sqcs: dict[str, str] = {}
    for ctg in fasta_seqs:
        fasta_sqcs[ctg] = str(fasta_seqs[ctg].seq)
    return fasta_sqcs
            

def juicer_post(asm_file: str, agp_file: str, fas_file: str) -> str:
    # Name output prefix
    out_pref = re.sub(r'\.fa.*', '.hic_scaffolds_rev', fas_file)
    # Name of post review files
    post_agp_file = out_pref + '.FINAL.agp'
    post_fas_file = out_pref + '.FINAL.fa'
    post_agp_rn_file = out_pref + '.agp'
    post_fas_rn_file = out_pref + '.fa'
    # Run juicer post if post review fasta file does not yet exist
    if not os.path.exists(post_fas_rn_file):
        # Run juicer post to produce post review fasta file
        os.system(
            f"juicer post -o {out_pref} {asm_file} {agp_file} {fas_file}")
        # Rename post output files
        os.system(f"mv {post_fas_file} {post_fas_rn_file}")
        os.system(f"mv {post_agp_file} {post_agp_rn_file}")
    # Return bed file
    return post_fas_rn_file


def muller_scaffolds(scaff_fas_file: str, ref_gnm_file: str, muller_dmel: dict) -> dict:
    # Muller scaffold dict
    muller_scaffs: dict = {}
    muller_scaffs_majstr: dict = {}
    # Name blast output file
    blast_out = scaff_fas_file + '.' + ref_gnm_file + '.bln'
    # Run juicer post if post review fasta file does not yet exist
    if not os.path.exists(blast_out):
        os.system(
            f"blastn -query {scaff_fas_file} -subject {ref_gnm_file} -out {blast_out} -outfmt \'7\'")
    # Get blast matches on scaffolds per chromosome/muller element
    chr_matches: dict = {}
    chr_matches_plus: dict = {}
    # Open blast output file to read
    with open(blast_out, "r") as blast:
        # Read the file line by line
        for line in blast:
            # Ignore comment lines
            if line.startswith('#'):
                continue
            # Strip whitespace from line end
            line = line.rstrip()
            # Get tab field content
            tab_data = line.split("\t")
            # Get scaffold and chromosome names, and alignment length
            scaff = tab_data[0]
            chrom = muller_dmel[tab_data[1]]
            a_len = tab_data[3]
            s_beg = tab_data[8]
            s_end = tab_data[9]
            # Init dict
            if chrom not in chr_matches:
                chr_matches[chrom] = {}
            if scaff not in chr_matches[chrom]:
                chr_matches[chrom][scaff] = 0
            # Add alignment length to match count
            chr_matches[chrom][scaff] += int(a_len)
            # Plus/plus alignment
            if s_beg < s_end:
                # Init dict
                if chrom not in chr_matches_plus:
                    chr_matches_plus[chrom] = {}
                if scaff not in chr_matches_plus[chrom]:
                    chr_matches_plus[chrom][scaff] = 0
                chr_matches_plus[chrom][scaff] += int(a_len)
    # Go through each chromosome/muller element
    for chrom in sorted(chr_matches):
        # Sort scaffolds by sum of matches
        for scaff, count in sorted(chr_matches[chrom].items(), key=lambda x: x[1], reverse=True):
            # Get stranded alignment counts
            c_plus = chr_matches_plus[chrom][scaff]
            c_minus = count-chr_matches_plus[chrom][scaff]
            #print(scaff, chrom, count, c_plus, c_minus)
            # Save muller element major alignment strand for scaffold
            muller_scaffs_majstr[scaff] = 'for' if c_plus >= c_minus else 'rev'
            # Save muller element name for scaffold
            muller_scaffs[scaff] = chrom
            break
        #print()
    # Return muller scaffold dict
    return muller_scaffs, muller_scaffs_majstr


def muller_orientations(scaff_fas_file: str, scaff_lens: dict, muller_scaffs: dict, muldir_dmel: dict):
    # Muller scaffold orientations
    muller_oris: dict = {}
    # Name repeatmasker output .out file
    rm_out_file = scaff_fas_file + '.out'
    # Run repeatmasker if repeat .out file does not yet exist
    if not os.path.exists(rm_out_file):
        # RM species repeat library
        spec = 'drosophila_flies_genus'
        # Run quick repeatmasker
        os.system(
            f"{repeatmasker} -nolow -no_is -norna -qq --species {spec} -pa {cores} {scaff_fas_file}")
    # Repeat position list dict
    rep_pos_list: dict = {}
    # Open repeatmasker output .out file to read
    with open(rm_out_file, "r") as rm_out:
        # Read the file line by line
        for line in rm_out:
            # Ignore title lines
            if not re.match(r'\s*\d', line):
                continue
            # Strip leading and trailing whitespace from line
            line = line.strip()
            # Get tab field content (split by whitespace)
            tab_data = line.split()
            # Get scaffold and repeat position
            scaff = tab_data[4]
            r_beg = int(tab_data[5])
            r_end = int(tab_data[6])
            r_cla = tab_data[10]
            # Skip satellites
            if r_cla == 'Satellite' or r_cla == 'Unknown':
                continue
            # Init repeat position list for scaffold
            if scaff not in rep_pos_list:
                rep_pos_list[scaff] = []
            # Add repeat position to list connected to scaffold
            for r_pos in range(r_beg, r_end+1):
                rep_pos_list[scaff].append(r_pos)
    # Go through each muller-allocated scaffold
    for scaff in muller_scaffs:
        # Get the relative median repeat position (0-1)
        rel_med_pos = int(statistics.median(
            rep_pos_list[scaff])) / scaff_lens[scaff]
        # Determine 'arm side' (L/R) according to relative repeat clustering
        arm_side = 'L' if rel_med_pos > 0.5 else 'R'
        # Determine scaffold orientation relative to dmel muller element
        mul_ori = 'for' if arm_side == muldir_dmel[muller_scaffs[scaff]] else 'rev'
        # Save scaffold orientation
        muller_oris[scaff] = mul_ori
        #print(scaff, muller_scaffs[scaff], arm_side, muldir_dmel[muller_scaffs[scaff]], mul_ori, rel_med_pos)
    # Return muller scaffold orientations dict
    return muller_oris


def chop_dna(dna_string: str, chop_len: int) -> list:
    return [dna_string[i:i+chop_len] for i in range(0, len(dna_string), chop_len)]


def reverse_complement(dna_sequence: str) -> str:
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'
    }
    return ''.join(complement[base] for base in reversed(dna_sequence))


if __name__ == '__main__':
    main()
