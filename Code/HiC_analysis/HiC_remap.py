#!/usr/bin/env python3
import sys
import subprocess
import os.path
import re
import argparse


# Options
cores = 10  # Cores to be used
map_q = 1  # Initial minimal map quality
mapper = 'bwa'

# Dependencies for this program
dependencies = ['bwa', 'samtools', 'bedtools', 'trim_galore', 'sambamba',
                'filter_five_end.pl', 'two_read_bam_combiner.pl']


def main():
    # # # Dependencies and arguments checks # # #
    # Check dependencies
    test_dependencies(dependencies)
    # Handle command line arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-1", "--read1", type=str, help="read1 fastq file", required=True, metavar='READ1')
    required.add_argument("-2", "--read2", type=str, help="read2 fastq file", required=True, metavar='READ2')
    required.add_argument("-g", "--genome", type=str, help="genome fasta file", required=True, metavar='GENOME')
    args = parser.parse_args()
    # Get input files from command line arguments
    genome_file = args.genome
    #masked_file = args.masked
    read1_fastq = args.read1
    read2_fastq = args.read2
    # Test if input files exist
    test_path_existence(genome_file)
    #test_path_existence(masked_file)
    test_path_existence(read1_fastq)
    test_path_existence(read2_fastq)
    # All input tests passed

    # # # Create genome indeces # # #
    if mapper == 'bowtie2':
        # Create bowtie2 index
        bt2_gnm_idx = re.sub(r'\.fa.*$', '', genome_file)
        if not os.path.exists(bt2_gnm_idx + '.rev.2.bt2'):
            subprocess.run(['bowtie2-build', genome_file, bt2_gnm_idx], capture_output=True)
    elif mapper == 'bwa':
        # Create bwa index
        if not os.path.exists(genome_file + '.sa'):
            subprocess.run(['bwa', 'index', genome_file], capture_output=True)
    # Create fasta index
    gnm_fai_file = genome_file + '.fai'
    if not os.path.exists(gnm_fai_file):
        subprocess.run(['samtools', 'faidx', genome_file], capture_output=True)
    # Create genome sizes file
    gnm_len_file = genome_file + '.sizes'
    if not os.path.exists(gnm_len_file):
        os.system(f"cut -f1,2 {gnm_fai_file} > {gnm_len_file}")

    # # # Map fastq files to genome and filter/sort resulting bam files # # #
    bam_fs_files = map_hic_fastq_pair(read1_fastq, read2_fastq, genome_file, cores, mapper)

    # # # Merge read1 and read2 bam files # # #
    merge_bam_file = merge_bam_files(bam_fs_files, gnm_fai_file, cores, map_q)
    
    # # # Remove duplicates from merged bam file # # #
    merge_dr_bam_file = deduplicate_bam_file(merge_bam_file, cores)
    
    # # # Convert merged bam file to bed file # # #
    merge_bed_file = bam_to_bed(merge_dr_bam_file)
    
    # # # Create assembly contact map # # #
    # Convert merged bed file to juicer pre format
    bed_pre_file = bed_to_pre(merge_bed_file)
    # Produce hic file
    hic_file = juicer_tools_hic(bed_pre_file, gnm_len_file)
    # Annotate TADs with given resolution
    bedpe_file_1 = juicer_tools_arrowhead(hic_file, cores, 5000)
    bedpe_file_2 = juicer_tools_arrowhead(hic_file, cores, 10000)
    #print(bedpe_file_1)
    #print(bedpe_file_2)


def test_path_existence(path: str):
    """Test whether given path exists, if not: exit"""
    if not os.path.exists(path):
        print("Error: File or directory '" + path + "' cannot be found! Exit.\n", file=sys.stderr)
        sys.exit()


def test_dependencies(depends: list):
    """Test whether dependencies are installed, if not: exit"""
    for depend in depends:
        result = subprocess.run(['which', depend], capture_output=True)
        if not result.returncode == 0:
            print("Error: Dependency '" + depend + "' cannot be found! Exit.\n", file=sys.stderr)
            sys.exit()


def map_hic_fastq_pair(read1_fq: str, read2_fq: str, genome_fa: str, procs: int, map_prog: str) -> list:
    # List for filtered and sorted bam file names (r1, r2)
    bam_fs_files: list[str] = []
    # Run pipeline for read1 and read2 files
    for fq_file in [read1_fq, read2_fq]:
        # Get the fastq file base name
        fq_base = re.sub(r'\.fq.*$', '', fq_file)
        fq_base = re.sub(r'\.fa.*$', '', fq_base)
        gnm_base = re.sub(r'\.fa.*$', '', genome_fa)
        # Mapper: bowtie2
        if map_prog == 'bowtie2':
            # Name bam file
            bam_file = fq_base + '.' + gnm_base + '.bt2.bam'
            # Get name of bowtie2 index
            bt2_gnm_idx = re.sub(r'\.fa.*$', '', genome_fa)
            # Create bam file if it does not yet exist
            if not os.path.exists(bam_file):
                os.system(f"bowtie2 -p {procs} -x {bt2_gnm_idx} -U {fq_file} "
                          f"| samtools view -Sb "
                          f"| samtools sort > {bam_file}")
        # Mapper: bwa
        elif map_prog == 'bwa':
            # Name bam file
            bam_file = fq_base + '.' + gnm_base + '.bwa.bam'
            # Create bam file if it does not yet exist
            if not os.path.exists(bam_file):
                os.system(f"bwa mem -t {procs} {genome_fa} {fq_file} "
                          f"| samtools view -h -F 256 -Sb "
                          f"| samtools sort > {bam_file}")
        else:
            bam_file = ''
        # Name filtered bam file
        bam_base = re.sub(r'\.bam.*$', '', bam_file)
        bam_filt_file = bam_base + '.filt.bam'
        # Create filtered bam file if it does not yet exist
        if not os.path.exists(bam_filt_file):
            os.system(
                f"samtools view -h {bam_file} | filter_five_end.pl | samtools view -h -Sb - > {bam_filt_file}")
        # Name filtered and sorted bam file
        bam_f_base = re.sub(r'\.bam.*$', '', bam_filt_file)
        bam_fs_file = bam_f_base + '.sort.bam'
        # Create filtered bam file if it does not yet exist
        if not os.path.exists(bam_fs_file):
            os.system(f"samtools sort -n {bam_filt_file} > {bam_fs_file}")
        # Save iltered bam file names
        bam_fs_files.append(bam_fs_file)
    return bam_fs_files


def merge_bam_files(bam_files: list, fai_file: str, procs: int, mapq: int) -> str:
    # Name merged bam file
    bam_f_base = re.sub(r'\.filt\.sort\.bam.*$', '', bam_files[0])
    bam_f_base = re.sub(r'_val_1\.', '.', bam_f_base)
    bam_f_base = re.sub(r'r_1\.', '', bam_f_base)
    bam_f_base = re.sub(r'_1\.', '.', bam_f_base)
    merge_bam_file = bam_f_base + '.fs.merge.bam'
    # Create merged bam file if it does not yet exist
    if not os.path.exists(merge_bam_file):
        os.system(f"two_read_bam_combiner.pl {bam_files[0]} {bam_files[1]} samtools {mapq} "
                  f"| samtools view -bS -t {fai_file} - "
                  f"| samtools sort -@ {procs} -o {merge_bam_file}")
    return merge_bam_file


def deduplicate_bam_file(bam_file: str, procs: int) -> str:
    # Name duplicate-removed merged bam file
    bam_dr_file = re.sub(r'\.bam$', '.dr.bam', bam_file)
    # Create duplicate-removed merged bam file if it does not yet exist
    if not os.path.exists(bam_dr_file):
        # Raise system limit of max open files
        os.system(f"ulimit -n 10240")
        # Remove duplicates with sambamba
        os.system(
            f"sambamba markdup -r -t {procs} --overflow-list-size 600000 --hash-table-size 786432 {bam_file} {bam_dr_file}")
    # Return duplicate-removed merged bam file
    return bam_dr_file


def bam_to_bed(bam_file: str) -> str:
    # Name merged bed file
    bed_file = re.sub(r'bam$', 'bed', bam_file)
    # Create merged bed file if it does not yet exist
    if not os.path.exists(bed_file):
        # Convert bam to bed
        os.system(f"bedtools bamtobed -i {bam_file} > {bed_file}")
        # Sort bed file
        os.system(f"sort -k 4 {bed_file} > tmp && mv tmp {bed_file}")
    # Return bed file
    return bed_file


def bed_to_pre(bed_file: str) -> str:
    # Name pre format file
    bed_pre_file = bed_file + '.pre'
    # Run juicer pre if assembly file does not yet exist
    if not os.path.exists(bed_pre_file):
        # Run juicer pre to produce assembly and merged_dedup files
        os.system(f"bed_to_pre_format_m.pl {bed_file}")
    # Return bed file
    return bed_pre_file


def juicer_tools_hic(bed_pre_file: str, gnm_len_file: str) -> str:
    # Name output file
    hic_file = bed_pre_file + '.hic'
    # Run juicer tools pre if assembly .hic file does not yet exist
    if not os.path.exists(hic_file):
        # Create bash script to run juicer tools pre
        bash_file = open('juicer_tools.sh', "w")
        print(
            f"java -Xmx36G -jar juicer_tools.1.9.9.jar pre {bed_pre_file} {hic_file} {gnm_len_file}", file=bash_file)
        bash_file.close()
        # Run juicer tools pre and remove script
        os.system(f"bash juicer_tools.sh")
        os.system(f"rm juicer_tools.sh")
    # Return bed file
    return hic_file


def juicer_tools_arrowhead(hic_file: str, procs: int, res: int) -> str:
    # Name output file
    out_dir = hic_file + '.tads'
    out_file = out_dir + '/' + str(res) + '_blocks.bedpe'
    # Run juicer tools pre if assembly .hic file does not yet exist
    if not os.path.exists(out_file):
        # Create bash script to run juicer tools pre
        bash_file = open('juicer_tools.sh', "w")
        print(
            f"java -Xmx36G -jar juicer_tools.1.9.9.jar arrowhead --ignore_sparsity --threads {procs} -r {res} {hic_file} {out_dir}", file=bash_file)
        bash_file.close()
        # Run juicer tools pre and remove script
        os.system(f"bash juicer_tools.sh")
        os.system(f"rm juicer_tools.sh")
    # Return bed file
    return out_file


if __name__ == '__main__':
    main()
