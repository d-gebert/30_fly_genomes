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
    """
    Main function to process Hi-C data, align reads, and generate contact maps.
    """
    test_dependencies(dependencies)
    
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-1", "--read1", type=str, help="read1 fastq file", required=True, metavar='READ1')
    required.add_argument("-2", "--read2", type=str, help="read2 fastq file", required=True, metavar='READ2')
    required.add_argument("-g", "--genome", type=str, help="genome fasta file", required=True, metavar='GENOME')
    args = parser.parse_args()
    
    genome_file = args.genome
    read1_fastq = args.read1
    read2_fastq = args.read2
    
    test_path_existence(genome_file)
    test_path_existence(read1_fastq)
    test_path_existence(read2_fastq)
    
    if mapper == 'bowtie2':
        bt2_gnm_idx = re.sub(r'\.fa.*$', '', genome_file)
        if not os.path.exists(bt2_gnm_idx + '.rev.2.bt2'):
            subprocess.run(['bowtie2-build', genome_file, bt2_gnm_idx], capture_output=True)
    elif mapper == 'bwa':
        if not os.path.exists(genome_file + '.sa'):
            subprocess.run(['bwa', 'index', genome_file], capture_output=True)
    
    gnm_fai_file = genome_file + '.fai'
    if not os.path.exists(gnm_fai_file):
        subprocess.run(['samtools', 'faidx', genome_file], capture_output=True)
    
    gnm_len_file = genome_file + '.sizes'
    if not os.path.exists(gnm_len_file):
        os.system(f"cut -f1,2 {gnm_fai_file} > {gnm_len_file}")

    bam_fs_files = map_hic_fastq_pair(read1_fastq, read2_fastq, genome_file, cores, mapper)
    merge_bam_file = merge_bam_files(bam_fs_files, gnm_fai_file, cores, map_q)
    merge_dr_bam_file = deduplicate_bam_file(merge_bam_file, cores)
    merge_bed_file = bam_to_bed(merge_dr_bam_file)
    
    bed_pre_file = bed_to_pre(merge_bed_file)
    hic_file = juicer_tools_hic(bed_pre_file, gnm_len_file)
    
    print(f"HiC contact map file created: '{hic_file}'\n")

def test_path_existence(path: str):
    """
    Test whether given path exists, if not: exit.
    :param path: Path to be tested
    """
    if not os.path.exists(path):
        print(f"Error: File or directory '{path}' cannot be found! Exit.\n", file=sys.stderr)
        sys.exit()

def test_dependencies(depends: list):
    """
    Test whether dependencies are installed, if not: exit.
    :param depends: List of dependencies to be tested
    """
    for depend in depends:
        result = subprocess.run(['which', depend], capture_output=True)
        if result.returncode != 0:
            print(f"Error: Dependency '{depend}' cannot be found! Exit.\n", file=sys.stderr)
            sys.exit()

def map_hic_fastq_pair(read1_fq: str, read2_fq: str, genome_fa: str, procs: int, map_prog: str) -> list:
    """
    Map Hi-C fastq pair to the genome and filter/sort resulting bam files.
    :param read1_fq: Path to read1 fastq file
    :param read2_fq: Path to read2 fastq file
    :param genome_fa: Path to genome fasta file
    :param procs: Number of processor cores to use
    :param map_prog: Mapping program to use ('bowtie2' or 'bwa')
    :return: List of filtered and sorted bam file names
    """
    bam_fs_files = []
    for fq_file in [read1_fq, read2_fq]:
        fq_base = re.sub(r'\.fq.*$', '', fq_file)
        fq_base = re.sub(r'\.fa.*$', '', fq_base)
        gnm_base = re.sub(r'\.fa.*$', '', genome_fa)
        
        if map_prog == 'bowtie2':
            bam_file = fq_base + '.' + gnm_base + '.bt2.bam'
            bt2_gnm_idx = re.sub(r'\.fa.*$', '', genome_fa)
            if not os.path.exists(bam_file):
                os.system(f"bowtie2 -p {procs} -x {bt2_gnm_idx} -U {fq_file} "
                          f"| samtools view -Sb "
                          f"| samtools sort > {bam_file}")
        elif map_prog == 'bwa':
            bam_file = fq_base + '.' + gnm_base + '.bwa.bam'
            if not os.path.exists(bam_file):
                os.system(f"bwa mem -t {procs} {genome_fa} {fq_file} "
                          f"| samtools view -h -F 256 -Sb "
                          f"| samtools sort > {bam_file}")
        else:
            bam_file = ''
        
        bam_base = re.sub(r'\.bam.*$', '', bam_file)
        bam_filt_file = bam_base + '.filt.bam'
        if not os.path.exists(bam_filt_file):
            os.system(
                f"samtools view -h {bam_file} | filter_five_end.pl | samtools view -h -Sb - > {bam_filt_file}")
        
        bam_f_base = re.sub(r'\.bam.*$', '', bam_filt_file)
        bam_fs_file = bam_f_base + '.sort.bam'
        if not os.path.exists(bam_fs_file):
            os.system(f"samtools sort -n {bam_filt_file} > {bam_fs_file}")
        
        bam_fs_files.append(bam_fs_file)
    return bam_fs_files

def merge_bam_files(bam_files: list, fai_file: str, procs: int, mapq: int) -> str:
    """
    Merge read1 and read2 bam files.
    :param bam_files: List of bam file paths
    :param fai_file: Path to the FAI file
    :param procs: Number of processor cores to use
    :param mapq: Mapping quality threshold
    :return: Path to the merged bam file
    """
    bam_f_base = re.sub(r'\.filt\.sort\.bam.*$', '', bam_files[0])
    bam_f_base = re.sub(r'_val_1\.', '.', bam_f_base)
    bam_f_base = re.sub(r'r_1\.', '', bam_f_base)
    bam_f_base = re.sub(r'_1\.', '.', bam_f_base)
    merge_bam_file = bam_f_base + '.fs.merge.bam'
    
    if not os.path.exists(merge_bam_file):
        os.system(f"two_read_bam_combiner.pl {bam_files[0]} {bam_files[1]} samtools {mapq} "
                  f"| samtools view -bS -t {fai_file} - "
                  f"| samtools sort -@ {procs} -o {merge_bam_file}")
    return merge_bam_file

def deduplicate_bam_file(bam_file: str, procs: int) -> str:
    """
    Remove duplicates from merged bam file.
    :param bam_file: Path to the bam file
    :param procs: Number of processor cores to use
    :return: Path to the duplicate-removed bam file
    """
    bam_dr_file = re.sub(r'\.bam$', '.dr.bam', bam_file)
    if not os.path.exists(bam_dr_file):
        os.system("ulimit -n 10240")
        os.system(
            f"sambamba markdup -r -t {procs} --overflow-list-size 600000 --hash-table-size 786432 {bam_file} {bam_dr_file}")
    return bam_dr_file

def bam_to_bed(bam_file: str) -> str:
    """
    Convert bam file to bed file.
    :param bam_file: Path to the bam file
    :return: Path to the bed file
    """
    bed_file = re.sub(r'bam$', 'bed', bam_file)
    if not os.path.exists(bed_file):
        os.system(f"bedtools bamtobed -i {bam_file} > {bed_file}")
        os.system(f"sort -k 4 {bed_file} > tmp && mv tmp {bed_file}")
    return bed_file

def bed_to_pre(bed_file: str) -> str:
    """
    Convert bed file to juicer pre format.
    :param bed_file: Path to the bed file
    :return: Path to the pre format file
    """
    bed_pre_file = bed_file + '.pre'
    if not os.path.exists(bed_pre_file):
        os.system(f"BED_to_PRE_format.py {bed_file}")
    return bed_pre_file

def juicer_tools_hic(bed_pre_file: str, gnm_len_file: str) -> str:
    """
    Produce hic file from pre format bed file.
    :param bed_pre_file: Path to the pre format bed file
    :param gnm_len_file: Path to the genome lengths file
    :return: Path to the hic file
    """
    hic_file = bed_pre_file + '.hic'
    if not os.path.exists(hic_file):
        with open('juicer_tools.sh', "w") as bash_file:
            print(
                f"java -Xmx36G -jar juicer_tools.1.9.9.jar pre {bed_pre_file} {hic_file} {gnm_len_file}", file=bash_file)
        os.system("bash juicer_tools.sh")
        os.system("rm juicer_tools.sh")
    return hic_file

if __name__ == '__main__':
    main()
