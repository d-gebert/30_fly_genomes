import argparse
import os
import subprocess

def main():
    """
    Main function to process Hi-C read pairs and generate .pre file.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process Hi-C read pairs BED file.")
    parser.add_argument('bed_file', type=str, help='Path to the Hi-C read pairs BED file')
    args = parser.parse_args()

    bed_file = args.bed_file

    # Open BED file
    with open(bed_file, 'r') as bed, open(f"{bed_file}.pre", 'w') as pre:
        process_bed_file(bed, pre)
    
    # Sort and finalize .pre file
    sort_and_finalize_pre_file(bed_file)

def process_bed_file(bed, pre):
    """
    Processes the BED file to extract read pairs and write to .pre file.
    :param bed: Open BED file handle
    :param pre: Open output .pre file handle
    """
    mq1 = 0
    for line in bed:
        # Get line data
        d = line.strip().split('\t')
        # Get pair id and read id
        pair_id, read_id = d[3].split('/')
        # Get strand
        str_flag = '0' if d[5] == '+' else '1'
        
        if read_id == '1':
            mq1 = d[4]
            pre.write(f"{pair_id} {str_flag} {d[0]} {d[1]} 0 ")
        elif read_id == '2':
            pre.write(f"{str_flag} {d[0]} {d[1]} 1 {mq1} {d[4]}\n")

def sort_and_finalize_pre_file(bed_file):
    """
    Sorts and finalizes the .pre file.
    :param bed_file: Path to the original BED file
    """
    pre_file = f"{bed_file}.pre"
    temp_file = f"{bed_file}.pre.temp"
    
    # Run the awk command
    awk_command = (
        f"awk '{{if ($3 > $7){{ print $1, $6, $7, $8, $9, $11, $2, $3, $4, $5, $10}}"
        f"else {{print}}}}' {pre_file} > {temp_file}"
    )
    subprocess.run(awk_command, shell=True, check=True)
    
    # Sort the .pre.temp file
    sort_command = f"sort -k3,3d -k7,7d {temp_file} > {pre_file}"
    subprocess.run(sort_command, shell=True, check=True)
    
    # Remove the temporary file
    os.remove(temp_file)

if __name__ == "__main__":
    main()
