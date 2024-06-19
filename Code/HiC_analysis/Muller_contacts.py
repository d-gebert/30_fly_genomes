import sys
import math

# Global constants
MAPQ_THRESHOLD = 10  # Minimum mapping quality threshold for reads to be considered
MULLERS = ['A', 'B', 'C', 'D', 'E', 'F']  # List of Muller elements

def main():
    """
    Main function to coordinate the processing of input files and calculation of results.
    """
    # Handle command line arguments
    if len(sys.argv) < 3:
        print_usage_and_exit(sys.argv[0])
	# Input file names
    bed_file = sys.argv[1]
    fai_file = sys.argv[2]

    # Extract chromosome lengths
    chr_lens = get_fai_seq_lens(fai_file)
    
    # Parse BED file to get read pairs
    read_pairs = parse_bed_file(bed_file)
    
    # Count contacts per Muller combination
    muller_contacts = count_muller_contacts(read_pairs)
    
    # Calculate and print the results
    calculate_and_print_results(muller_contacts, chr_lens)

def get_fai_seq_lens(filename: str) -> dict:
    """
    Parses a .fai file to get chromosome lengths.
    :param filename: Path to the .fai file
    :return: Dictionary of chromosome lengths
    """
    chr_lens = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chr_lens[parts[0]] = int(parts[1])
    return chr_lens

def parse_bed_file(filename: str) -> dict:
    """
    Parses BED file data to extract read pairs.
    :param filename: Path to the BED file
    :return: Dictionary of read pairs organized by pair ID
    """
    read_pairs = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()  # Remove any trailing whitespace
            dat = line.split('\t')  # Split the line into columns based on tab delimiter
            chr, rid, mpq = dat[0], dat[3], int(dat[4])  # Extract chromosome, read ID, and mapping quality
            pid = rid.rsplit('/', 1)[0]  # Extract pair ID from read ID
            num = int(rid.rsplit('/', 1)[1])  # Extract read number from read ID
            if pid not in read_pairs:
                read_pairs[pid] = {}
            read_pairs[pid][num] = (chr, mpq)  # Store chromosome and mapping quality in the dictionary
    return read_pairs

def count_muller_contacts(read_pairs: dict) -> dict:
    """
    Counts contacts per Muller combination based on read pairs.
    :param read_pairs: Dictionary of read pairs
    :return: Dictionary of Muller contacts
    """
    muller_contacts = {}
    for pid, reads in read_pairs.items():
        # Skip read pairs with mapping quality below the threshold
        if reads[1][1] < MAPQ_THRESHOLD or reads[2][1] < MAPQ_THRESHOLD:
            continue
        chr1, chr2 = reads[1][0], reads[2][0]  # Get chromosomes for the read pair
        if chr1 not in muller_contacts:
            muller_contacts[chr1] = {}
        if chr2 not in muller_contacts[chr1]:
            muller_contacts[chr1][chr2] = 0
        muller_contacts[chr1][chr2] += 1  # Increment contact count for the chromosome pair
    return muller_contacts

def calculate_and_print_results(muller_contacts: dict, chr_lens: dict):
    """
    Calculates and prints the normalized contact values for each Muller combination.
    :param muller_contacts: Dictionary of Muller contacts
    :param chr_lens: Dictionary of chromosome lengths
    """
    for mul_x in MULLERS:
        for mul_y in MULLERS:
            lens = sorted([chr_lens[mul_x], chr_lens[mul_y]])  # Get and sort the chromosome lengths
            len_val = math.sqrt(lens[0] * lens[1])  # Calculate the geometric mean of the lengths
            contact_count = muller_contacts.get(mul_x, {}).get(mul_y, 0)  # Get contact count or 0 if not present
            val = round(contact_count / len_val * 1000, 2)  # Normalize the contact count and round it
            print(f"{val}\t", end="")  # Print the value followed by a tab
        print("")  # Print a newline after each row

def print_usage_and_exit(program_name: str):
    """
    Prints the usage message and exits the program.
    :param program_name: Name of the program (usually sys.argv[0])
    """
    usage = f"Usage: python {program_name} <merged_hic.bed> <genome.fai>\n"
    sys.exit(usage)

if __name__ == "__main__":
    main()
