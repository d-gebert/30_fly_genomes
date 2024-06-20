import argparse
import pandas as pd # type: ignore
import numpy as np # type: ignore
from collections import defaultdict
from sklearn.decomposition import PCA # type: ignore
from sklearn.preprocessing import StandardScaler # type: ignore
from sklearn.impute import SimpleImputer # type: ignore

# Define the size of bins and minimum contig length
DEFAULT_BIN_SIZE = 100000
DEFAULT_MIN_CONTIG_LENGTH = 100000

def main():
    """
    Process a BED file and create a read pair count matrix.
    """
    parser = argparse.ArgumentParser(description="Process BED file to create a read pair count matrix.")
    parser.add_argument('file_path', type=str, help='Path to the BED file')
    parser.add_argument('fai_path', type=str, help='Path to the FAI file')
    parser.add_argument('-b', '--bin_size', type=int, default=DEFAULT_BIN_SIZE, help='Size of the bins (default: 100000)')
    parser.add_argument('-m', '--min_contig_length', type=int, default=DEFAULT_MIN_CONTIG_LENGTH, help='Minimum contig length to be considered (default: 100000)')
    parser.add_argument('-t', '--transformation', type=str, choices=['log', 'sqrt', 'none'], default='log', help='Transformation to apply to the matrix (default: log)')
    parser.add_argument('-p', '--pca_stage', type=str, choices=['before', 'after'], default='before', help='Stage to perform PCA (default: before)')
    args = parser.parse_args()

    # Ensure min_contig_length is at least as high as bin_size
    if args.min_contig_length < args.bin_size:
        args.min_contig_length = args.bin_size

    print("Processing FAI file...")
    contig_lengths = process_fai_file(args.fai_path, args.min_contig_length)
    print(f"Contig lengths: {contig_lengths}")

    print("Processing BED file...")
    matrix, bin_names = process_bed_file(args.file_path, contig_lengths, args.bin_size)
    print(f"Matrix shape: {matrix.shape}")

    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("The processed matrix is empty. Please check the input data.")

    print("Performing ICE normalization...")
    balanced_matrix = ice_normalization(matrix)
    print(f"Balanced matrix shape: {balanced_matrix.shape}")

    # Apply transformation if specified
    if args.transformation == 'log':
        print("Applying log transformation...")
        transformed_matrix = np.log10(balanced_matrix + 1)
    elif args.transformation == 'sqrt':
        print("Applying square root transformation...")
        transformed_matrix = np.sqrt(balanced_matrix)
    else:
        transformed_matrix = balanced_matrix

    # Round the transformed matrix to one decimal place
    transformed_matrix = np.round(transformed_matrix, 1)
    print("Transformation applied.")

    # Perform PCA at the specified stage
    if args.pca_stage == 'before':
        print("Performing PCA before transformation...")
        pc_values = perform_pca(balanced_matrix)
    else:
        print("Performing PCA after transformation...")
        pc_values = perform_pca(transformed_matrix)

    print("PCA completed.")

    # Save the matrix and PCA values to CSV files
    print("Saving matrix and PCA values...")
    save_matrix(transformed_matrix, bin_names)
    save_pca_values(pc_values, bin_names)
    print("Files saved.")

def get_bin(position: int, bin_size: int) -> int:
    """
    Returns the bin index for a given position.
    :param position: The position in the genome
    :param bin_size: The size of each bin
    :return: The bin index
    """
    return position // bin_size

def process_fai_file(fai_path: str, min_contig_length: int) -> dict:
    """
    Processes the FAI file to get contigs/chromosome lengths.
    :param fai_path: Path to the FAI file
    :param min_contig_length: Minimum contig length to be considered
    :return: Dictionary with contig names as keys and their lengths as values
    """
    contig_lengths = {}
    with open(fai_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            contig_name = parts[0]
            contig_length = int(parts[1])
            if contig_length >= min_contig_length:
                contig_lengths[contig_name] = contig_length
    return contig_lengths

def process_bed_file(file_path: str, contig_lengths: dict, bin_size: int) -> tuple:
    """
    Processes the BED file to create a read pair count matrix.
    :param file_path: Path to the BED file
    :param contig_lengths: Dictionary with contig names as keys and their lengths as values
    :param bin_size: The size of each bin
    :return: A tuple containing the matrix of read pair counts and a list of bin names
    """
    counts = defaultdict(list)
    chromosome_positions = defaultdict(list)
    chromosome_bins = {}

    with open(file_path, 'r') as file:
        for line in file:
            chrom, start, end, read_id, _, strand = line.strip().split()
            if chrom not in contig_lengths:
                continue
            start, end = int(start), int(end)
            read_pair = read_id.split('/')[0]
            bin_start = get_bin(start, bin_size)
            bin_end = get_bin(end, bin_size)
            
            counts[read_pair].append((chrom, bin_start, bin_end))
            chromosome_positions[chrom].append(bin_start)
            chromosome_positions[chrom].append(bin_end)
    
    # Determine bin range for each chromosome
    current_bin_index = 0
    bin_names = []
    for chrom in sorted(chromosome_positions):
        unique_bins = sorted(set(chromosome_positions[chrom]))
        chromosome_bins[chrom] = {bin_idx: current_bin_index + i for i, bin_idx in enumerate(unique_bins)}
        bin_names.extend([f"{chrom}:{bin_idx*bin_size}-{(bin_idx+1)*bin_size}" for bin_idx in unique_bins])
        current_bin_index += len(unique_bins)
    
    # Create a full matrix for the genome
    total_bins = current_bin_index
    matrix = np.zeros((total_bins, total_bins), dtype=int)

    for read_pair, locations in counts.items():
        if len(locations) == 2:
            (chrom1, bin_start1, bin_end1), (chrom2, bin_start2, bin_end2) = locations
            bin_idx1_start = chromosome_bins[chrom1][bin_start1]
            bin_idx1_end = chromosome_bins[chrom1][bin_end1]
            bin_idx2_start = chromosome_bins[chrom2][bin_start2]
            bin_idx2_end = chromosome_bins[chrom2][bin_end2]

            matrix[bin_idx1_start, bin_idx2_start] += 1
            matrix[bin_idx2_start, bin_idx1_start] += 1  # Assuming symmetry

            if bin_idx1_start != bin_idx1_end:
                matrix[bin_idx1_end, bin_idx2_start] += 1
                matrix[bin_idx2_start, bin_idx1_end] += 1  # Assuming symmetry

            if bin_idx2_start != bin_idx2_end:
                matrix[bin_idx1_start, bin_idx2_end] += 1
                matrix[bin_idx2_end, bin_idx1_start] += 1  # Assuming symmetry

            if bin_idx1_start != bin_idx1_end and bin_idx2_start != bin_idx2_end:
                matrix[bin_idx1_end, bin_idx2_end] += 1
                matrix[bin_idx2_end, bin_idx1_end] += 1  # Assuming symmetry

    return matrix, bin_names

def ice_normalization(matrix: np.ndarray, max_iter: int = 100, epsilon: float = 1e-5) -> np.ndarray:
    """
    Performs ICE normalization on the Hi-C contact matrix.
    :param matrix: The contact matrix to be normalized
    :param max_iter: Maximum number of iterations
    :param epsilon: Convergence threshold
    :return: The normalized (balanced) matrix
    """
    bias = np.ones(matrix.shape[0])
    for iteration in range(max_iter):
        bias_new = np.divide(1, np.nanmean(matrix * bias, axis=0), where=np.nanmean(matrix * bias, axis=0) != 0)
        bias_new /= np.nanmean(bias_new, where=~np.isnan(bias_new))  # Normalize bias to have mean 1

        if np.nanmax(np.abs(bias - bias_new), initial=0, where=~np.isnan(np.abs(bias - bias_new))) < epsilon:
            break
        bias = bias_new
    
    matrix_balanced = matrix * bias[:, np.newaxis] * bias[np.newaxis, :]
    return matrix_balanced

def perform_pca(matrix: np.ndarray, n_components: int = 3) -> np.ndarray:
    """
    Performs PCA on the Hi-C contact matrix.
    :param matrix: The matrix to perform PCA on
    :param n_components: Number of principal components to return
    :return: The first principal component values
    """
    imputer = SimpleImputer(strategy='mean')
    matrix_imputed = imputer.fit_transform(matrix)
    
    scaler = StandardScaler()
    z_matrix = scaler.fit_transform(matrix_imputed)
    corr_matrix = pd.DataFrame(z_matrix).corr()

    print("Correlation matrix details:")
    print(f"Shape: {corr_matrix.shape}")
    print(f"Any NaNs: {corr_matrix.isna().any().any()}")
    print(f"Number of NaNs: {corr_matrix.isna().sum().sum()}")

    corr_matrix = corr_matrix.fillna(0)

    pca = PCA(n_components=n_components)
    matrix_pca = pca.fit_transform(corr_matrix)
    return matrix_pca[:, 0]  # Return the first principal component values

def save_matrix(matrix: np.ndarray, bin_names: list) -> None:
    """
    Saves the matrix to a CSV file with bin names.
    :param matrix: The matrix to be saved
    :param bin_names: List of bin names corresponding to the rows and columns of the matrix
    """
    df = pd.DataFrame(matrix, index=bin_names, columns=bin_names)
    df.to_csv("genome_matrix_balanced_log10.csv", index=True, header=True)

def save_pca_values(pc_values: np.ndarray, bin_names: list) -> None:
    """
    Saves the PCA values to a CSV file with bin names.
    :param pc_values: The PCA values to be saved
    :param bin_names: List of bin names corresponding to the PCA values
    """
    df = pd.DataFrame(pc_values, index=bin_names, columns=["PC1"])
    df.to_csv("pca_values.csv", index=True, header=True)

if __name__ == "__main__":
    main()
