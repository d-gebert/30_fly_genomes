import argparse
import re
from pathlib import Path

def main():
    """
    Main function to rename gene IDs in a GFF file based on an ortholog table.
    """
    args = parse_arguments()

    gff_file = args.gff_file
    tsv_file = args.tsv_file
    out_file = f"{gff_file}.rn.gff"

    # Get a dictionary with gene IDs and names
    gid_dict = tsv_to_dict(tsv_file)

    gene2orth = {}
    orth2gene = {}

    # Modify gene orthologs dictionary
    for gid in sorted(gid_dict.keys()):
        orthologs = gid_dict[gid]
        n_orths = len(orthologs)
        n_gcids = sum(1 for orth in orthologs if orth.startswith("CG"))

        # Remove CG IDs if there are other IDs available
        if n_orths > n_gcids:
            gid_dict[gid] = [orth for orth in orthologs if not orth.startswith("CG")]

        # Set ortholog name
        orth_name = "/".join(gid_dict[gid])
        orth2gene.setdefault(orth_name, {})[gid] = 1
        gene2orth[gid] = orth_name

    # Handle cases where an ortholog is assigned to more than one gene ID
    for ort in sorted(orth2gene.keys()):
        if len(orth2gene[ort]) > 1:
            for i, gid in enumerate(sorted(orth2gene[ort].keys()), 1):
                gene2orth[gid] = f"{ort}.{i}"

    # Open GFF file and output file
    with open(gff_file, 'r') as gff, open(out_file, 'w') as out:
        gid = ''
        for line in gff:
            line = line.rstrip()  # better chomp

            # Skip comment lines
            if line.startswith("#"):
                continue

            # Get line data
            fields = line.split("\t")
            typ = fields[2]
            inf = fields[8]

            # Case: gene type
            if typ == 'gene':
                match = re.search(r"Name=(\S+)", inf)
                if match:
                    gid = match.group(1)
                    if gid in gene2orth:
                        print(f"{gid}\t{gene2orth[gid]}")
                    else:
                        gid = ''
                else:
                    gid = ''

            if gid:
                line = line.replace(gid, gene2orth[gid])

            out.write(line + "\n")

def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="Rename gene IDs in a GFF file based on an OrthoFinder ortholog table.")
    parser.add_argument("gff_file", type=str, help="Input GFF file (maker_annotation.gff)")
    parser.add_argument("tsv_file", type=str, help="Input ortholog table TSV file (maker_annotation.proteins.tsv)")
    return parser.parse_args()

def tsv_to_dict(filename):
    """
    Parse a TSV file to create a dictionary with gene IDs and their orthologs.
    :param filename: Path to the TSV file
    :return: Dictionary with gene IDs as keys and list of orthologs as values
    """
    data_dict = {}
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.rstrip()  # better chomp

            # Skip title line
            if line.startswith("Orthogroup"):
                continue

            # Split line
            fields = line.split("\t")

            # Get target fields
            ort_group = fields[0]
            maker_ids = fields[2]
            ort_prots = fields[3]

            # Split maker IDs list and remove mRNA tags
            maker_ids = [id.split("-mRNA")[0] for id in maker_ids.split(", ")]

            # Split ortholog proteins list and remove protein tags to get (unique) gene names
            ort_genes = set()
            for ort_prot in ort_prots.split(", "):
                gene_name = ort_prot.split("-P")[0]
                ort_genes.add(gene_name)

            # Link each maker ID to ortholog genes
            for id in maker_ids:
                if id not in data_dict:
                    data_dict[id] = []
                data_dict[id].extend(ort_genes)

    return data_dict

if __name__ == "__main__":
    main()
