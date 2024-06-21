import os
import argparse
from pathlib import Path
import shutil
import subprocess
import re

def main():
    """
    Main function to merge outputs from Maker datastore.
    """
    args = parse_arguments()
    input_dirs = find_maker_dirs(args.input)

    if not input_dirs:
        print("No maker output directory found, exiting...")
        return

    if len(input_dirs) == 1:
        print(f"Found maker output directory: {input_dirs[0]}")
    else:
        print(f"Found {len(input_dirs)} maker output directories:")
        for maker_dir in input_dirs:
            print(f"    {maker_dir}")

    output_dir = Path(args.ctlfolder) / args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    for dir_path in input_dirs:
        print(f"\nDealing with {dir_path}:")
        genome_name = Path(dir_path).name.split('.')[0]
        datastore_file = find_datastore_file(dir_path)

        if datastore_file:
            print(f"Datastore file found in {dir_path}, merging annotations...")
        else:
            print(f"Could not find datastore index in {dir_path}, skipping...")
            continue

        outfolder = get_output_folder(args.output, genome_name, len(input_dirs))
        outfolder.mkdir(parents=True, exist_ok=True)

        if any(outfolder.glob("*.fasta")) or any(outfolder.glob("*.gff")):
            print("Output fasta/gff file already exists. Skipping gathering step.")
        else:
            print("Collecting gff and fasta files...")
            file_hds = {}
            collect_recursive(file_hds, datastore_file.parent, outfolder, genome_name, "maker_annotation", "maker_mix")

            for handler in file_hds.values():
                handler.close()

            add_gff_header(outfolder)

        print("Done.")

def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="Merge outputs from Maker datastore.")
    parser.add_argument('-i', '--input', type=str, help='Input directory')
    parser.add_argument('-o', '--output', type=str, help='Output directory', default='annotations')
    parser.add_argument('-c', '--ctlfolder', type=str, help='Control folder', default=os.getcwd())
    return parser.parse_args()

def find_maker_dirs(input_dir):
    """
    Identify Maker output directories.
    """
    if input_dir:
        if not os.path.isdir(input_dir):
            print(f"The input directory {input_dir} doesn't exist.")
            return []
        return [input_dir]
    else:
        current_dir = Path(os.getcwd())
        return [str(entry) for entry in current_dir.iterdir() if entry.is_dir() and 'processed' not in entry.name and entry.name.endswith('.maker.output')]

def find_datastore_file(maker_dir_path):
    """
    Identify the correct datastore index file.
    """
    for file in Path(maker_dir_path).iterdir():
        if file.name.endswith('_master_datastore_index.log'):
            return file
    return None

def get_output_folder(output_base, genome_name, num_dirs):
    """
    Determine the output folder name based on the number of directories and the genome name.
    """
    if num_dirs == 1:
        return Path(output_base)
    else:
        return Path(f"{output_base}_{genome_name}")

def collect_recursive(file_hds, full_path, out, genome_name, maker_annotation_prefix, maker_mix_prefix):
    """
    Recursively collect files from datastore.
    """
    full_path = Path(full_path)
    if full_path.is_file():
        process_file(file_hds, full_path, out, genome_name, maker_annotation_prefix, maker_mix_prefix)
    elif full_path.is_dir():
        if "theVoid" not in full_path.name:
            for sub in full_path.iterdir():
                collect_recursive(file_hds, sub, out, genome_name, maker_annotation_prefix, maker_mix_prefix)

def process_file(file_hds, file_path, out, genome_name, maker_annotation_prefix, maker_mix_prefix):
    """
    Process individual file based on its suffix and type.
    """
    suffix = file_path.suffix
    name = file_path.stem

    if suffix == ".fasta":
        process_fasta(file_hds, file_path, name, out, genome_name, maker_annotation_prefix)
    elif suffix == ".gff":
        process_gff(file_path, out, maker_annotation_prefix, maker_mix_prefix)

def process_fasta(file_hds, file_path, name, out, genome_name, maker_annotation_prefix):
    """
    Process FASTA files and append to appropriate output file.
    """
    key, file_type = None, None

    if re.search(r'\.transcripts', name):
        key, file_type = re.search(r'([^\.]+)\.transcripts', name).groups()[0], "transcripts"
    elif re.search(r'\.proteins', name):
        key, file_type = re.search(r'([^\.]+)\.proteins', name).groups()[0], "proteins"
    elif re.search(r'\.noncoding', name):
        key, file_type = re.search(r'([^\.]+)\.noncoding', name).groups()[0], "noncoding"

    if key:
        output_name = f"{maker_annotation_prefix}.{file_type}.fasta" if key == "maker" else f"{genome_name}.all.maker.{key}.{file_type}.fasta"
        if output_name not in file_hds:
            file_hds[output_name] = open(out / output_name, 'w')
        with open(file_path, 'r') as fh:
            shutil.copyfileobj(fh, file_hds[output_name])

def process_gff(file_path, out, maker_annotation_prefix, maker_mix_prefix):
    """
    Process GFF files and append to appropriate output file.
    """
    subprocess.run(["gawk", "-F", "\t", f'NF==9{{print $0 >> "{out}/{maker_mix_prefix}.gff"}}', file_path])
    subprocess.run(["gawk", '{if($2 ~ /[a-zA-Z]+/) if($2=="maker") { print $0 >> "' + str(out) + '/' + maker_annotation_prefix + '.gff" } else { OFS="\t"; gsub(/:/, "_", $2); print $0 >> "' + str(out) + '/" $2 ".gff" } }', file_path])

def add_gff_header(outfolder):
    """
    Add GFF version header to all GFF files in the output folder.
    """
    for gff_file in outfolder.glob("*.gff"):
        subprocess.run(["sed", "-i", "", "1s/^/##gff-version 3\\\n/", gff_file])

if __name__ == "__main__":
    main()
