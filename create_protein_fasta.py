"""
Convert ORF amino acid sequences from CSV to FASTA format
"""

import argparse
import os
import csv


# Set up the argument parser
parser = argparse.ArgumentParser(
    description="Convert ORF amino acid sequences from CSV to FASTA format."
)

parser.add_argument(
    "-i",
    "--input",
    required=True,
    type=str,
    help="Path to the input CSV file containing ORF sequences.",
)

parser.add_argument(
    "-g",
    "--gene-column",
    default=None,
    help="Name of the column containing gene names.",
)

parser.add_argument(
    "-o",
    "--orf-column",
    default=None,
    help="Name of the column containing ORF names.",
)

parser.add_argument(
    "-s",
    "--sequence-column",
    required=True,
    type=str,
    help="Name of the column containing amino acid sequences.",
)

args = parser.parse_args()

# Check if the input file exists
if not os.path.isfile(args.input):
    raise FileNotFoundError(f"Input file {args.input} does not exist.")

# Read the CSV file and write to FASTA format
input_file = args.input.replace(".csv", ".fasta")
print(f"Writing FASTA file to {input_file}")

gene_column = args.gene_column
orf_column = args.orf_column

with open(args.input, "r") as csvfile, open(input_file, "w") as fasta_file:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if gene_column and not orf_column:
            name = row[gene_column]
        elif not gene_column and orf_column:
            name = row[orf_column]
        elif gene_column and orf_column:
            name = f"{row[gene_column]}_{row[orf_column]}"
        else:
            raise ValueError("Either gene and/or ORF columns must be specified.")

        sequence = row[args.sequence_column]

        # Write the sequence in FASTA format
        fasta_file.write(f">{name}\n")
        fasta_file.write(f"{sequence}\n")

print("Done!")
