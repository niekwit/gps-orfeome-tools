"""
Convert ORF amino acid sequences from CSV to FASTA format
"""

import os
import argparse
import pandas as pd


def create_fasta(df, csv, stability, args):
    """
    Create a FASTA file from a DataFrame containing sequence names and sequences.
    """
    # Merge the filtered GPSW DataFrame with the CSV DataFrame
    df = pd.merge(df, csv, on=["orf_id", "gene"], how="inner")

    # Create a new column 'name' combining gene and ORF names
    df["name"] = df["gene"] + "_" + df["orf_id"]

    # Only keep name and sequence columns
    df = df[["name", "sequence"]]

    # Define the output file name based on stability
    base_name = os.path.basename(args.gpsw).replace("_gene.summary.csv", "")
    output_file = os.path.join(args.outdir, f"{base_name}_{stability}.fasta")

    # Prepend > to each sequence name
    df["name"] = ">" + df["name"]

    # Remove any non-alphanumeric characters from the sequences
    df["sequence"] = df["sequence"].str.replace(r"[^A-Za-z0-9]", "", regex=True)

    # Write to FASTA file
    df.to_csv(output_file, sep="\n", index=False, header=False)


def main(args):

    # Open CSV and GPSW files
    csv = pd.read_csv(args.input)
    gpsw = pd.read_csv(args.gpsw)

    # Rename CSV columns for consistency
    csv.rename(
        columns={
            args.gene_column: "gene",
            args.orf_column: "orf_id",
            args.sequence_column: "sequence",
        },
        inplace=True,
    )

    # Keep only relevant columns in csv
    csv = csv[["orf_id", "gene", "sequence"]]

    ## Prepare data for FASTA files based on stability conditions
    _list = ["stabilised", "destabilised"]

    for stability in _list:
        # Filter the GPSW DataFrame based on the stability condition
        gpsw_filtered = gpsw[gpsw[stability] == True]

        # Create FASTA file for the current stability condition
        create_fasta(gpsw_filtered, csv, stability, args)

    ## Prepare data for background sequences
    # Absolute dPSI should be less than the background cutoff
    df = gpsw[abs(gpsw["delta_PSI_mean"]) < args.background_dpsi_cutoff]

    # Create FASTA file for background sequences
    create_fasta(df, csv, "background", args)

    print("Done!")


if __name__ == "__main__":
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
        required=True,
        help="Name of the column containing gene names.",
    )

    parser.add_argument(
        "-o",
        "--orf-column",
        required=True,
        help="Name of the column containing ORF names.",
    )

    parser.add_argument(
        "-s",
        "--sequence-column",
        required=True,
        type=str,
        help="Name of the column containing amino acid sequences.",
    )

    parser.add_argument(
        "--gpsw",
        required=True,
        type=str,
        help="GPSW gene summary file",
    )

    parser.add_argument(
        "-d",
        "--dpsi-cutoff",
        type=float,
        default=1.0,
        help="dPSI cutoff value for filtering sequences",
    )

    parser.add_argument(
        "-b",
        "--background-dpsi-cutoff",
        type=float,
        default=0.1,
        help="dPSI cutoff value for filtering background sequences",
    )

    parser.add_argument(
        "--outdir",
        default=".",
        type=str,
        help="Output directory for FASTA files (default: current directory)",
    )

    args = parser.parse_args()

    main(args)
