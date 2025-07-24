import sys
import csv
import os
import re
import argparse
import logging
import datetime
import requests
import json

from Bio import Align
from Bio.Align import substitution_matrices
import pandas as pd

VERSION = "0.1.0"


def initialise_logging():
    """
    Initialises the logging configuration for the script.

    Creates a log file named with the current date and time, and sets up logging to write
    DEBUG and higher level messages to this file. The log format includes the log level,
    timestamp, and message.

    Returns:
        None
    """
    log = f"check_sequences_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s:%(message)s",
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log)],
    )

    # Disable DEBUG logging for 'requests' and 'urllib3' (used by requests)
    # Set their loggers to a higher level (e.g., INFO, WARNING, ERROR, CRITICAL)
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    # The 'http.client' module is often the one generating the "GET /..." messages
    logging.getLogger("http.client").setLevel(logging.WARNING)


def detect_csv_delimiter(filepath, possible_delimiters=";,\t|", num_lines_to_sniff=5):
    """
    Detects the delimiter of a CSV-like file by sniffing a sample of its content.

    Args:
        filepath (str): The path to the file.
        possible_delimiters (str): A string containing characters that could be delimiters.
                                   The sniffer will test these in order of preference.
                                   Common ones: ';', ',', '\t', '|'.
        num_lines_to_sniff (int): The number of initial lines to read to sample the file
                                  for delimiter detection. More lines provide a better sample,
                                  but for very large files, it's better to keep it small.

    Returns:
        str or None: The detected delimiter character (e.g., ',', ';', '\t', '|')
                     or None if the delimiter cannot be determined or the file cannot be opened.
    """
    try:
        if not os.path.exists(filepath):
            print(f"Error: File not found at '{filepath}'")
            return None

        with open(filepath, "r", newline="") as f:
            sample_data = ""
            for _ in range(num_lines_to_sniff):
                line = f.readline()
                if not line:  # End of file reached
                    break
                sample_data += line

            if (
                not sample_data.strip()
            ):  # Check if sample is empty after stripping whitespace
                print(
                    f"Error: File '{filepath}' is empty or contains only whitespace. Cannot sniff delimiter."
                )
                return None

            sniffer = csv.Sniffer()
            try:
                dialect = sniffer.sniff(sample_data, delimiters=possible_delimiters)
                return dialect.delimiter
            except csv.Error as e:
                print(f"Error: Could not determine CSV delimiter for '{filepath}': {e}")
                return None

    except Exception as e:
        print(f"An unexpected error occurred while detecting delimiter: {e}")
        return None


def get_protein_sequence(gene_symbol, orf_id, organism_name):
    """
    Retrieves the amino acid sequence of the canonical protein variant
    for a given gene symbol and organism from the UniProt REST API.

    Args:
        gene_symbol (str): The gene symbol (e.g., "TP53").
        organism_name (str): The common or scientific name of the organism (e.g., "Human", "Homo sapiens").

    Returns:
        tuple: A tuple containing (sequence_string, uniprot_accession_id) if found,
               otherwise (None, None).
    """
    # Sometimes the annotation may contain gene symbols as zero or missing value, etc
    # Skip these
    if gene_symbol == "0" or gene_symbol == 0 or gene_symbol == "":
        logging.warning(f" Skipping invalid gene symbol: {orf_id} {gene_symbol}")
        return None, None

    # Also check for gene symbols that have been converted to a date by Excel
    # Correct the gene name and log this
    if re.search("^[0-9]+\-Mar$", gene_symbol):
        old_gene_symbol = gene_symbol
        # Extract numerical component from gene_symbol
        num = re.sub("^(\\d+)-Mar$", "\\1", gene_symbol)

        # Strip leading zeros
        num = num.lstrip("0")
        gene_symbol = f"MARCHF{num}"

        logging.info(f" Corrected gene symbol: {old_gene_symbol} to {gene_symbol}")
    elif re.search("^[0-9]+\-Sep$", gene_symbol):
        old_gene_symbol = gene_symbol
        # Extract numerical component from gene_symbol
        num = re.sub("^(\\d+)-Sep$", "\\1", gene_symbol)

        # Strip leading zeros
        num = num.lstrip("0")
        gene_symbol = f"SEPTIN{num}"

        logging.info(f" Corrected gene symbol: {old_gene_symbol} to {gene_symbol}")

    # Base URL for UniProtKB search endpoint
    UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

    # Query parameters:
    # 1. 'query': Search for the gene name AND organism name.
    #    - 'gene:<symbol>' searches the gene name field.
    #    - 'organism_name:<name>' searches by organism name.
    #    - 'reviewed:true' specifically requests Swiss-Prot (canonical/reviewed) entries.
    # 2. 'format': Request JSON format.
    # 3. 'fields': Specify which fields to retrieve (accession, sequence).
    # 4. 'size': Limit to 1 result, as we expect the canonical to be the primary hit.

    params = {
        "query": f"gene:{gene_symbol} AND organism_name:{organism_name} AND reviewed:true",
        "format": "json",
        "fields": "accession,sequence",
        "size": 1,
    }

    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)

        data = response.json()

        if data and data.get("results"):
            # The first result should be the most relevant canonical entry
            entry = data["results"][0]
            uniprot_accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value")

            if uniprot_accession and sequence:
                return sequence, uniprot_accession
            else:
                logging.warning(
                    f" Missing accession or sequence in UniProt response for: {gene_symbol}."
                )
                return None, None
        else:
            logging.warning(
                f" No canonical UniProt entry found for: '{gene_symbol}' in '{organism_name}'."
            )
            return None, None

    except requests.exceptions.HTTPError as http_err:
        logging.error(
            f" HTTP error occurred: {http_err} - {response.text}"
        )  # Print response text for more details
        return None, None
    except requests.exceptions.ConnectionError as conn_err:
        logging.error(f" Connection error occurred: {conn_err}")
        return None, None
    except requests.exceptions.Timeout as timeout_err:
        logging.error(f" Timeout error occurred: {timeout_err}")
        return None, None
    except requests.exceptions.RequestException as req_err:
        logging.error(f" UniProt request error: {req_err}")
        return None, None
    except json.JSONDecodeError as json_err:
        logging.error(f" Decoding error JSON response from UniProt: {json_err}")
        logging.error(
            f" Problematic response content: {response.text[:500]}..."
        )  # Show part of the response
        return None, None
    except Exception as e:
        logging.error(f" An unexpected error occurred: {e}")
        return None, None


def main():

    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="GPSW: A tool for analysing and processing Global Protein Stability Profiling data.",
        prog="gpsw",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=VERSION,
    )

    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="Annotation file containing gene symbols and amino acid sequences to check",
    )

    parser.add_argument(
        "-g",
        "--gene-column",
        type=str,
        help="Column name that contains gene symbols to retrieve the canonical protein sequence for (from annotation file)",
    )

    parser.add_argument(
        "-o",
        "--orf-column",
        type=str,
        help="Column name that contains ORF identifiers to retrieve the canonical protein sequence for (from annotation file)",
    )

    parser.add_argument(
        "-a",
        "--aminoacid-column",
        type=str,
        help="Column name that contains the ORF library amino acid sequences to retrieve the canonical protein sequence for (from annotation file)",
    )

    parser.add_argument(
        "-c",
        "--csv-file",
        type=str,
        help="Gene summary CSV (GPSW output)",
    )

    parser.add_argument(
        "--organism",
        type=str,
        default="Human",
        help="Organism name (e.g., 'Human', 'Mouse')",
    )

    # Parse the command line arguments
    args = parser.parse_args()

    # Load annotation file
    delimiter = detect_csv_delimiter(args.file)
    if not delimiter:
        logging.error(" Could not detect CSV delimiter. Exiting.")
        sys.exit(1)
    annotation = pd.read_csv(args.file, delimiter=delimiter)

    # Get all gene symbols/orf ids from the annotation file
    if args.gene_column not in annotation.columns:
        logging.error(
            f" Gene column '{args.gene_column}' not found in annotation file."
        )
        sys.exit(1)

    gene_symbols = annotation[args.gene_column]

    if args.orf_column not in annotation.columns:
        logging.error(f" ORF column '{args.orf_column}' not found in annotation file.")
        sys.exit(1)

    orf_ids = annotation[args.orf_column]

    # Retrieve the canonical protein sequences for each gene symbol that
    # relates to the ORF id
    canonical_sequences = []
    uniprot_accessions = []
    for gene_symbol, orf_id in zip(gene_symbols, orf_ids):
        # Get the canonical protein sequence for the gene symbol
        sequence, uniprot_accession = get_protein_sequence(
            gene_symbol, orf_id, args.organism
        )
        canonical_sequences.append(sequence)
        uniprot_accessions.append(uniprot_accession)

    # Convert to a DataFrame
    annotation_df = pd.DataFrame(
        {
            "orf_id": orf_ids,
            "canonical_sequence": canonical_sequences,
            "uniprot_accession": uniprot_accessions,
        }
    )

    # Get the amino acid sequences for each ORF id/gene symbol
    if args.aminoacid_column not in annotation.columns:
        logging.error(
            f" Amino acid column '{args.aminoacid_column}' not found in annotation file."
        )
        sys.exit(1)

    orf_amino_acid_sequences = annotation[args.aminoacid_column].tolist()
    orf_amino_acid_sequences = [seq.upper() for seq in orf_amino_acid_sequences]

    # Remove any non-alphabetic characters from the sequences
    orf_amino_acid_sequences = [
        re.sub(r"[^A-Za-z]", "", seq).upper() for seq in orf_amino_acid_sequences
    ]

    # Check if the canonical sequences match the ORF amino acid sequences
    results = []
    for canonical_sequence, amino_acid_sequence in zip(
        annotation_df["canonical_sequence"],
        orf_amino_acid_sequences,
    ):
        if canonical_sequence == amino_acid_sequence:
            results.append(True)
        else:
            results.append(False)

    # Create a DataFrame for the results: orf_id and whether the sequence matches
    results_df = pd.DataFrame(
        {
            "orf_id": annotation_df["orf_id"],
            "uniprot_accession": annotation_df["uniprot_accession"],
            "canonical_sequence": annotation_df["canonical_sequence"],
            "orf_amino_acid_sequence": orf_amino_acid_sequences,
            "match": results,
        }
    )

    # Load GPSW gene summary CSV file
    if not os.path.exists(args.csv_file):
        logging.error(f" Gene summary CSV file '{args.csv_file}' not found.")
        sys.exit(1)
    gene_summary_df = pd.read_csv(args.csv_file)

    # Add to gene summary DataFrame
    gene_summary_df = pd.merge(
        gene_summary_df,
        results_df,
        on="orf_id",
        how="left",
    )

    # Save the results to a new CSV file
    outfile = args.csv_file.replace(".csv", "_annotated.csv")
    gene_summary_df.to_csv(outfile, index=False, na_rep="NA")

    logging.info(f" Results saved to '{outfile}'.")
    print("Done!")

if __name__ == "__main__":
    initialise_logging()  # Initialize logging

    # Get full command line arguments and write to log
    logging.info(f" {' '.join(sys.argv)}")

    # Call the main function to execute the script
    main()
