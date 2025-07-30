import sys
import csv
import os
import re
import argparse
import logging
import datetime
import requests
import json
import itertools
import collections
import xml.etree.ElementTree as ET

from tqdm import tqdm
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices


# Custom filter to allow specific levels to pass to the console handler
class ConsoleLevelFilter(logging.Filter):
    def filter(self, record):
        # Allow INFO, ERROR, and CRITICAL to pass
        return record.levelno in [logging.INFO, logging.ERROR, logging.CRITICAL]


# TQDM logging handler
# This ensures that log messages printed during a tqdm loop
# do not mess up the progress bar display.
class TqdmLoggingHandler(logging.StreamHandler):
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(
                msg, file=sys.stdout
            )  # Use tqdm's write to not interfere with progress bar
            self.flush()
        except (IOError, BrokenPipeError):
            # Handle situations where the pipe to stdout might be broken
            pass
        except Exception:
            self.handleError(record)


def initialise_logging():
    """
    Initialises the logging configuration for the script.

    Sets up logging to write:
    - All messages (DEBUG and higher) to a time-stamped log file.
    - INFO, ERROR, and CRITICAL messages to the console (stdout).
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # CRITICAL FIX: Clear all existing handlers to prevent duplicates
    if logger.handlers:
        for handler in list(
            logger.handlers
        ):  # Use list() to avoid issues with modifying the list while iterating
            logger.removeHandler(handler)

    log_file_name = (
        f"check_sequences_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
    )

    # 1. File Handler: for all messages (DEBUG and up)
    file_handler = logging.FileHandler(log_file_name)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        "%(levelname)s:%(asctime)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # 2. TQDM Logging Handler: for console output (Filtered to INFO, ERROR, CRITICAL)
    tqdm_handler = TqdmLoggingHandler()
    tqdm_handler.setLevel(logging.INFO)
    tqdm_formatter = logging.Formatter("%(message)s")
    tqdm_handler.setFormatter(tqdm_formatter)
    tqdm_handler.addFilter(ConsoleLevelFilter())
    logger.addHandler(tqdm_handler)

    # Disable verbose DEBUG logging from external libraries
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("http.client").setLevel(logging.WARNING)
    logging.getLogger("chardet").setLevel(logging.WARNING)
    logging.getLogger("charset_normalizer").setLevel(logging.WARNING)


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


def get_all_uniprot_isoforms(gene_symbol, orf_id, organism_name, max_isoforms=5):
    """
    Retrieve up to `max_isoforms` UniProt isoform sequences for a given gene symbol and organism.
    Sequences are sorted by isoform ID (e.g., Q9UH92-1, Q9UH92-2, ...).
    """
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"gene:{gene_symbol} AND organism_name:{organism_name} AND reviewed:true"

    params = {
        "query": query,
        "format": "json",
        "fields": "accession",
        "size": 1,
    }

    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()

        if not data.get("results"):
            logging.warning(f"No UniProt entry found for {gene_symbol} ({orf_id})")
            return {}

        canonical_accession = data["results"][0]["primaryAccession"]

        # Get isoforms
        # https://stackoverflow.com/questions/46621982/using-biopython-to-retrieve-isoform-sequences-of-a-swissprot-entry
        # key: accesion number, value: sequence
        isoforms = dict()

        # Make a call to EBI API
        r = requests.get(
            "https://www.ebi.ac.uk/proteins/api/proteins/{}/isoforms".format(
                canonical_accession
            )
        )

        # Parse the returned XML
        uniprot = ET.fromstring(r.text)

        for isoform in uniprot:
            # Get the sequence
            seq = isoform.find("{https://uniprot.org/uniprot}sequence")

            # Get the accession number
            iso_accession = isoform.find("{https://uniprot.org/uniprot}accession")

            # Store the sequence and accession number in the dictionary
            if (
                seq is not None
                and iso_accession is not None
                and seq.text
                and iso_accession.text
            ):
                isoforms[iso_accession.text] = seq.text

        # Sort dictionary by accession number
        sorted_isoforms = collections.OrderedDict(sorted(isoforms.items()))

        # Limit to the first `max_isoforms` entries
        sorted_isoforms = dict(list(sorted_isoforms.items())[:max_isoforms])

        return sorted_isoforms

    except Exception as e:
        logging.error(f"Error fetching isoforms for {gene_symbol} ({orf_id}): {e}")
        return {}


def get_protein_cigar(seq1, seq2, orf_id):
    """
    Compares two amino acid sequences and returns an extended CIGAR-like string.

    This function performs a global alignment using the Needleman-Wunsch algorithm
    with a common protein substitution matrix (BLOSUM62) and affine gap penalties.

    Args:
        seq1 (str): The first amino acid sequence.
        seq2 (str): The second amino acid sequence.
        orf_id (str): The identifier for the ORF (Open Reading Frame) being compared.

    Returns:
        The extended CIGAR-like string representing the alignment of the two sequences.

        Returns ("N/A") if no valid alignment can be found.
    """
    if not seq1 or not seq2:
        return "NA"

    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10  # Penalty for opening a gap in seq1
    aligner.extend_gap_score = -0.5  # Penalty for extending a gap in seq1
    aligner.target_open_gap_score = -10  # Penalty for opening a gap in seq2
    aligner.target_extend_gap_score = -0.5  # Penalty for extending a gap in seq2
    aligner.mode = "global"

    # Perform the alignment
    # aligner(seq1, seq2) returns an Alignments object which is iterable.
    # It might return multiple equally-scoring alignments. We usually take the first.
    alignments = list(aligner.align(seq1, seq2))

    if not alignments:
        logging.warning(f"No valid alignments found for {orf_id}.")
        return "NA"

    # Get the best alignment (usually the first one returned)
    best_alignment = alignments[0]

    # Manually construct the CIGAR-like string
    gapped_seq1 = str(best_alignment[0])  # Get the gapped query sequence
    gapped_seq2 = str(best_alignment[1])  # Get the gapped target sequence

    extended_cigar_ops_raw = []
    for i in range(len(gapped_seq1)):
        char1 = gapped_seq1[i]
        char2 = gapped_seq2[i]

        if char1 != "-" and char2 != "-":
            if char1 == char2:
                extended_cigar_ops_raw.append("=")  # Exact match
            else:
                extended_cigar_ops_raw.append("X")  # Mismatch
        elif char1 == "-":
            extended_cigar_ops_raw.append(
                "I"
            )  # Insertion in query (seq1 has gap, seq2 has char)
        elif char2 == "-":
            extended_cigar_ops_raw.append(
                "D"
            )  # Deletion from query (seq1 has char, seq2 has gap)

    # Consolidate consecutive operations into CIGAR format (e.g., MMM -> 3M)
    cigar_string = ""
    if extended_cigar_ops_raw:
        for k, g in itertools.groupby(extended_cigar_ops_raw):
            cigar_string += str(len(list(g))) + k

    return cigar_string


def main():
    # Load GPSW gene summary CSV file
    logging.info(f"Loading gene summary CSV file: {args.csv_file}")
    if not os.path.exists(args.csv_file):
        logging.error(f"Gene summary CSV file '{args.csv_file}' not found.")
        sys.exit(1)
    gene_summary_df = pd.read_csv(args.csv_file)

    if not args.all_orfs:
        logging.info(
            "Filtering gene summary DataFrame for stabilised/destabilised ORFs."
        )
        gene_summary_df = gene_summary_df[
            (gene_summary_df["stabilised"] == True)
            | (gene_summary_df["destabilised"] == True)
        ]
        orfs_to_keep = gene_summary_df["orf_id"].unique().tolist()

    logging.info(f"Loading annotation file: {args.file}")
    delimiter = detect_csv_delimiter(args.file)
    if not delimiter:
        logging.error("Could not detect CSV delimiter. Exiting.")
        sys.exit(1)
    annotation = pd.read_csv(args.file, delimiter=delimiter)

    if not args.all_orfs:
        logging.info(
            "Filtering annotation DataFrame for all ranked ORFs in GPSW gene summary."
        )
        annotation = annotation[annotation["orf_id"].isin(orfs_to_keep)]

    if args.gene_column not in annotation.columns:
        logging.error(f"Gene column '{args.gene_column}' not found in annotation file.")
        sys.exit(1)
    if args.orf_column not in annotation.columns:
        logging.error(f"ORF column '{args.orf_column}' not found in annotation file.")
        sys.exit(1)
    if args.aminoacid_column not in annotation.columns:
        logging.error(
            f"Amino acid column '{args.aminoacid_column}' not found in annotation file."
        )
        sys.exit(1)

    gene_symbols = annotation[args.gene_column].tolist()
    orf_ids = annotation[args.orf_column].tolist()
    orf_sequences_raw = annotation[args.aminoacid_column].tolist()

    # Clean amino acid sequences
    orf_amino_acid_sequences = [
        re.sub(r"[^A-Za-z]", "", seq.upper()) if isinstance(seq, str) else ""
        for seq in orf_sequences_raw
    ]

    results = []

    logging.info(
        f"Retrieving UniProt isoform sequences and comparing for {len(gene_symbols)} entries."
    )

    for gene_symbol, orf_id, orf_seq in tqdm(
        zip(gene_symbols, orf_ids, orf_amino_acid_sequences), total=len(gene_symbols)
    ):
        isoform_dict = get_all_uniprot_isoforms(
            gene_symbol, orf_id, args.organism, args.max_isoforms
        )

        result_entry = {
            "orf_id": orf_id,
            "orf_amino_acid_sequence": orf_seq,
        }

        for idx, (acc, iso_seq) in enumerate(isoform_dict.items(), start=1):
            cigarx = get_protein_cigar(iso_seq, orf_seq, orf_id)
            result_entry[f"uniprot_accession_{idx}"] = acc
            result_entry[f"variant_sequence_{idx}"] = iso_seq
            result_entry[f"cigarx_string_{idx}"] = cigarx

        results.append(result_entry)

    # Combine all the results into a dataframe
    results_df = pd.DataFrame(results)

    # Merge results back into the gene summary DataFrame
    logging.info("Merging results with gene summary data.")
    gene_summary_df = pd.merge(
        gene_summary_df,
        results_df,
        on="orf_id",
        how="left",
    )

    # Save to new file
    outfile = args.csv_file.replace(".csv", "_annotated.csv")
    gene_summary_df.to_csv(outfile, index=False, na_rep="NA")
    logging.info(f"Results saved to '{outfile}'.")
    logging.info("Done!")


if __name__ == "__main__":
    # Initialize logging
    initialise_logging()

    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Check ORF amino acid sequences against canonical protein sequences from UniProt."
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
        help="Column name that contains the ORF library amino acid sequences (from annotation file)",
    )

    parser.add_argument(
        "-c",
        "--csv-file",
        type=str,
        help="Gene summary CSV (GPSW output)",
    )

    parser.add_argument(
        "--all_orfs",
        type=bool,
        default=False,
        help="Annotate all ORFS or only those that are stabilised/destabilised",
    )

    parser.add_argument(
        "--max-isoforms",
        type=int,
        default=5,
        help="Maximum number of UniProt isoforms to compare against (default: 5)",
    )

    parser.add_argument(
        "--organism",
        type=str,
        default="Human",
        help="Organism name (e.g., 'Human', 'Mouse')",
    )

    # Parse the command line arguments
    args = parser.parse_args()

    # Log initial command
    logging.debug(f"Command executed: {' '.join(sys.argv)}")

    # Call the main function to execute the script
    main()
