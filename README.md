# gps-orfeome-tools
Extra tools for Global Protein Stability profiling analysis

## Installation of software dependencies

```bash
$ conda env create -f environment.yml
```

# Instructions

## check_sequences.py

### Description

This script checks ORF amino acid sequences against canonical protein isoform sequences from UniProt. It retrieves the isoform protein sequences for a given gene symbol, aligns them to the ORF amino acid sequence and generates CIGARx strings to visualize the alignment.

### Usage

```console
$ python check_sequences.py --help
usage: check_sequences.py [-h] [-f FILE] [-g GENE_COLUMN] [-o ORF_COLUMN] [-a AMINOACID_COLUMN] [-c CSV_FILE] [--organism ORGANISM]

Check ORF amino acid sequences against canonical protein sequences from UniProt.

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Annotation file containing gene symbols and amino acid sequences to check
  -g GENE_COLUMN, --gene-column GENE_COLUMN
                        Column name that contains gene symbols to retrieve the canonical protein sequence for (from annotation file)
  -o ORF_COLUMN, --orf-column ORF_COLUMN
                        Column name that contains ORF identifiers to retrieve the canonical protein sequence for (from annotation file)
  -a AMINOACID_COLUMN, --aminoacid-column AMINOACID_COLUMN
                        Column name that contains the ORF library amino acid sequences (from annotation file)
  -c CSV_FILE, --csv-file CSV_FILE
                        Gene summary CSV (GPSW output)
  --organism ORGANISM   Organism name (e.g., 'Human', 'Mouse')
```

### Output

The output will include the following files:
- A CSV file containing the alignment results, including CIGARx strings for each ORF.
- A log file with detailed information about the alignment process and any errors encountered.


## plot_profiles.R

### Description

This script generates barcode profiles for specified ORFs from a CSV file in a grid layout.

### Usage

```r
Rscript plot_profiles.R <barcode_summary.csv> <orf_names.txt> <dimensions> <twinpeaks>
```
- `<barcode_summary.csv>`: Path to the CSV file containing barcode summary data.
- `<orf_names.txt>`: Path to a text file containing ORF names (one per line).
- `<dimensions>`: String indicating the dimensions (rows x columns) of the plot (e.g., "2x3").
- `<twinpeaks>`: String indicating whether to include barcodes with twin peaks ("keeptwinpeaks=TRUE" or "keeptwinpeaks=FALSE").

### Output

The script will generate a PDF file named `barcode_profiles.pdf` containing the barcode profiles for the specified ORFs. The profiles will be arranged in a grid layout based on the provided dimensions.

## create_protein_fasta.py

### Description

This script creates a FASTA file containing protein sequences for each ORF in the provided CSV file. 

### Usage

```console
$ python create_protein_fasta.py --help
usage: create_protein_fasta.py [-h] -i INPUT [-g GENE_COLUMN] [-o ORF_COLUMN] -s SEQUENCE_COLUMN

Convert ORF amino acid sequences from CSV to FASTA format.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input CSV file containing ORF sequences.
  -g GENE_COLUMN, --gene-column GENE_COLUMN
                        Name of the column containing gene names.
  -o ORF_COLUMN, --orf-column ORF_COLUMN
                        Name of the column containing ORF names.
  -s SEQUENCE_COLUMN, --sequence-column SEQUENCE_COLUMN
                        Name of the column containing amino acid sequences.
```

When both GENE_COLUMN and ORF_COLUMN are provided, the script will use both columns to create headers in the FASTA file. If either column name is parsed, only that column will be used for the header.

### Output

The script will generate a FASTA file named `orf_sequences.fasta` containing the protein sequences for each ORF in the provided CSV file. Each sequence will be prefixed with a header line containing the ORF name and gene name.