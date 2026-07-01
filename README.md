# Sequence Analysis Tool

A Python script that parses a multi-sequence FASTA file and computes 
per-sequence statistics including length, GC content, and ORF detection.

## What it does
- Reads any `.fna` or `.fasta` file with multiple sequences
- Calculates GC% for each sequence
- Detects Open Reading Frames (ORFs) across all 3 forward reading frames
- Exports results to a CSV file

## Input
A FASTA format file (e.g. TP53 CDS sequences downloaded from NCBI)

## Output
A `.csv` file with the following columns:
- `ID` — sequence identifier
- `Length` — sequence length in base pairs
- `GC %` — GC content percentage
- `ORF count` — number of ORFs >= 100 nucleotides
- `Longest` — length of the longest ORF found (0 if none)

## How to run
```bash
python sequence_analysis.py
```

## Dependencies
- Biopython
- pandas

Install with:
```bash
pip install biopython pandas
```

## Data source
Fetched directly from NCBI Nucleotide database via Biopython Entrez API.
No manual download required.
