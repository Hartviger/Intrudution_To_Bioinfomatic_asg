Overview

This R script aligns a given query protein sequence against two datasets (human.fa and mouse.fa). It uses Biostrings and DECIPHER for alignments and crayon for colorized output.

What It Does

Reads a query AA sequence and two datasets of reference sequences.
Performs sequence alignment (local by default) using the BLOSUM62 matrix and specified gap penalties.
Calculates alignment scores and percentage identities.
Identifies top 5 best and worst alignments.
Prints colored alignments to the console:
Matches in green
Mismatches in red
Gaps in yellow
Requirements

R packages: Biostrings, DECIPHER, crayon.
FASTA files: human.fa, mouse.fa.
Usage

Install required packages:

install.packages("crayon")
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")

Place align_sequences.R in the same directory as human.fa and mouse.fa.
Run:
source("align_sequences.R")
Customization

Adjust alignment_type ("local", "global", etc.).
Change gap penalties in the pairwiseAlignment() call.
