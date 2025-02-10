Bacterial Ori Finder
This repository contains a Python implementation to identify the ori) in bacterial genomes by analyzing nucleotide composition patterns. The script also detects DNaA boxes, which are sequences that serve as binding sites for the DnaA protein, a key regulator of bacterial DNA replication initiation. Bacterial chromosomes are typically circular, and their replication starts at a specific location called the Ori. This program determines the Ori by analyzing nucleotide frequency skews (G-C skews), C frequency distributions, and by identifying the most frequent k-mers in the predicted Ori region.

How It Works?
It works by calculating the G-C skew across the genome to pinpoint the location where the minimum skew occurs.
The minimum skew position is the putative Ori location, where the transition from the reverse strand to the forward strand occurs.
To detect DNaAboxes, it extracts a 500 bp window around the predicted Ori region and identifies the most frequently occurring k-mers.
It finds the reverse complements of candidate sequences and allows for mismatches using Hamming distance to detect imperfect DNaA boxes.

This code uses the following functions:

1. patcount(symbol, genome)
Counts occurrences of a given nucleotide in the genome.

2. symbolarray(genome, symbol)
Creates a frequency map for a given nucleotide across the genome using a sliding window.

3. skew(genome)
Computes the cumulative G-C skew across the genome.

4. min_skew(genome)
Finds positions with minimum skew, which are likely candidates for the Ori.

5. freqmap(seq, k)
Generates a frequency map of k-mers in a sequence.

6. getfreqword(seq, k)
Finds the most frequent k-mer in a sequence.

7. revcomp(seq)
Finds the reverse complement of a DNA sequence.

8. ham_dist(p, q)
Calculates the Hamming distance, which measures mismatches between two sequences.

9. app_pat_match(text, pat, ham)
Finds approximate matches of a pattern in a sequence, allowing for mismatches up to the Hamming threshold.

10. find_ori_and_dnaboxes(genome, k, ham)
Main function that identifies Ori region using minimum skew.
Extracts a 500 bp region around the Ori.
Finds frequent k-mers as potential DNaA boxes.
Checks for reverse complements.
Searches for approximate matches allowing for mismatches.
