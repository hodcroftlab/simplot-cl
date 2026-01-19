# SimPlot-CL: A Command-Line Similarity Plot Generator

SimPlot-CL is a simplified, Python-based command-line reimplementation of the classic SimPlot program (Lole et al., 1999). It lets you generate similarity plots (SimPlots) and pairwise similarity tables directly from sequence alignments, without the need for a GUI.

I built this mainly for my own viral genomics work, where I wanted to run lots of SimPlot analyses automatically instead of clicking through the Windows interface a hundred times.
If you have the same problem, this might save you some time too.

## What this does

Given one or more aligned viral genome sequences, the script:
1. Splits the alignment into overlapping windows of a chosen size.
2. Calculates the pairwise similarity between a query and other sequences in each window.
3. Produces:
    a) a plot showing how similarity changes along the genome, and
    b) a CSV table with similarity values (optional).

![example simplot](simplots/OP137282.1_PV3_simplot.png)


## How similarity is calculated

For each sliding window:
- Count how many positions differ between the query and a reference sequence (Hamming distance).
- Divide that by the number of valid positions (ignoring gaps and Ns depending on your settings).
- That gives the p-distance: $p = \text{differences} / \text{valid positions}$.
- Then $\text{similarity} = 1 - p$.
That’s what’s plotted along the genome.

## Requirements
Requires Python ≥ 3.9 and the following libraries:

```
pandas
numpy
biopython
matplotlib
argcomplete
```

## Usage

The script can run in two main ways:

1️⃣ **Using one alignment file (specify a query ID)**

```python simplot.py -a alignment.fasta -q Query1 ```

2️⃣ **Using separate query and reference alignments**

```python simplot.py -a queries.fasta -r references.fasta ```

Then you can tweak window/step size, output directories, metadata, colors, etc.

## Arguments

| Flag                             | Description                                                                    |
| -------------------------------- | ------------------------------------------------------------------------------ |
| `-a`, `--alignment`              | Path to the main alignment (FASTA).                                            |
| `-q`, `--query-id`               | ID of the query sequence within that alignment (mutually exclusive with `--reference-alignment`).                                |
| `-r`, `--reference-alignment`    | Path to a separate reference alignment (must be the same nucleotide length; mutually exclusive with `--query-id`).   |
| `-n`, `--no-align`    | If set, skip MAFFT alignment before similarity plotting. Else, align sequences before plotting using `mafft --auto`.  |
| `-t`, `--threads`    | Number of threads to use for MAFFT alignment (default: 1).  |
| `-m`, `--metadata`               | Optional CSV/TSV file with sequence info (mapping accessions to genotypes). If provided, genotype information will be added to the output plots.   |
| `-mi`, `--metadata-id-col`       | Column name in metadata for sequence IDs (default: `Accession`).                           |
| `-mg`, `--metadata-genotype-col` | Column name in metadata for genotype info (default: `Genotype`).                           |
| `-mm`, `--metadata-mode`         | Whether metadata applies to `query`, `reference`, or `both` (default: `both`). |
| `-c`, `--colors`                 | Optional file mapping genotypes to colors (`tsv` or `csv`).                    |
| `-w`, `--windowsize`             | Window size (default: 100).                                                    |
| `-s`, `--stepsize`               | Step size between windows (default: 50).                                       |
| `-g`, `--gaps`                   | How to treat gaps: 0 = skip position if one or both sequences have a gap, 1 = mismatch if one has a gap, match if both have a gap, 2 = mismatch if one has a gap, skip position if both have a gap.  |
| `-p`, `--outplots`               | Directory for plot outputs (default: `simplots/`).                             |
| `-o`, `--outcsv`                 | Directory for CSV outputs (optional; if not provided, tables will not be saved).                                          |
| `-f`, `--outformat`              | Plot format: `png`, `pdf`, `svg`, or `jpg` (default: `png`).                                      |


## Output

Each run creates:
- A similarity plot (`simplots/simplot_<query>.png`)
- A similarity table (if `--outcsv` is set)

Plots show:
- Genome position on the x-axis
- Similarity (1 − p-distance) on the y-axis
- One line per reference sequence (colored by genotype if available)

## Examples

**Simple run**

Compare one query to all others in the same alignment:

```python simplot.py -a demo_data/query_alignment.fasta -q OP137282.1 -w 200 -s 50 -p plots -o tables```

**With a separate reference alignment**

```python simplot.py -a demo_data/query_alignment.fasta -r demo_data/reference_alignment.fasta -w 300 -s 100 -p plots```

**With metadata and custom colors**

```
python simplot.py \
    -a demo_data/query_alignment.fasta \
    -r demo_data/reference_alignment.fasta \
    -m demo_data/metadata.csv \
    -c demo_data/colors.tsv 
```

Example `metadata.csv`:
```
Accession, Genotype
Query1, A
Ref1, B
Ref2, C
```

Example `colors.tsv`:
```
A	#1f77b4
B	#ff7f0e
C	#2ca02c
```

## References

Original SimPlot software: https://sray.med.som.jhmi.edu/SCRoftware/SimPlot/ <br>
Lole KS, Bollinger RC, Paranjape RS, Gadkari D, Kulkarni SS, Novak NG, Ingersoll R, Sheppard HW, Ray SC. Full-length human immunodeficiency virus type 1 genomes from subtype C-infected seroconverters in India, with evidence of intersubtype recombination. J Virol. 1999 Jan;73(1):152-60.

SimPlot++ (modern GUI version): https://github.com/Stephane-S/Simplot_PlusPlus <br>
Stéphane Samson, Étienne Lord, Vladimir Makarenkov, SimPlot++: a Python application for representing sequence similarity and detecting recombination, Bioinformatics, Volume 38, Issue 11, June 2022, Pages 3118–3120, https://doi.org/10.1093/bioinformatics/btac287