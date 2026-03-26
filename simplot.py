#!/usr/bin/env python

### Similarity plot generator v1.0.5

# Import required packages
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import argparse 
import argcomplete
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
mpl.use("agg") # Use the Agg backend to save plots without displaying them

# Set font (Arial if available, otherwise DejaVu Sans)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


# Define arguments
def get_args():
    parser = argparse.ArgumentParser(description="Similarity plot generator v1.0.0")
    parser.add_argument("-s", "--sequences", required=True, help="Path to input query sequences (fasta).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-q", "--query-id", nargs="+", help="ID/Accession of query sequence(s) in the query fasta.")
    group.add_argument("-r", "--reference-sequences", help="Path to input reference sequences (fasta; if aligned, must have same nucleotide length as query sequences).")
    parser.add_argument("-i", "--include-queries-as-refs", action="store_true", help="If set, treat other --query-id sequences as references for each query (default: excluded).")
    parser.add_argument("-n", "--no-align", action="store_true", help="If set, skip MAFFT alignment before similarity plotting.")
    parser.add_argument("-m", "--metadata", default=None, help="Path to input metadata file (tsv/csv). If provided, genotype information will be added to the output plot.")
    parser.add_argument("-mi", "--metadata-id-col", default="Accession", help="Column name in metadata file that contains sequence IDs (default: Accession).")
    parser.add_argument("-mg", "--metadata-genotype-col", default="Genotype", help="Column name in metadata file that contains genotype/grouping information (default: Genotype).")
    parser.add_argument("-mm", "--metadata-mode", choices=["reference", "query", "both"], default="both", help="Which sequences the metadata applies to (default: both):\nOptions: 'reference' = metadata applies only to reference sequences\n'query' = metadata applies only to query sequences\n'both' = metadata includes both query and reference sequences.")
    parser.add_argument("-c", "--colors", default=None, help="Path to input colors file (tsv/csv). If provided, colors will be used for each genotype in the output plot.")
    parser.add_argument("-ws", "--windowsize", type=int, default=100, help="Window size for similarity plots (default: 100).")
    parser.add_argument("-ss", "--stepsize", type=int, default=50, help="Step size for similarity plots (default: 50).")
    parser.add_argument("-dm", "--distance-model", default="pdist",
                        choices=["pdist", "jc69", "k80", "hky", "tn93"],
                        help="Distance model to use for similarity calculations (default: pdist):\n"
                             " pdist  = p-distance (raw proportion of differing sites)\n"
                             " jc69   = Jukes-Cantor 1969 (equal base frequencies, one rate)\n"
                             " k80    = Kimura 1980 (one transition rate, one transversion rate)\n"
                             " hky    = Hasegawa-Kishino-Yano 1984/85 (empirical freqs, single ts rate)\n"
                             " tn93   = Tamura-Nei 1993 (empirical freqs, separate purine/pyrimidine ts rates)\n"
                             "NOTE: For hky and tn93, base frequencies are estimated from the full\n"
                             "alignment across all sequences by default (see --local-freqs).")
    parser.add_argument("-lf", "--local-freqs", action="store_true",
                        help="If set, estimate base frequencies per window rather than from the full\n"
                             "alignment (only affects hky and tn93). This captures local variation in\n"
                             "base composition along the genome, but estimates may be noisy in short\n"
                             "windows. By default, frequencies are estimated once from the full alignment.")
    parser.add_argument("-mgf", "--max-gap-frequency", type=float, default=0.1,
                        help="Maximum allowed proportion of gap/ambiguous positions in a window (0.0-1.0).\n"
                             "If the proportion of stripped positions exceeds this threshold, the window is\n"
                             "skipped for that pair and reported as missing in the output (optional; default=0.1).\n"
                             "Example: --max-gap-frequency 0.5 skips any window where >50%% of sites were gaps.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for MAFFT alignment (default: 1).")
    parser.add_argument("-f", "--outformat", default="png", help="Output file format for the plots (png/jpg/pdf/svg, default: png).")
    parser.add_argument("-p", "--outplots", default="simplots", help="Output directory for similarity plots (default: simplots).")
    parser.add_argument("-o", "--outcsv", default=None, help="Output directory for tables with similarity results for each query (optional). If not provided, tables will not be saved.")
    parser.add_argument("-oa", "--outaln", default=None, help="Output file for alignment in fasta format (optional). If not provided, the alignment will not be saved.")
    # Plot size customization (axes width and figure height) in inches
    parser.add_argument("-wd", "--width", type=float, default=14.0, help="Width of the plotting axes area in inches (default: 14.0).")
    parser.add_argument("-ht", "--height", type=float, default=5.0, help="Height of the entire figure in inches (default: 5.0).")
    
    # Register autocompletion
    argcomplete.autocomplete(parser)
    
    return parser.parse_args()

# Function to to align sequences using MAFFT
def run_mafft(input_fasta, output_fasta, threads=1):
    mafft_cmd = ["mafft", "--auto", "--thread", str(threads), input_fasta]
    with open(output_fasta, "w") as out_f:
        subprocess.run(mafft_cmd, stdout=out_f)

# Helper: normalize sequence case and convert U->T (preserve Seq type)
def normalize_records(records):
    for rec in records:
        rec.seq = Seq(str(rec.seq).upper().replace("U", "T"))
    return records

# Function to split the alignment into windows of the given window size and step size
def split_alignment(alignment, windowsize, stepsize):
    windows = {}    # Initialize a dictionary to store the windows
    sequence_length = len(alignment[0].seq)

    # Only emit windows that fit completely within the sequence (no edge truncation)
    # This means every plotted point represents exactly windowsize sites
    first_center = windowsize // 2  # first position where a full window fits against the left edge
    last_center  = sequence_length - (windowsize - windowsize // 2)  # last position where a full window fits against the right edge

    if first_center > last_center:
        # Sequence is shorter than the window; no valid windows exist
        return windows

    for center in range(first_center, last_center + 1, stepsize):
        start = center - windowsize // 2
        end   = start + windowsize

        start = max(0, start)
        end   = min(sequence_length, end)

        if start >= end:
            continue

        window_alignment = []
        for record in alignment:
            sub_seq    = record.seq[start:end]
            new_record = record[:]
            new_record.seq = sub_seq
            window_alignment.append(new_record)

        windows[center] = window_alignment

    return windows


# ---------------------------------------------------------------------------
# Distance model functions
#
# Each function receives:
#   seq1, seq2  : 1-D numpy arrays of single characters, already filtered to
#                 valid (unambiguous nucleotide) positions
#   base_freqs  : dict of {A, C, G, T} -> float, pre-computed from the full
#                 window alignment (used by HKY85 and TN93). When None, freqs
#                 are estimated from seq1 and seq2 alone (fallback for the
#                 simple one-query / one-reference case).
#
# A return value of np.nan signals that the formula could not be evaluated
# (e.g. a log argument is <= 0 due to saturation), in which case the caller
# falls back to p-distance and emits a warning.
# ---------------------------------------------------------------------------

def _compute_base_freqs(arrays):
    """Estimate base frequencies from one or more 1-D nucleotide arrays.

    Parameters
    ----------
    arrays : list of np.ndarray
        One or more arrays of nucleotide characters (A/C/G/T only).
        Pass all sequences in the window alignment for alignment-wide
        estimation, or just [seq1, seq2] for per-pair estimation.

    Returns
    -------
    dict  {base: frequency}  with keys A, C, G, T.
    """
    combined = np.concatenate(arrays)
    n = len(combined)
    return {base: np.sum(combined == base) / n for base in ("A", "C", "G", "T")}


def _freqs_from_records(records):
    """Estimate base frequencies from a list of SeqRecords.

    Only unambiguous nucleotide positions (A/C/G/T) are counted; gaps and
    ambiguous characters are excluded. Used to compute full-alignment
    frequencies before the windowing loop when --global-freqs is set.

    Parameters
    ----------
    records : list of SeqRecord

    Returns
    -------
    dict  {base: frequency}  with keys A, C, G, T, or None if no valid sites.
    """
    NUCLEOTIDES = {"A", "C", "G", "T"}
    arrays = []
    for rec in records:
        seq = np.array(list(rec.seq))
        valid = seq[np.isin(seq, list(NUCLEOTIDES))]
        if len(valid) > 0:
            arrays.append(valid)
    if not arrays:
        return None
    return _compute_base_freqs(arrays)


def dist_pdist(seq1, seq2, base_freqs=None):
    """p-distance: proportion of differing sites.

    p = x / n
    where x is the number of differing sites and n is the sequence length.
    (Decoding Genomes, Stadler et al. 2024, eq. 5.69)

    base_freqs is accepted but not used (present for a consistent signature).
    """
    return np.sum(seq1 != seq2) / len(seq1)


def dist_jc69(seq1, seq2, base_freqs=None):
    """Jukes-Cantor 69 distance.

    d_JC69 = -(3/4) * log(1 - (4/3)*p)
    where p = x/n is the proportion of differing sites.
    (Decoding Genomes, Stadler et al. 2024, eq. 5.68-5.69)

    Assumes equal base frequencies, so alignment-wide base_freqs are not used.
    Returns np.nan if the log argument is non-positive (p >= 0.75).
    """
    p = dist_pdist(seq1, seq2)
    arg = 1.0 - (4.0 / 3.0) * p
    if arg <= 0:
        return np.nan
    return -0.75 * np.log(arg)


def dist_k80(seq1, seq2, base_freqs=None):
    """Kimura 1980 (K80 / K2P) distance.

    d_K80 = -(1/2) * log(1 - 2S - V) - (1/4) * log(1 - 2V)
    where:
      S = proportion of sites with transitional differences (A<->G or C<->T)
      V = proportion of sites with transversional differences
    (Decoding Genomes, Stadler et al. 2024, eq. 5.70)

    Assumes equal base frequencies, so alignment-wide base_freqs are not used.
    Returns np.nan if either log argument is non-positive.
    """
    n = len(seq1)
    ts = (
        ((seq1 == "A") & (seq2 == "G")) | ((seq1 == "G") & (seq2 == "A")) |
        ((seq1 == "C") & (seq2 == "T")) | ((seq1 == "T") & (seq2 == "C"))
    )
    diff = seq1 != seq2
    S = np.sum(ts) / n
    V = np.sum(diff & ~ts) / n

    arg1 = 1.0 - 2.0 * S - V
    arg2 = 1.0 - 2.0 * V
    if arg1 <= 0 or arg2 <= 0:
        return np.nan
    return -0.5 * np.log(arg1) - 0.25 * np.log(arg2)


def dist_hky(seq1, seq2, base_freqs=None):
    """Hasegawa-Kishino-Yano distance.

    Uses empirical base frequencies and a single pooled transition proportion S.
    HKY differs from TN93 in that it uses one combined S (not S1/S2).

    d_HKY = 2*(pi_T*pi_C/(pi_T+pi_C) + pi_A*pi_G/(pi_A+pi_G)) * a
            - 2*(pi_T*pi_C*(pi_A+pi_G)/(pi_T+pi_C)
                 + pi_A*pi_G*(pi_T+pi_C)/(pi_A+pi_G)
                 - (pi_T+pi_C)*(pi_A+pi_G)) * b

    a = -log(1 - S / (2*(pi_T*pi_C/(pi_T+pi_C) + pi_A*pi_G/(pi_A+pi_G)))
               - ((pi_T*pi_C*(pi_A+pi_G)/(pi_T+pi_C)
                   + pi_A*pi_G*(pi_T+pi_C)/(pi_A+pi_G)) * V)
                 / (2*(pi_T*pi_C*(pi_A+pi_G) + pi_A*pi_G*(pi_T+pi_C))))
    b = -log(1 - V / (2*(pi_T+pi_C)*(pi_A+pi_G)))

    where S = total proportion of transitional differences,
          V = proportion of transversional differences.
    (Decoding Genomes, Stadler et al. 2024, eq. 5.71-5.73)

    base_freqs : if provided, these alignment-wide frequencies are used
                 instead of estimating from seq1/seq2 alone.
    Returns np.nan if any log argument is non-positive, or if purine or
    pyrimidine frequencies are zero.
    """
    f = base_freqs if base_freqs is not None else _compute_base_freqs([seq1, seq2])
    pA, pC, pG, pT = f["A"], f["C"], f["G"], f["T"]
    pR = pA + pG
    pY = pC + pT

    if pR == 0 or pY == 0:
        return np.nan

    n = len(seq1)
    ts = (
        ((seq1 == "A") & (seq2 == "G")) | ((seq1 == "G") & (seq2 == "A")) |
        ((seq1 == "C") & (seq2 == "T")) | ((seq1 == "T") & (seq2 == "C"))
    )
    diff = seq1 != seq2
    S = np.sum(ts) / n
    V = np.sum(diff & ~ts) / n

    tc_over_y = pT * pC / pY
    ag_over_r = pA * pG / pR
    coeff_a   = 2.0 * (tc_over_y + ag_over_r)
    v_num     = tc_over_y * pR + ag_over_r * pY
    v_den     = 2.0 * (pT * pC * pR + pA * pG * pY)
    if v_den == 0:
        return np.nan
    coeff_b = 2.0 * (pT * pC * pR / pY + pA * pG * pY / pR - pY * pR)

    arg_a = 1.0 - S / coeff_a - (v_num * V) / v_den
    arg_b = 1.0 - V / (2.0 * pR * pY)
    if arg_a <= 0 or arg_b <= 0:
        return np.nan

    return coeff_a * (-np.log(arg_a)) - coeff_b * (-np.log(arg_b))


def dist_tn93(seq1, seq2, base_freqs=None):
    """Tamura-Nei 1993 distance.

    Extends HKY by using separate transition proportions for pyrimidines
    (S1: C<->T) and purines (S2: A<->G), each with their own rate parameter.

    d_TN93 = (2*pi_T*pi_C / (pi_T+pi_C)) * (a1 - (pi_A+pi_G)*b)
           + (2*pi_A*pi_G / (pi_A+pi_G)) * (a2 - (pi_T+pi_C)*b)
           + 2*(pi_T+pi_C)*(pi_A+pi_G)*b

    a1 = -log(1 - (pi_T+pi_C)*S1 / (2*pi_T*pi_C) - V / (2*(pi_T+pi_C)))
    a2 = -log(1 - (pi_A+pi_G)*S2 / (2*pi_A*pi_G) - V / (2*(pi_A+pi_G)))
    b  = -log(1 - V / (2*(pi_T+pi_C)*(pi_A+pi_G)))

    where:
      S1 = proportion of sites with C<->T differences  (pyrimidine transitions)
      S2 = proportion of sites with A<->G differences  (purine transitions)
      V  = proportion of sites with transversional differences
    (Decoding Genomes, Stadler et al. 2024, eq. 5.74-5.77)

    base_freqs : if provided, these alignment-wide frequencies are used
                 instead of estimating from seq1/seq2 alone.
    Returns np.nan if any log argument is non-positive, or if a frequency
    product needed in a denominator is zero.
    """
    f = base_freqs if base_freqs is not None else _compute_base_freqs([seq1, seq2])
    pA, pC, pG, pT = f["A"], f["C"], f["G"], f["T"]
    pR = pA + pG 
    pY = pC + pT

    if pR == 0 or pY == 0 or pA * pG == 0 or pC * pT == 0:
        return np.nan

    n = len(seq1)
    S1 = np.sum(
        ((seq1 == "C") & (seq2 == "T")) | ((seq1 == "T") & (seq2 == "C"))
    ) / n
    S2 = np.sum(
        ((seq1 == "A") & (seq2 == "G")) | ((seq1 == "G") & (seq2 == "A"))
    ) / n
    ts = (
        ((seq1 == "C") & (seq2 == "T")) | ((seq1 == "T") & (seq2 == "C")) |
        ((seq1 == "A") & (seq2 == "G")) | ((seq1 == "G") & (seq2 == "A"))
    )
    V = np.sum((seq1 != seq2) & ~ts) / n

    arg_a1 = 1.0 - (pY * S1) / (2.0 * pT * pC) - V / (2.0 * pY)
    arg_a2 = 1.0 - (pR * S2) / (2.0 * pA * pG) - V / (2.0 * pR)
    arg_b  = 1.0 - V / (2.0 * pY * pR)
    if arg_a1 <= 0 or arg_a2 <= 0 or arg_b <= 0:
        return np.nan

    a1 = -np.log(arg_a1)
    a2 = -np.log(arg_a2)
    b  = -np.log(arg_b)

    return (
        (2.0 * pT * pC / pY) * (a1 - pR * b) +
        (2.0 * pA * pG / pR) * (a2 - pY * b) +
        2.0 * pY * pR * b
    )



# Map model name -> distance function
DISTANCE_MODELS = {
    "pdist": dist_pdist,
    "jc69":  dist_jc69,
    "k80":   dist_k80,
    "hky":   dist_hky,
    "tn93":  dist_tn93,
}



# Function to calculate pairwise distances between the query sequence and all reference sequences in the alignment (query sequence should be the first sequence in the alignment)
def calculate_pairwise_distances(alignment, current_step, model="pdist", max_gap_frequency=None, global_base_freqs=None):

    # Full set of unambiguous nucleotides; gaps, Ns, and any other characters
    # are always stripped for all models.
    NUCLEOTIDES = {"A", "C", "G", "T"}

    # Get query sequence from alignment
    query_seq = np.array(list(alignment[0].seq))
    query_id  = alignment[0].id

    # Remove the query sequence from the alignment to get the reference sequences
    reference_sequences = [record for record in alignment if record.id != query_id]

    # ------------------------------------------------------------------
    # Base frequency estimation for HKY and TN93
    #
    # By default, global_base_freqs is pre-computed from the full alignment
    # in main() and passed in here, giving stable estimates consistent with
    # MEGA. If --local-freqs is set, global_base_freqs is None and frequencies
    # are estimated per window from all sequences in that window instead.
    # ------------------------------------------------------------------
    alignment_base_freqs = None
    if model in ("hky", "tn93"):
        if global_base_freqs is not None:
            alignment_base_freqs = global_base_freqs
        else:
            all_valid_seqs = []
            for record in alignment:
                seq = np.array(list(record.seq))
                all_valid_seqs.append(seq[np.isin(seq, list(NUCLEOTIDES))])
            if all_valid_seqs:
                alignment_base_freqs = _compute_base_freqs(all_valid_seqs)

    # Intialize results list
    results = []

    for record in reference_sequences:
        reference_seq = np.array(list(record.seq))
        window_len = len(reference_seq)

        # ------------------------------------------------------------------
        # Step 1 – strip any position where either sequence is not an
        # unambiguous nucleotide (gaps, Ns, or any other character).
        # This rule is applied consistently for all distance models.
        # ------------------------------------------------------------------
        valid_positions = (
            np.isin(query_seq, list(NUCLEOTIDES)) &
            np.isin(reference_seq, list(NUCLEOTIDES))
        )

        # ------------------------------------------------------------------
        # Step 2 – apply the gap-frequency threshold
        # ------------------------------------------------------------------
        seq_len_valid  = int(np.sum(valid_positions))
        proportion_valid = seq_len_valid / window_len
        gap_frequency  = 1.0 - proportion_valid

        # Skip the window if the gap frequency exceeds the threshold
        if gap_frequency > max_gap_frequency:
            print(f"        └── Skipping {query_id} vs {record.id} at step {current_step}: "
                  f"gap frequency {gap_frequency*100:.1f}% exceeds --max-gap-frequency "
                  f"{max_gap_frequency*100:.1f}%.")
            continue

        # ------------------------------------------------------------------
        # Step 3 – subset to valid positions and compute the distance
        # ------------------------------------------------------------------
        reference_valid = reference_seq[valid_positions]
        query_valid     = query_seq[valid_positions]

        # Look up the requested distance function
        dist_fn = DISTANCE_MODELS.get(model, dist_pdist)

        # Pass alignment-wide base frequencies to models that use them (HKY, TN93).
        # For all other models the argument is accepted but ignored.
        dist = dist_fn(query_valid, reference_valid, base_freqs=alignment_base_freqs)

        # If the distance formula is undefined for this window (e.g. log argument
        # <= 0 due to saturation), record NaN for both distance and similarity.
        # Matplotlib will leave a gap at this position in the plot, and the CSV
        # will contain NaN so the user can identify which windows were affected.
        if dist is None or np.isnan(dist) or dist < 0:
            print(f"        └── [WARN] {model.upper()} distance undefined for "
                  f"{query_id} vs {record.id} at step {current_step} "
                  f"(window will appear as a gap in the plot).")
            dist       = np.nan
            similarity = np.nan
        else:
            dist       = round(float(dist), 4)
            similarity = max(0.0, 1 - dist)

        # Append the result as a tuple to the results list
        results.append((query_id, record.id, current_step, dist, similarity, proportion_valid))

    return results


# Function to assign colors to the results dataframe based on metadata and/or colors mapping
def assign_colors(results_df, metadata=None, metadata_id_col=None, metadata_genotype_col=None, colors=None, metadata_mode="both"):
    """
    Assigns colors to results_df depending on metadata and/or colors mapping,
    with support for metadata_mode = reference | query | both.

    Handle five cases:
      1. metadata_mode includes reference, colors provided → color by genotype
      2. metadata_mode includes reference, no colors → default colors by genotype
      3. metadata_mode excludes reference → skip genotype merge entirely, default color by sequence ID
      4. colors only → warn and default color by sequence ID
      5. no metadata, no colors → default color by sequence ID
    """

    default_colors = plt.colormaps["tab20"]

    # --- CASE 1 & 2: metadata present and includes reference sequences ---
    if metadata is not None and metadata_mode in ["reference", "both"]:
        # Merge metadata to get genotypes
        md = metadata.rename(
            columns={metadata_id_col: "seq2", metadata_genotype_col: "genotype"}
        )
        results_df = results_df.merge(md[["seq2", "genotype"]], on="seq2", how="left")

        # --- CASE 1: colors + metadata ---
        if colors is not None and "genotype" in results_df.columns:
            results_df = results_df.merge(colors, on="genotype", how="left")

            # Fill in missing colors if any genotype not in color map
            if results_df["color"].isnull().any():
                # Drop NaNs and ensure all genotypes are strings before printing
                missing = (
                    results_df.loc[results_df["color"].isnull(), "genotype"]
                    .dropna()
                    .unique()
                )
                missing_str = [str(x) for x in missing]

                if len(missing_str) > 0:
                    print(f"        └── Missing colors for genotypes: {', '.join(missing_str)}. Assigning default palette.")

                for i, g in enumerate(missing):
                    default_color = default_colors(i % default_colors.N)
                    results_df.loc[
                        (results_df["genotype"] == g) & (results_df["color"].isnull()), "color",] = mpl.colors.to_hex(default_color)
                    
            # If any colors are missing due to missing genotype metadata, assign default colors by sequence
            if results_df["genotype"].isnull().any():
                # Find reference sequences with missing genotype
                missing_refs = results_df.loc[results_df["genotype"].isnull(), "seq2"].unique()
                missing_refs_str = [str(x) for x in missing_refs]
                if len(missing_refs_str) > 0:
                    print(f"        └── Missing genotypes for reference sequences: {', '.join(missing_refs_str)}. Assigning default palette by sequence.")
                for i, seq in enumerate(missing_refs):
                    default_color = default_colors(i % default_colors.N)
                    results_df.loc[
                        (results_df["seq2"] == seq) & (results_df["genotype"].isnull()), "color",] = mpl.colors.to_hex(default_color)


        # --- CASE 2: metadata only (no colors file) ---
        elif "genotype" in results_df.columns:
            print("        └── No colors file provided. Using default colors for genotypes.")
            unique_gts = results_df["genotype"].dropna().unique()
            for i, g in enumerate(unique_gts):
                default_color = default_colors(i % default_colors.N)
                results_df.loc[
                    results_df["genotype"] == g, "color"
                ] = mpl.colors.to_hex(default_color)

            # If any colors are missing due to missing genotype metadata, assign default colors by sequence
            if results_df["genotype"].isnull().any():
                # Find reference sequences with missing genotype
                missing_refs = results_df.loc[results_df["genotype"].isnull(), "seq2"].unique()
                missing_refs_str = [str(x) for x in missing_refs]
                if len(missing_refs_str) > 0:
                    print(f"        └── Missing genotypes for reference sequences: {', '.join(missing_refs_str)}. Assigning default palette by sequence.")
                for i, seq in enumerate(missing_refs):
                    default_color = default_colors(i % default_colors.N)
                    results_df.loc[
                        (results_df["seq2"] == seq) & (results_df["genotype"].isnull()), "color",] = mpl.colors.to_hex(default_color)

        else:
            print("        └── Metadata merge produced no genotype column. Using default colors by sequence.")
            reference_seqs = results_df["seq2"].unique()
            for i, seq in enumerate(reference_seqs):
                default_color = default_colors(i % default_colors.N)
                results_df.loc[
                    results_df["seq2"] == seq, "color"
                ] = mpl.colors.to_hex(default_color)

    # --- CASE 3: metadata present but mode excludes references ---
    elif metadata is not None and metadata_mode == "query":
        print(
            "        └── Metadata mode 'query': skipping genotype merge for reference sequences. Using default colors."
        )
        reference_seqs = results_df["seq2"].unique()
        for i, seq in enumerate(reference_seqs):
            default_color = default_colors(i % default_colors.N)
            results_df.loc[
                results_df["seq2"] == seq, "color"
            ] = mpl.colors.to_hex(default_color)

    # --- CASE 4: colors file only ---
    elif metadata is None and colors is not None:
        print(
            "        └── Colors file provided but no metadata file. Using default colors by sequence."
        )
        reference_seqs = results_df["seq2"].unique()
        for i, seq in enumerate(reference_seqs):
            default_color = default_colors(i % default_colors.N)
            results_df.loc[
                results_df["seq2"] == seq, "color"
            ] = mpl.colors.to_hex(default_color)

    # --- CASE 5: neither metadata nor colors ---
    else:
        print(
            "        └── No metadata or colors provided. Using default colors by sequence."
        )
        reference_seqs = results_df["seq2"].unique()
        for i, seq in enumerate(reference_seqs):
            default_color = default_colors(i % default_colors.N)
            results_df.loc[
                results_df["seq2"] == seq, "color"
            ] = mpl.colors.to_hex(default_color)

    return results_df


# Function to generate and save the SimPlots
def plot_simplot(results_df, outdir, outformat, query_genotype=None, windowsize=None, stepsize=None, axes_width_in=14.0, fig_height_in=5.0, model="pdist"):

    # Base margins and paddings (in inches)
    base_left_margin_in = 0.6
    right_padding_in = 0.4
    top_margin_in = 0.3
    bottom_margin_in = 0.6

    # Legend gap (extra space between axes and legend) in inches
    legend_gap_in = 0.4  # increase this to add more space between axes and legend

    # Compute axes physical height
    axes_height_in = fig_height_in - top_margin_in - bottom_margin_in
    if axes_height_in <= 0:
        axes_height_in = fig_height_in * 0.8

    # Start with a temporary extra space for legend (will be adjusted)
    temp_legend_space_in = 2.0
    initial_fig_width = base_left_margin_in + axes_width_in + temp_legend_space_in + right_padding_in

    fig = plt.figure(figsize=(initial_fig_width, fig_height_in))
    ax = fig.add_axes([
        base_left_margin_in / initial_fig_width,
        bottom_margin_in / fig_height_in,
        axes_width_in / initial_fig_width,
        axes_height_in / fig_height_in,
    ])

    # Get IDs of reference sequences (plotting one line for each)
    reference_seqs = results_df["seq2"].unique()

    # Get the query sequence ID and genotype (if available)
    query_seq = results_df["seq1"].values[0]
    if query_genotype:
        query_seq = f"{query_seq} ({query_genotype})"

    # Keep track of how many times each color has been used
    color_counts = {}
    line_styles = ['-', '--', '-.', ':']  # If the same color is used multiple times, use different line styles

    # Plot similarity for each reference sequence
    for seq in reference_seqs:
        seq_results = results_df[results_df["seq2"] == seq]
        color = seq_results["color"].values[0]

        # Check if genotype column exists; if so, append genotype to label
        if "genotype" in seq_results.columns:
            genotype = seq_results["genotype"].values[0]
            if pd.notna(genotype):
                label = f"{seq} ({genotype})"
            else:
                label = seq
        else:
            label = seq

        # How many times has this color been used so far?
        count = color_counts.get(color, 0)

        # Pick a line style based on the count
        linestyle = line_styles[count % len(line_styles)]

        # Increment the counter for this color
        color_counts[color] = count + 1

        ax.plot(seq_results["step"], seq_results["similarity"], label=label, color=color, linestyle=linestyle)

    ax.set_title(f"Query Sequence: {query_seq}", fontsize=20)
    ax.set_xlabel("Position", fontsize=20)
    ax.set_ylabel("Similarity", fontsize=20)

    ncol = 2 if len(reference_seqs) > 14 else 1

    # Create a temporary legend to measure its size
    legend = ax.legend(
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=12,
        ncol=ncol,
        borderaxespad=0,
        frameon=False,
    )

    # Draw canvas to compute sizes in pixels
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Measure legend width in inches
    legend_bb = legend.get_window_extent(renderer=renderer)
    legend_width_in = legend_bb.width / fig.dpi

    # Measure y-axis label and tick label widths to ensure left margin is sufficient
    try:
        ylabel_bb = ax.yaxis.get_label().get_window_extent(renderer=renderer)
        ylabel_width_in = ylabel_bb.width / fig.dpi
    except Exception:
        ylabel_width_in = 0.0

    # Some tick labels may be empty; consider only visible ones
    ytick_bboxes = [t.get_window_extent(renderer=renderer) for t in ax.get_yticklabels() if t.get_text() != ""]
    max_ytick_width_in = max((bb.width for bb in ytick_bboxes), default=0.0) / fig.dpi

    # Compute required left margin: ensure there's enough space for ticks + ylabel + a small padding
    required_label_space_in = max(ylabel_width_in, max_ytick_width_in)
    desired_left_margin_in = max(base_left_margin_in, required_label_space_in + 0.35)  # 0.35in padding

    # Compute the new figure width so that axes keep their physical size and legend has enough space
    new_fig_width = desired_left_margin_in + axes_width_in + legend_gap_in + legend_width_in + right_padding_in

    # Resize the figure and recompute axes position so the axes keep their physical size
    fig.set_size_inches(new_fig_width, fig_height_in)

    # Recompute axes position in figure-relative coordinates and set it
    new_left = desired_left_margin_in / new_fig_width
    new_axes_width_frac = axes_width_in / new_fig_width
    ax.set_position([new_left, bottom_margin_in / fig_height_in, new_axes_width_frac, axes_height_in / fig_height_in])

    # Remove and recreate legend so it is placed correctly after resizing and with extra gap
    legend.remove()
    legend = ax.legend(
        loc="center left",
        bbox_to_anchor=(1.0 + (legend_gap_in / new_fig_width), 0.5),
        fontsize=12,
        ncol=ncol,
        borderaxespad=0,
        frameon=False,
    )

    # Plot parameter choices in the bottom left corner of the plot (window size, step size, distance model)
    if windowsize and stepsize:
        ax.text(0.01, -0.15, f"Window size: {windowsize} | Step size: {stepsize} | Distance model: {model.upper()}", transform=ax.transAxes, fontsize=12, va="top", ha="left")

    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=16)
    y_min = np.nanmin(results_df["similarity"]) - 0.02
    ax.set_ylim(y_min, 1.02)

    # Clean file-safe query name
    query_seq_fname = query_seq.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
    output_fname = f"{outdir}/{query_seq_fname}_{model}_simplot.{outformat}"

    # Save using the current figure size (which accounts for legend)
    plt.savefig(output_fname, dpi=fig.dpi, bbox_inches="tight")
    print(f"        └── SimPlot saved to {output_fname}")
    ax.clear()
    plt.close(fig)


# Main function
def main():
    args = get_args()

    # Read query sequences
    query_sequences = list(SeqIO.parse(args.sequences, "fasta"))
    
    # Normalize sequence case and convert U to T
    query_sequences = normalize_records(query_sequences)

    # Read metadata if provided
    if args.metadata:
        if not os.path.exists(args.metadata):
            raise ValueError(f"Metadata file {args.metadata} does not exist.")
        
        # Read metadata
        if args.metadata.endswith(".csv"):
            metadata = pd.read_csv(args.metadata)
        elif args.metadata.endswith(".tsv"):
            metadata = pd.read_csv(args.metadata, sep="\t")
        else:
            raise ValueError("Please provide a valid metadata file (csv/tsv).")
        
        # Check if metadata_id_col and metadata_genotype_col (default: "Accession" & "Genotype") are in the metadata file
        if args.metadata_id_col not in metadata.columns or args.metadata_genotype_col not in metadata.columns:
            raise ValueError("Please provide --metadata_id_col and --metadata_genotype_col arguments or ensure that the metadata file contains 'Accession' and 'Genotype' columns.")
        
        # Check if metadata covers all query sequences
        if args.metadata_mode in ["query", "both"]:
            query_ids = [record.id for record in query_sequences]
            missing_queries = [qid for qid in query_ids if qid not in metadata[args.metadata_id_col].values]
            if len(missing_queries) > 0:
                print(f"[WARN] The following query IDs are missing from the metadata file: {', '.join(missing_queries)}")
    else:
        metadata = None
        
    # Read colors if provided
    if args.colors:
        if not os.path.exists(args.colors):
            raise ValueError(f"Colors file {args.colors} does not exist.")
        if args.colors.endswith(".csv"):
            colors = pd.read_csv(args.colors, header=None)
        elif args.colors.endswith(".tsv"):
            colors = pd.read_csv(args.colors, sep="\t", header=None)
        else:
            raise ValueError("Please provide a valid colors file (csv/tsv).")
        if colors.shape[1] != 2:
            raise ValueError("Colors file must have exactly two columns: 1st for genotypes, 2nd for colors.")
        colors.columns = ["genotype", "color"]
    else:
        colors = None

    # Check if output directories exist, if not create them
    if not os.path.exists(args.outplots):
        os.makedirs(args.outplots)
    if args.outcsv and not os.path.exists(args.outcsv):
        os.makedirs(args.outcsv)

    # Check if reference sequences were provided
    if args.reference_sequences:
        
        print(f"[INFO] Using reference sequences: {args.reference_sequences}")
        reference_sequences = list(SeqIO.parse(args.reference_sequences, "fasta"))
        reference_sequences = normalize_records(reference_sequences) # Normalize sequence case and convert U to T

        # If no-align flag is set, skip MAFFT alignment (assume sequences are already aligned)
        if args.no_align:
            print(f"[INFO] Alignment skipped (--no-align specified). Assuming sequences are already aligned.")
            # Check if query and reference alignments have the same length
            if len(query_sequences[0].seq) != len(reference_sequences[0].seq):
                raise ValueError("Query and reference alignments must be of the same length.")
            
            if args.outaln:
                print("[WARN] --outaln specified but --no-align flag is set. Skipping alignment output since no alignment was performed.")

        else:
            print(f"[INFO] Aligning query and reference sequences using MAFFT ...")

            # Combine query and reference sequences into a single, temporary fasta file
            combined_fasta = "temp_combined_sequences.fasta"
            with open(combined_fasta, "w") as out_f:
                SeqIO.write(query_sequences + reference_sequences, out_f, "fasta")

            # Run MAFFT
            # Check if user provided an output file for the alignment; if not, use a temporary file that will be deleted after reading the aligned sequences
            if args.outaln:
                aligned_fasta = args.outaln
            else:
                aligned_fasta = "temp_aligned_sequences.fasta"
            run_mafft(combined_fasta, aligned_fasta, threads=args.threads)

            # Read aligned sequences
            aligned_sequences = list(SeqIO.parse(aligned_fasta, "fasta"))

            # Split aligned sequences back into query and reference sequences
            query_ids = set(record.id for record in query_sequences)
            query_sequences = [record for record in aligned_sequences if record.id in query_ids]
            reference_sequences = [record for record in aligned_sequences if record.id not in query_ids]


            print(f"[INFO] Alignment completed.")
            # Remove temporary files
            os.remove(combined_fasta)
            if not args.outaln:
                os.remove(aligned_fasta)
            else:
                print(f"[INFO] Aligned sequences saved to {args.outaln}.")

        
        # Check if metadata covers all reference sequences
        if metadata is not None:
            if args.metadata_mode in ["reference", "both"]:
                ref_ids = [record.id for record in reference_sequences]
                missing_refs = [rid for rid in ref_ids if rid not in metadata[args.metadata_id_col].values]
                if len(missing_refs) > 0:
                    print(f"[WARN] The following reference IDs are missing from the metadata file: {', '.join(missing_refs)}")

        # Estimate base frequencies from the full alignment once, before the per-query
        # loop. The pool of sequences is the same for every query (all queries +
        # all references), so there is no reason to recompute this per query.
        # Ignored for pdist, jc69, and k80, which do not use base frequencies.
        global_base_freqs = None
        if args.distance_model in ("hky", "tn93") and not args.local_freqs:
            global_base_freqs = _freqs_from_records(query_sequences + reference_sequences)
            print(f"[INFO] Full-alignment base frequencies: "
                  f"A={global_base_freqs['A']:.4f}, C={global_base_freqs['C']:.4f}, "
                  f"G={global_base_freqs['G']:.4f}, T={global_base_freqs['T']:.4f}")
        elif args.distance_model in ("hky", "tn93") and args.local_freqs:
            print(f"[INFO] --local-freqs specified: base frequencies will be estimated per window from the sequences in each window.")

        # Loop through each sequence in the query alignment
        for query_record in query_sequences:
            query_id = query_record.id
            print(f"[INFO] Processing query sequence: {query_id}")
            
            # Create an alignment with the query sequence as the first sequence, followed by all reference sequences
            final_alignment = [query_record] + reference_sequences

            # Split the combined alignment into windows
            print(f"    └── Splitting alignment into windows (window size: {args.windowsize}, step size: {args.stepsize})")
            windows = split_alignment(final_alignment, args.windowsize, args.stepsize)
            
            # Initialize a list to store the results
            final_results = []
            
            # Calculate pairwise distances for each window
            print(f"    └── Calculating pairwise distances for each window")
            for step, aln in windows.items():
                window_results = calculate_pairwise_distances(alignment=aln, current_step=step, model=args.distance_model, max_gap_frequency=args.max_gap_frequency, global_base_freqs=global_base_freqs)
                final_results.extend(window_results)

            # Convert the list of results to a dataframe
            results_df = pd.DataFrame(final_results, columns=["seq1", "seq2", "step", "distance", "similarity", "proportion_valid"])
            
            # Save the results as a CSV file if output directory is provided
            if args.outcsv:
                results_df.to_csv(f"{args.outcsv}/{query_id}_{args.distance_model}_similarity_results.csv", index=False)
            
            # Assign colors to the results dataframe
            print(f"    └── Assigning colors for plotting")
            results_df = assign_colors(results_df, metadata=metadata, metadata_id_col=args.metadata_id_col, metadata_genotype_col=args.metadata_genotype_col, colors=colors, metadata_mode=args.metadata_mode)

            # If metadata is provided, get the genotype of the query sequence
            query_genotype = None
            if args.metadata and args.metadata_mode in ["query", "both"]:
                q_match = metadata.loc[metadata[args.metadata_id_col] == query_id]
                if not q_match.empty:
                    query_genotype = q_match[args.metadata_genotype_col].values[0]
                else:
                    print(f"[WARN] Query ID {query_id} not found in metadata (mode={args.metadata_mode}). Proceeding without query genotype label.")
            elif args.metadata and args.metadata_mode == "reference":
                print(f"    └── Metadata mode 'reference': skipping query genotype lookup.")
                
            # Plot the SimPlot
            print(f"    └── Creating SimPlot")
            plot_simplot(results_df, args.outplots, args.outformat, query_genotype, args.windowsize, args.stepsize, axes_width_in=args.width, fig_height_in=args.height, model=args.distance_model)

            print(f"[INFO] Finished processing query sequence: {query_id}\n============================================================")


    elif args.query_id:

        # If no-align flag is set, skip MAFFT alignment (assume sequences are already aligned)
        if args.no_align:
            print(f"[INFO] Alignment skipped (--no-align specified). Assuming sequences are already aligned.")

            if args.outaln:
                print("[WARN] --outaln specified but --no-align flag is set. Skipping alignment output since no alignment was performed.")
        else:
            print(f"[INFO] Aligning sequences using MAFFT ...")

            # Run MAFFT
            if args.outaln:
                aligned_fasta = args.outaln
            else:
                aligned_fasta = "temp_aligned_sequences.fasta"
            run_mafft(args.sequences, aligned_fasta, threads=args.threads)

            # Read aligned query sequences
            query_sequences = list(SeqIO.parse(aligned_fasta, "fasta"))

            print(f"[INFO] Alignment completed.")
            if not args.outaln:
                os.remove(aligned_fasta) # Remove temporary file
            else:
                print(f"[INFO] Aligned sequences saved to {args.outaln}.")

        # Normalize args.query_id to a list of IDs
        query_ids = args.query_id if isinstance(args.query_id, (list, tuple)) else [args.query_id]

        # Decide whether other query IDs should be used as references
        if args.include_queries_as_refs:
            # Reference set = all sequences except the current query (keeps other query IDs as refs)
            reference_pool = lambda qs: [record for record in query_sequences if record.id != qs]
        else:
            # Reference set = all sequences that were not requested as queries (exclude all other query IDs)
            reference_pool = lambda qs: [record for record in query_sequences if record.id not in query_ids]

        # Estimate base frequencies from the full alignment once, before the per-query
        # loop. All sequences are pooled together regardless of which is the current
        # query, so this is the same for every iteration.
        # Ignored for pdist, jc69, and k80, which do not use base frequencies.
        global_base_freqs = None
        if args.distance_model in ("hky", "tn93") and not args.local_freqs:
            global_base_freqs = _freqs_from_records(query_sequences)
            print(f"[INFO] Full-alignment base frequencies: "
                  f"A={global_base_freqs['A']:.4f}, C={global_base_freqs['C']:.4f}, "
                  f"G={global_base_freqs['G']:.4f}, T={global_base_freqs['T']:.4f}")
        elif args.distance_model in ("hky", "tn93") and args.local_freqs:
            print(f"[INFO] --local-freqs specified: base frequencies will be estimated per window from the sequences in each window.")

        # Process each requested query ID separately
        for query_id in query_ids:
            print(f"[INFO] Processing query sequence: {query_id}")
            query_record = next((record for record in query_sequences if record.id == query_id), None)
            if query_record is None:
                raise ValueError(f"Query ID {query_id} not found in sequences. Please provide valid query ID(s) that are present in the input fasta.")

            # Reorder the alignment to position the query sequence as the first sequence
            final_alignment = [query_record] + reference_pool(query_id)

            # Split the combined alignment into windows
            print(f"    └── Splitting alignment into windows (window size: {args.windowsize}, step size: {args.stepsize})")
            windows = split_alignment(final_alignment, args.windowsize, args.stepsize)
            
            # Initialize a list to store the results
            final_results = []
            
            # Calculate pairwise distances for each window
            print(f"    └── Calculating pairwise distances for each window")
            for step, aln in windows.items():
                window_results = calculate_pairwise_distances(alignment=aln, current_step=step, model=args.distance_model, max_gap_frequency=args.max_gap_frequency, global_base_freqs=global_base_freqs)
                final_results.extend(window_results)

            # Convert the list of results to a dataframe
            results_df = pd.DataFrame(final_results, columns=["seq1", "seq2", "step", "distance", "similarity", "proportion_valid"])
            
            # Save the results as a CSV file if output directory is provided
            if args.outcsv:
                results_df.to_csv(f"{args.outcsv}/{query_id}_{args.distance_model}_similarity_results.csv", index=False)

            # Assign colors to the results dataframe
            print(f"    └── Assigning colors for plotting")
            results_df = assign_colors(results_df, metadata=metadata, metadata_id_col=args.metadata_id_col, metadata_genotype_col=args.metadata_genotype_col, colors=colors, metadata_mode=args.metadata_mode)
            
            # If metadata is provided, get the genotype of the query sequence
            query_genotype = None
            if args.metadata and args.metadata_mode in ["query", "both"]:
                q_match = metadata.loc[metadata[args.metadata_id_col] == query_id]
                if not q_match.empty:
                    query_genotype = q_match[args.metadata_genotype_col].values[0]
                else:
                    print(f"[WARN] Query ID {query_id} not found in metadata (mode={args.metadata_mode}). Proceeding without query genotype label.")
            elif args.metadata and args.metadata_mode == "reference":
                print(f"    └── Metadata mode 'reference': skipping query genotype lookup.")

            # Plot the SimPlot
            print(f"    └── Creating SimPlot")
            plot_simplot(results_df, args.outplots, args.outformat, query_genotype, args.windowsize, args.stepsize, axes_width_in=args.width, fig_height_in=args.height, model=args.distance_model)

            print(f"[INFO] Finished processing query sequence: {query_id}\n============================================================")
    
    else:
        raise ValueError("Please provide either a reference alignment file (--reference-alignment) or a query ID (--query-id).")


if __name__ == "__main__":
    main()