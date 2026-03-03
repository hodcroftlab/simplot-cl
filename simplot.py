#!/usr/bin/env python

### Similarity plot generator v1.0.0

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
    parser.add_argument("-g", "--gaps", type=int, default=0, help="How to deal with gaps (default: 0):\n 0 = skip position if one or both sequences have a gap\n 1 = mismatch if one has a gap, match if both have a gap\n 2 = mismatch if one has a gap, skip position if both have a gap.")
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
    
    # Split the alignment into windows of the given window size, centered around the step size * i
    # At edges, the window size will be smaller
    for i in range(1, len(alignment[0].seq) // stepsize + 1):
        center = i * stepsize
        start = center - windowsize // 2
        end = center + windowsize // 2
        
        if start < 0:
            start = 0
        
        if end > sequence_length:
            end = sequence_length
            
        # Initialize alignment for each window
        window_alignment = []
            
        for record in alignment:
            sub_seq = record.seq[start:end]   # Slice the sequence
            
            # Create a new SeqRecord object with the sliced sequence
            new_record = record[:]     # Make a shallow copy of the record
            new_record.seq = sub_seq   # Assign the sliced sequence
            
            # Append the new record to the window alignment
            window_alignment.append(new_record)
            
        # Store the windowed alignment in the dictionary with the center position as the key
        windows[center] = window_alignment
        
    return windows


# Function to calculate pairwise distances between the query sequence and all reference sequences in the alignment (query sequence should be the first sequence in the alignment)
def calculate_pairwise_distances(alignment, current_step, gaps):

    # Get query sequence from alignment
    query_seq = np.array(list(alignment[0].seq))
    query_id = alignment[0].id
    
    # Remove the query sequence from the alignment to get the reference sequences
    reference_sequences = [record for record in alignment if record.id != query_id]
    
    # Intialize results list
    results = []
    
    for record in reference_sequences:
        reference_seq = np.array(list(record.seq))

        if gaps == 0:
            # Mask all positions where either of the two sequences has a gap or an N
            valid_positions = (reference_seq != "-") & (reference_seq != "N") & (reference_seq != "n") & (query_seq != "-") & (query_seq != "N") & (query_seq != "n")
        
        elif gaps == 1:
            # Mask all positions where either of the two sequences has an N
            valid_positions = (reference_seq != "N") & (reference_seq != "n") & (query_seq != "N") & (query_seq != "n")

        elif gaps == 2:
            # Mask all positions where both sequences have a gap or either sequence has an N
            valid_positions = ~((reference_seq == "-") & (query_seq == "-")) & (reference_seq != "N") & (reference_seq != "n") & (query_seq != "N") & (query_seq != "n")

        else:
            raise ValueError("Please provide a valid option for --gaps: 0, 1, or 2 (check --help for more detailed information on the options).")

        # Subset the sequences to only valid positions
        reference_valid = reference_seq[valid_positions]
        query_valid = query_seq[valid_positions]

        # Update the sequence length
        seq_len_valid = len(reference_valid)

        # Calculate the proportion of valid positions
        proportion_valid = seq_len_valid / len(reference_seq)
        if proportion_valid < 0.1:
            # If less than 10% of positions are valid, skip this comparison
            print(f"        └── Skipping comparison between {query_id} and {record.id} at step {current_step} due to insufficient valid positions ({proportion_valid*100:.2f}% valid; likely caused by gaps).")
            continue

        # Calculate the number of differing positions
        nd = np.sum(query_valid != reference_valid)
        
        # Calculate the p-distance / Hamming distance
        dist = round((nd / seq_len_valid), 4)  # Round to 4 decimal places
        similarity = 1 - dist
        
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
def plot_simplot(results_df, outdir, outformat, query_genotype=None, windowsize=None, stepsize=None, axes_width_in=14.0, fig_height_in=5.0):

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

    # Plot parameter choices in the bottom left corner of the plot (window size, step size)
    if windowsize and stepsize:
        ax.text(0.01, -0.15, f"Window size: {windowsize} | Step size: {stepsize}", transform=ax.transAxes, fontsize=12, va="top", ha="left")

    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=16)
    y_min = results_df["similarity"].min() - 0.02
    ax.set_ylim(y_min, 1.02)

    # Clean file-safe query name
    query_seq_fname = query_seq.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
    output_fname = f"{outdir}/{query_seq_fname}_simplot.{outformat}"

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
                window_results = calculate_pairwise_distances(alignment=aln, current_step=step, gaps=args.gaps)
                final_results.extend(window_results)

            # Convert the list of results to a dataframe
            results_df = pd.DataFrame(final_results, columns=["seq1", "seq2", "step", "distance", "similarity", "proportion_valid"])
            
            # Save the results as a CSV file if output directory is provided
            if args.outcsv:
                results_df.to_csv(f"{args.outcsv}/{query_id}_similarity_results.csv", index=False)
            
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
            plot_simplot(results_df, args.outplots, args.outformat, query_genotype, args.windowsize, args.stepsize, axes_width_in=args.width, fig_height_in=args.height)

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
                window_results = calculate_pairwise_distances(alignment=aln, current_step=step, gaps=args.gaps)
                final_results.extend(window_results)

            # Convert the list of results to a dataframe
            results_df = pd.DataFrame(final_results, columns=["seq1", "seq2", "step", "distance", "similarity", "proportion_valid"])
            
            # Save the results as a CSV file if output directory is provided
            if args.outcsv:
                results_df.to_csv(f"{args.outcsv}/{query_id}_similarity_results.csv", index=False)

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
            plot_simplot(results_df, args.outplots, args.outformat, query_genotype, args.windowsize, args.stepsize, axes_width_in=args.width, fig_height_in=args.height)

            print(f"[INFO] Finished processing query sequence: {query_id}\n============================================================")
    
    else:
        raise ValueError("Please provide either a reference alignment file (--reference-alignment) or a query ID (--query-id).")


if __name__ == "__main__":
    main()