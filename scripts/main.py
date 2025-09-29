### Similarity plot generator v1.0.0

# Import required packages
import pandas as pd
import numpy as np
from Bio import SeqIO
#from Bio.Seq import Seq
import argparse 
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("agg") # Use the Agg backend to save plots without displaying them

# Set font to Arial
plt.rcParams["font.family"] = "Arial"


# Define arguments
def get_args():
    parser = argparse.ArgumentParser(description="Similarity plot generator v1.0.0")
    parser.add_argument("-q", "--alignment", required=True, help="Path to input query alignment (fasta).")
    parser.add_argument("-r", "--reference", required=True, help="Path to input reference alignment (fasta, must be same length as query alignment) -OR- ID/Accession of query sequence in the query alignment.")
    parser.add_argument("-m", "--metadata", default=None, help="Path to input metadata file (tsv/csv). If provided, genotype information will be added to the output plot.")
    parser.add_argument("-mi", "--metadata_id_col", default="Accession", help="Column name in metadata file that contains sequence IDs (default: Accession).")
    parser.add_argument("-mg", "--metadata_genotype_col", default="Genotype", help="Column name in metadata file that contains genotype information (default: Genotype).")
    parser.add_argument("-c", "--colors", default=None, help="Path to input colors file (tsv/csv). If provided, colors will be used for each genotype in the output plot.")
    parser.add_argument("-w", "--windowsize", type=int, default=100, help="Window size for similarity plots (default: 100).")
    parser.add_argument("-s", "--stepsize", type=int, default=50, help="Step size for similarity plots (default: 50).")
    parser.add_argument("-g", "--gaps", type=int, default=0, help="How to deal with gaps (default: 0):\n 0 = skip position if one or both sequences have a gap\n 1 = mismatch if one has a gap, match if both have a gap\n 2 = mismatch if one has a gap, skip position if both have a gap.")
    parser.add_argument("-f", "--outformat", default="png", help="Output file format for the plots (png/jpg/pdf/svg, default: png).")
    parser.add_argument("-p", "--outplots", default="simplots", help="Output directory for similarity plots (default: simplots).")
    parser.add_argument("-o", "--outcsv", default=None, help="Output directory for tables with similarity results for each query; if not provided, tables will not be saved.")
    return parser.parse_args()


# Main function
def main():
    args = get_args()

    # Read query alignment
    query_alignment = list(SeqIO.parse(args.alignment, "fasta"))

    # Check if reference is a fasta file
    if args.reference.endswith(".fasta") or args.reference.endswith(".fa") or args.reference.endswith(".fas"):
        reference_alignment = list(SeqIO.parse(args.reference, "fasta"))

        # Check if query and reference alignments have the same length
        if len(query_alignment[0].seq) != len(reference_alignment[0].seq):
            raise ValueError("Query and reference alignments must be of the same length.")
        
        # Loop through each sequence in the query alignment
        for query_record in query_alignment:
            query_id = query_record.id
            print(f"[INFO] Processing query sequence: {query_id}")
            
            # Create an alignment with the query sequence as the first sequence, followed by all reference sequences
            final_alignment = [query_record] + reference_alignment
            
            # Split the combined alignment into windows
            print(f"    └── Splitting alignment into windows (window size: {args.windowsize}, step size: {args.stepsize})")
            windows = split_alignment(final_alignment, args.windowsize, args.stepsize)
            
            # Initialize a list to store the results
            final_results = []
            
            # Calculate pairwise distances for each window
            print(f"    └── Calculating pairwise distances for each window")
            for step, aln in windows.items():
                #window_results = calculate_pairwise_distances(alignment=aln, current_step=step, gaps=args.gaps)
                window_results = calculate_pairwise_distances(alignment=aln, current_step=step, gaps=0)
                final_results.extend(window_results)

            # Convert the list of results to a dataframe
            results_df = pd.DataFrame(final_results, columns=["seq1", "seq2", "step", "distance", "similarity", "proportion_valid"])
            
            # Check if output directories exist, if not create them
            if not os.path.exists(args.outplots):
                os.makedirs(args.outplots)
            if args.outcsv and not os.path.exists(args.outcsv):
                os.makedirs(args.outcsv)
            
            # Save the results as a CSV file if output directory is provided
            if args.outcsv:
                results_df.to_csv(f"{args.outcsv}/{query_id}_similarity_results.csv", index=False)
            
            # Plot the SimPlot
            print(f"    └── Creating SimPlot for {query_id}")
            plot_simplot(results_df, args.outplots, args.outformat, args.metadata, args.colors, args.metadata_id_col, args.metadata_genotype_col)

            print(f"[INFO] Finished processing query sequence: {query_id}\n============================================================")

    else:
        # If reference is not a file, assume it is the ID/Accession of the query sequence
        query_id = args.reference
        print(f"[INFO] Processing query sequence: {query_id}")
        query_record = [record for record in query_alignment if record.id == query_id]
        if len(query_record) == 0:
            raise ValueError(f"Reference ID {query_id} not found in query alignment. Please provide a valid reference alignment file (.fasta/.fa/.fas) or a valid query ID.")
        
        # Reorder the alignment to position the query sequence as the first sequence
        final_alignment = query_record + [record for record in query_alignment if record.id != query_id]

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
        
        # Check if output directories exist, if not create them
        if not os.path.exists(args.outplots):
            os.makedirs(args.outplots)
        if args.outcsv and not os.path.exists(args.outcsv):
            os.makedirs(args.outcsv)
        
        # Save the results as a CSV file if output directory is provided
        if args.outcsv:
            results_df.to_csv(f"{args.outcsv}/{query_id}_similarity_results.csv", index=False)
        
        # Plot the SimPlot
        print(f"    └── Creating SimPlot for {query_id}")
        plot_simplot(results_df, args.outplots, args.outformat, args.metadata, args.colors, args.metadata_id_col, args.metadata_genotype_col)

        print(f"[INFO] Finished processing query sequence: {query_id}\n============================================================")


# Function to split the alignment into windows of the given window size and step size
def split_alignment(alignment, windowsize, stepsize):
    windows = {} # Initialize a dictionary to store the windows
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
            new_record = record[:]    # Make a shallow copy of the record
            new_record.seq = sub_seq  # Assign the sliced sequence
            
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
    
    # Remove the query sequence from the alignment
    reference_alignment = [record for record in alignment if record.id != query_id]
    
    # Intialize results list
    results = []
    
    for record in reference_alignment:
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
            print(f"        └── Skipping comparison between {query_id} and {record.id} at step {current_step} due to insufficient valid positions ({proportion_valid*100:.2f}% valid).")
            continue

        # Calculate the number of differing positions
        nd = np.sum(query_valid != reference_valid)
        
        # Calculate the p-distance / Hamming distance
        dist = round((nd / seq_len_valid), 4)  # Round to 4 decimal places
        similarity = 1 - dist
        
        # Append the result as a tuple to the results list
        results.append((query_id, record.id, current_step, dist, similarity, proportion_valid))
    
    return results


# Function to generate and save the SimPlots
def plot_simplot(results_df, outdir, outformat, metadata=None, colors=None, metadata_id_col="Accession", metadata_genotype_col="Genotype"):

    default_colors = plt.colormaps["tab20"]

    if metadata:
        # Read metadata file
        if metadata.endswith(".csv"):
            meta_df = pd.read_csv(metadata)
        elif metadata.endswith(".tsv"):
            meta_df = pd.read_csv(metadata, sep="\t")
        else:
            raise ValueError("Please provide a valid metadata file (csv/tsv).")
        
        # Check if required columns are present in the metadata file
        if metadata_id_col not in meta_df.columns:
            raise ValueError(f"Metadata file does not contain a column named '{metadata_id_col}' (for sequence IDs).")
        if metadata_genotype_col not in meta_df.columns:
            raise ValueError(f"Metadata file does not contain a column named '{metadata_genotype_col}' (for genotype information).")
        
        # Merge metadata with results dataframe
        meta_df.rename(columns={metadata_id_col: "seq2", metadata_genotype_col: "genotype"}, inplace=True)
        results_df = results_df.merge(meta_df[["seq2", "genotype"]], on="seq2", how="left")

        # If color-genotype mapping is provided, add colors to the output plot; else use default colors for each genotype
        if colors:
            if colors.endswith(".csv"):
                colors_df = pd.read_csv(colors, header=None)
            elif colors.endswith(".tsv"):
                colors_df = pd.read_csv(colors, sep="\t", header=None)
            else:
                raise ValueError("Please provide a valid colors file (csv/tsv).")
        
            # Check if colors file has two columns
            if colors_df.shape[1] != 2:
                raise ValueError("Colors file must have exactly two columns: genotype and color (hex code or color name).")
            colors_df.columns = ["genotype", "color"]

            # Merge colors with results dataframe
            results_df = results_df.merge(colors_df, on="genotype", how="left")

            # Check for any genotypes without a color assigned
            if results_df["color"].isnull().any():
                missing_colors = results_df[results_df["color"].isnull()]["genotype"].unique()
                missing_colors_str = ', '.join(missing_colors)
                print(f"        └── The following genotypes do not have a color assigned: {missing_colors_str}. Assigning default colors.")
                # Assign default colors to missing genotypes
                unique_genotypes = results_df[results_df["color"].isnull()]["genotype"].unique()
                for i, genotype in enumerate(unique_genotypes):
                    default_color = default_colors(i % default_colors.N)
                    results_df.loc[(results_df["genotype"] == genotype) & (results_df["color"].isnull()), "color"] = mpl.colors.to_hex(default_color)

    # If neither metadata nor colors are provided, use default colors for each reference sequence
    else:
        reference_seqs = results_df["seq2"].unique()
        for i, seq in enumerate(reference_seqs):
            default_color = default_colors(i % default_colors.N)
            results_df.loc[results_df["seq2"] == seq, "color"] = mpl.colors.to_hex(default_color)


    # Make the similarity plot
    fig, ax = plt.subplots(figsize=(18, 5))
    
    # Get IDs of reference sequences (plotting one line for each)
    reference_seqs = results_df["seq2"].unique()
    
    # Get ID of query sequence
    query_seq = results_df["seq1"].unique()[0]
    if metadata:   # If metadata is available, add genotype information to the query sequence
        query_genotype = meta_df.loc[meta_df["seq2"] == query_seq, "genotype"].values[0]
        if pd.isna(query_genotype):
            query_seq = query_seq
        else:
            query_seq = f"{query_seq} ({query_genotype})"

    # Keep track of how many times each color has been used
    color_counts = {}
    line_styles = ['-', '--', '-.', ':'] # If the same color is used multiple times, use different line styles
    
    # Plot similarity for each reference sequence
    for seq in reference_seqs:
        seq_results = results_df[results_df["seq2"] == seq]
        color = seq_results["color"].values[0]

        # How many times has this color been used so far?
        count = color_counts.get(color, 0)

        # Pick a line style based on the count
        linestyle = line_styles[count % len(line_styles)]

        # Increment the counter for this color
        color_counts[color] = count + 1

        ax.plot(seq_results["step"], seq_results["similarity"], label=seq, color=color, linestyle=linestyle)
    
    ax.set_title(f"Query Sequence: {query_seq}", fontsize=20)
    ax.set_xlabel("Position", fontsize=20)
    ax.set_ylabel("Similarity", fontsize=20)
    ncol = 2 if len(reference_seqs) > 10 else 1
    ax.legend(
        loc="center",                   # Center the legend box itself
        bbox_to_anchor=(1.18, 0.5),     # Position relative to the right side of the plot
        fontsize=12,
        ncol=ncol,
        borderaxespad=0,                # Space between the legend and the axes
        frameon=False                   # Remove the frame for a cleaner look
    )
    
    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=16)
    y_min = results_df["similarity"].min() - 0.02
    ax.set_ylim(y_min, 1.02)
    query_seq = query_seq.replace(" ", "_").replace("(", "").replace(")", "")
    plt.tight_layout()
    plt.savefig(f"{outdir}/simplot_{query_seq}.{outformat}")
    ax.clear()
    plt.close(fig)
    

if __name__ == "__main__":
    main()