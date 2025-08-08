#!/usr/bin/env python3
import sys
import os
import argparse

def parse_fasta_for_lengths(fasta_filename):
    """
    Parses a FASTA file to get the name and length of each sequence.
    This function does not load the sequences into memory, making it
    efficient for large genomes.

    Args:
        fasta_filename (str): The path to the input FASTA file.

    Returns:
        list: A list of tuples, where each tuple contains the
              sequence name (str) and its length (int).
        int: The total length of all sequences combined.
    """
    sequences = []
    total_length = 0
    current_seq_name = None
    current_seq_len = 0

    try:
        with open(fasta_filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # If we were tracking a sequence, save it before starting a new one.
                    if current_seq_name:
                        sequences.append((current_seq_name, current_seq_len))
                        total_length += current_seq_len

                    # Start tracking the new sequence. Take the first word after '>' as the name.
                    current_seq_name = line[1:].split()[0]
                    current_seq_len = 0
                else:
                    # This is a sequence line, add its length to the current sequence's count.
                    current_seq_len += len(line)

            # After the loop, make sure to save the very last sequence in the file.
            if current_seq_name:
                sequences.append((current_seq_name, current_seq_len))
                total_length += current_seq_len

    except FileNotFoundError:
        print(f"Error: The file '{fasta_filename}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    if not sequences:
        print("Error: No sequences found in the FASTA file.", file=sys.stderr)
        sys.exit(1)

    return sequences, total_length

def main():
    """
    Main function to parse arguments and orchestrate the BED file creation.
    """
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into N approximately equal-sized BED files. "
                    "The regions in the BED files can span across multiple sequences."
    )
    parser.add_argument(
        "fasta_file",
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "n",
        type=int,
        help="The number of BED files to create."
    )
    parser.add_argument(
        "--out_dir",
        default=".",
        help="Output directory for the BED files. Defaults to the current directory."
    )
    parser.add_argument(
        "--prefix",
        help="Prefix for the output BED files. If not provided, it's derived "
             "from the input FASTA file name."
    )

    args = parser.parse_args()

    if args.n <= 0:
        print("Error: The number of files (n) must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    # --- 1. Parse FASTA to get sequence lengths ---
    print(f"Reading sequences from '{args.fasta_file}'...")
    sequences, total_length = parse_fasta_for_lengths(args.fasta_file)
    print(f"Found {len(sequences)} sequences with a total length of {total_length} bp.")

    if total_length < args.n:
        print(f"Warning: Total sequence length ({total_length} bp) is less than the number of splits ({args.n}). "
              f"Some output files may be empty.", file=sys.stderr)

    # --- 2. Calculate chunk sizes ---
    # Distribute the total length as evenly as possible among n chunks.
    base_chunk_size, remainder = divmod(total_length, args.n)
    chunk_sizes = [base_chunk_size + 1] * remainder + [base_chunk_size] * (args.n - remainder)

    # --- 3. Generate BED files ---
    # Determine the prefix if not provided by the user.
    prefix = args.prefix or os.path.splitext(os.path.basename(args.fasta_file))[0]

    # Create the output directory if it doesn't already exist.
    os.makedirs(args.out_dir, exist_ok=True)

    seq_idx = 0       # Index for the current sequence in the `sequences` list
    pos_in_seq = 0    # 0-based position in the current sequence

    print(f"Generating {args.n} BED files in '{args.out_dir}' with prefix '{prefix}'...")

    for i, chunk_size in enumerate(chunk_sizes):
        # Construct the new filename and path
        bed_filename = f"{prefix}.splitted.{i+1}.bed"
        output_path = os.path.join(args.out_dir, bed_filename)
        remaining_in_chunk = chunk_size

        with open(output_path, 'w') as bed_file:
            # Continue writing lines until this chunk is complete
            while remaining_in_chunk > 0 and seq_idx < len(sequences):
                current_seq_name, current_seq_len = sequences[seq_idx]

                # Space available from the current position to the end of the sequence
                available_in_seq = current_seq_len - pos_in_seq

                # Determine how much to write in this BED line
                to_write = min(remaining_in_chunk, available_in_seq)

                if to_write > 0:
                    bed_start = pos_in_seq
                    bed_end = bed_start + to_write
                    bed_file.write(f"{current_seq_name}\t{bed_start}\t{bed_end}\n")

                # Update our pointers
                remaining_in_chunk -= to_write
                pos_in_seq += to_write

                # If we've reached the end of the current sequence, move to the next one
                if pos_in_seq >= current_seq_len:
                    seq_idx += 1
                    pos_in_seq = 0

    print("Done.")

if __name__ == "__main__":
    main()
