import argparse

# Convert sto (stockholm) to a3m file format for MSA.
# Primary logic from https://github.com/Dapid/stockholm_reformat/

def sto_to_a3m(sto_filename, a3m_filename):
    """
    Convert a Stockholm (.sto) file to an A3M (.a3m) file.
    
    The function assumes the Stockholm file is in standard format:
      - It may contain header lines (starting with '#' or empty lines)
      - Sequence blocks may be split across multiple sections
      - The first sequence encountered is used as the reference.
    
    In the conversion:
      - Only alignment columns where the reference sequence has a non-gap character
        (i.e. not '-' or '.') are retained.
      - The resulting alignment is written in FASTA (A3M) format.
    
    Parameters:
        sto_filename (str): Path to the input Stockholm file.
        a3m_filename (str): Path to the output A3M file.
    """
    sequences = {}
    
    # Read and parse the Stockholm file
    with open(sto_filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            # Expect sequence lines with at least two parts: name and sequence
            parts = line.split()
            if len(parts) < 2:
                continue
            name, seq_fragment = parts[0], parts[1]
            # Accumulate fragments in case the alignment is split into blocks
            sequences[name] = sequences.get(name, "") + seq_fragment
    
    if not sequences:
        raise ValueError("No sequences were found in the provided Stockholm file.")

    # Choose the first sequence as the reference (query)
    ref_seq = next(iter(sequences.values()))
    
    # Determine the columns to keep: where the reference does not have a gap
    keep_positions = [i for i, char in enumerate(ref_seq) if char not in "-."]
    
    # Reformat all sequences: keep only the characters in the positions where the reference is non-gap
    new_sequences = {}
    for name, seq in sequences.items():
        new_seq = ''.join(seq[i] for i in keep_positions)
        new_sequences[name] = new_seq
    
    # Write the reformatted alignment to the A3M file in FASTA format
    with open(a3m_filename, 'w') as out_file:
        for name, seq in new_sequences.items():
            out_file.write(f'>{name}\n{seq}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a Stockholm (.sto) file to an A3M (.a3m) file."
    )
    parser.add_argument("input_sto", help="Path to the input Stockholm (.sto) file")
    parser.add_argument("output_a3m", help="Path for the output A3M (.a3m) file")
    args = parser.parse_args()

    sto_to_a3m(args.input_sto, args.output_a3m)
    print("Conversion complete!")
