import os
import re
import sys
import argparse

# For subsampleAF2 to convert the template PDB to a sequence automatically
# Primary parsing code from https://github.com/kad-ecoli/pdb2fasta/blob/master/pdb2fasta.py

def get_sequence_from_pdb(known_pdb_path, get_seq_from_pdb, fallbackSequence, generate_fasta=False):
    """
    Convert a PDB file to a colon-separated sequence string and optionally a FASTA formatted string,
    or return a fallback sequence.
    
    Parameters:
        known_pdb_path (str): Path to the PDB file.
        get_seq_from_pdb (bool): Flag to determine whether to extract the sequence from the PDB file.
        fallbackSequence (str): Fallback sequence to return if get_seq_from_pdb is False.
        generate_fasta (bool): Flag to determine whether to generate FASTA formatted output.
    
    Returns:
        tuple: (colon_separated_sequence, fasta_output)
               If generate_fasta is False, fasta_output is an empty string.
    
    Exits the program if get_seq_from_pdb is False and fallbackSequence is empty.
    """
    if not get_seq_from_pdb:
        if not fallbackSequence:
            print("Error: Fallback sequence must be specified if get_seq_from_pdb is False.")
            sys.exit(1)
        else:
            colon_output = fallbackSequence
            fasta_output = ">fallback\n" + fallbackSequence + "\n" if generate_fasta else ""
            return colon_output, fasta_output
    else:
        # Mapping of three-letter codes to one-letter codes
        aa3to1 = {
           'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
           'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
           'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
           'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
           'MSE':'M',
        }

        # Pattern to match CA atoms in ATOM and HETATM records
        ca_pattern = re.compile(
            r"^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|"
            r"^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"
        )

        # Prepare to collect sequences from different chains
        chain_dict = {}
        chain_list = []

        # Get filename (without extension) for header information
        filename = os.path.basename(known_pdb_path).split('.')[0]

        with open(known_pdb_path, 'r') as fp:
            for line in fp:
                if line.startswith("ENDMDL"):
                    break
                match_list = ca_pattern.findall(line)
                if match_list:
                    # Concatenate groups to get the three-letter residue code
                    resn = match_list[0][0] + match_list[0][2]
                    # Concatenate groups to determine the chain identifier
                    chain = match_list[0][1] + match_list[0][3]
                    if chain in chain_dict:
                        chain_dict[chain] += aa3to1[resn]
                    else:
                        chain_dict[chain] = aa3to1[resn]
                        chain_list.append(chain)

        # Generate colon-separated sequence output (order preserved by chain_list)
        sequences = [chain_dict[chain] for chain in chain_list]
        colon_output = ":".join(sequences)

        # Optionally generate FASTA formatted output for each chain
        if generate_fasta:
            fasta_lines = []
            for chain in chain_list:
                fasta_lines.append(f">{filename}:{chain}")
                fasta_lines.append(chain_dict[chain])
            fasta_output = "\n".join(fasta_lines) + "\n"
        else:
            fasta_output = ""
        return colon_output, fasta_output

def main():
    parser = argparse.ArgumentParser(
        description="Extract sequence from a PDB file as a colon-separated string, "
                    "and optionally write the FASTA file."
    )
    parser.add_argument("pdb_path", help="Path to the PDB file")
    parser.add_argument("--fallback", type=str, default="", help="Fallback sequence if not extracting from PDB")
    # By default, get_seq_from_pdb is True; use --no_seq to set it to False.
    parser.add_argument("--no_seq", action="store_false", dest="get_seq_from_pdb",
                        help="Do not extract sequence from pdb file; use fallback sequence instead")
    # New optional argument to control FASTA generation. Default is False.
    parser.add_argument("--generate_fasta", action="store_true", default=False,
                        help="Generate FASTA formatted output (default is false)")
    parser.add_argument("--fasta_out", type=str, default=None,
                        help="Path to output the FASTA formatted sequence(s) (if not provided, FASTA is not written to file)")
    args = parser.parse_args()

    colon_output, fasta_output = get_sequence_from_pdb(
        args.pdb_path, args.get_seq_from_pdb, args.fallback, args.generate_fasta
    )
    
    # Print the colon-separated sequence string
    print(colon_output)

    # If FASTA generation is enabled and an output file is specified, write the FASTA formatted output
    if args.generate_fasta and args.fasta_out:
        try:
            with open(args.fasta_out, 'w') as f:
                f.write(fasta_output)
            print(f"FASTA file written to {args.fasta_out}")
        except IOError as e:
            print(f"Error writing FASTA file: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()
