import os
import re
import sys
import argparse

# For subsampleAF2 to convert the template PDB to a sequence automatically
# Need to improve functionality to choose only certain chains in the PDB (if needed)

def get_sequence_from_pdb(known_pdb_path, get_seq_from_pdb, fallbackSequence):
    """
    Convert a PDB file to a colon-separated sequence string or return a fallback sequence.
    
    Parameters:
        known_pdb_path (str): Path to the PDB file.
        get_seq_from_pdb (bool): Flag to determine whether to extract the sequence from the PDB file.
        fallbackSequence (str): Fallback sequence to return if get_seq_from_pdb is False.
    
    Returns:
        str: Colon-separated sequence string (if extraction is enabled) or fallbackSequence.
    
    Exits the program if get_seq_from_pdb is False and fallbackSequence is empty.
    """
    if not get_seq_from_pdb:
        if not fallbackSequence:
            print("Error: Fallback sequence must be specified if get_seq_from_pdb is False.")
            sys.exit(1)
        else:
            return fallbackSequence
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
        sequences = []
        chain_dict = {}
        chain_list = []
        
        # Use the provided pdb file path
        pdb_file = known_pdb_path
        filename = os.path.basename(pdb_file).split('.')[0]

        with open(pdb_file, 'r') as fp:
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

        # Collect sequences in the order chains were found
        for chain in chain_list:
            sequences.append(chain_dict[chain])
        
        # Join sequences with a colon separator and return the result.
        output = ":".join(sequences)
        return output

def main():
    parser = argparse.ArgumentParser(
        description="Extract sequence from a PDB file as a colon-separated string, or use a fallback sequence."
    )
    parser.add_argument("pdb_path", help="Path to the PDB file")
    parser.add_argument("--fallback", type=str, default="", help="Fallback sequence if not extracting from PDB")
    # By default, get_seq_from_pdb is True; use --no_seq to set it to False.
    parser.add_argument("--no_seq", action="store_false", dest="get_seq_from_pdb",
                        help="Do not extract sequence from pdb file, use fallback sequence instead")
    args = parser.parse_args()

    result = get_sequence_from_pdb(args.pdb_path, args.get_seq_from_pdb, args.fallback)
    print(result)

if __name__ == "__main__":
    main()
