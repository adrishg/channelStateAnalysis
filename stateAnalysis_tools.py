import math
import re

# Function to extract B-factors (which are pLDDT scores in AlphaFold models) from a PDB file
def extractBFactor(pdbFile, chainID):
    patternATOM = re.compile(
        r'^ATOM\s+'
        r'\d+\s+'
        r'CA\s+'
        r'([A-Z]{3})\s+'
        + re.escape(chainID) +
        r'\s?(\d+)\s+'
        r'[-+]?\d*\.\d+\s+'
        r'[-+]?\d*\.\d+\s+'
        r'[-+]?\d*\.\d+\s+'
        r'\d*\.\d+\s+'
        r'(\d+\.\d+)'
    )
    bFactors = {}

    with open(pdbFile, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and " CA " in line:
                atomInfo = re.search(patternATOM, line)
                if atomInfo:
                    resName = atomInfo.group(1)
                    resNumber = int(atomInfo.group(2))
                    bFactor = float(atomInfo.group(3))
                    bFactors[resName + ' ' + str(resNumber)] = bFactor

    return bFactors

# Function to calculate the average pLDDT from the extracted B-factors
def calculate_average_plddt(bFactors):
    if not bFactors:
        return None
    return sum(bFactors.values()) / len(bFactors)


# Define the side chain atoms for each amino acid, including CA
sidechain_atoms = {
    'ALA': ['CA', 'CB'], 'ARG': ['CA', 'NH2', 'NH1', 'CZ', 'NE', 'CD', 'CG', 'CB'],
    'ASN': ['CA', 'ND2', 'OD1', 'CG', 'CB'], 'ASP': ['CA', 'OD2', 'OD1', 'CG', 'CB'],
    'CYS': ['CA', 'SG', 'CB'], 'GLN': ['CA', 'NE2', 'OE1', 'CD', 'CG', 'CB'],
    'GLU': ['CA', 'OE2', 'OE1', 'CD', 'CG', 'CB'], 'GLY': ['CA'],
    'HIS': ['CA', 'NE2', 'CE1', 'CD2', 'ND1', 'CG', 'CB'], 'ILE': ['CA', 'CD1', 'CG2', 'CG1', 'CB'],
    'LEU': ['CA', 'CD2', 'CD1', 'CG', 'CB'], 'LYS': ['CA', 'NZ', 'CE', 'CD', 'CG', 'CB'],
    'MET': ['CA', 'CE', 'SD', 'CG', 'CB'], 'PHE': ['CA', 'CZ', 'CE2', 'CE1', 'CD2', 'CD1', 'CG', 'CB'],
    'PRO': ['CA', 'CD', 'CG', 'CB'], 'SER': ['CA', 'OG', 'CB'],
    'THR': ['CA', 'OG1', 'CG2', 'CB'], 'TRP': ['CA', 'CH2', 'CZ3', 'CZ2', 'CE3', 'CE2', 'NE1', 'CD2', 'CD1', 'CG', 'CB'],
    'TYR': ['CA', 'OH', 'CZ', 'CE2', 'CE1', 'CD2', 'CD1', 'CG', 'CB'], 'VAL': ['CA', 'CG2', 'CG1', 'CB']
}

# Define the extreme side chain atoms for each amino acid
extreme_sidechain_atoms = {
    'ALA': ['CB'], 'ARG': ['NH2', 'NH1'], 'ASN': ['ND2', 'OD1'], 'ASP': ['OD2', 'OD1'],
    'CYS': ['SG'], 'GLN': ['NE2', 'OE1'], 'GLU': ['OE2', 'OE1'], 'GLY': ['CA'],
    'HIS': ['NE2', 'CE1'], 'ILE': ['CD1', 'CG2'], 'LEU': ['CD2', 'CD1'], 'LYS': ['NZ'],
    'MET': ['CE'], 'PHE': ['CZ'], 'PRO': ['CD'], 'SER': ['OG'], 'THR': ['OG1'],
    'TRP': ['CH2'], 'TYR': ['OH'], 'VAL': ['CG2', 'CG1']
}

def parse_pdb(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def get_atom_coordinates(lines, residue_type, residue_number, atom_name):
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_residue_type = line[17:20].strip()
            atom_residue_number = int(line[22:26].strip())
            atom_name_field = line[12:16].strip()
            if atom_residue_type == residue_type and atom_residue_number == residue_number and atom_name_field == atom_name:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                return (x, y, z)
    return None

def calculate_distance(coord1, coord2):
    return math.sqrt(sum([(a - b) ** 2 for a, b in zip(coord1, coord2)]))

def measure_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number):
    lines = parse_pdb(pdb_file)

    distances = []

    for atom1 in sidechain_atoms[residue1_type]:
        coord1 = get_atom_coordinates(lines, residue1_type, residue1_number, atom1)
        if coord1:
            for atom2 in sidechain_atoms[residue2_type]:
                coord2 = get_atom_coordinates(lines, residue2_type, residue2_number, atom2)
                if coord2:
                    distance = calculate_distance(coord1, coord2)
                    distances.append((distance, residue1_type, residue1_number, atom1, residue2_type, residue2_number, atom2))
                    print(f"Distance between {residue1_type}{residue1_number} ({atom1}) and {residue2_type}{residue2_number} ({atom2}): {distance:.2f} Å")

    return distances

def measure_extreme_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number):
    lines = parse_pdb(pdb_file)

    distances = []

    for atom1 in extreme_sidechain_atoms[residue1_type]:
        coord1 = get_atom_coordinates(lines, residue1_type, residue1_number, atom1)
        if coord1:
            for atom2 in extreme_sidechain_atoms[residue2_type]:
                coord2 = get_atom_coordinates(lines, residue2_type, residue2_number, atom2)
                if coord2:
                    distance = calculate_distance(coord1, coord2)
                    distances.append((distance, residue1_type, residue1_number, atom1, residue2_type, residue2_number, atom2))

    return distances

def measure_shortest_distance(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number):
    distances = measure_extreme_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)
    if distances:
        shortest_distance = min(distances, key=lambda x: x[0])
        distance, res1_type, res1_num, atom1, res2_type, res2_num, atom2 = shortest_distance
        formatted_distance = round(distance, 2)  # Round to two decimal places
        print(f"\nShortest distance between {res1_type}{res1_num} ({atom1}) and {res2_type}{res2_num} ({atom2}): {formatted_distance:.2f} Å")
        return formatted_distance  # Return the rounded distance
    else:
        print("Could not find specified atoms in the PDB file.")
        return None  # Return None if no distances were found



# Example usage:
if __name__ == "__main__":
    pdb_file = 'your_protein.pdb'
    residue1_type = 'ARG'
    residue1_number = 10
    residue2_type = 'GLU'
    residue2_number = 20

    print("All distances:")
    measure_distances(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)

    print("\nShortest distance:")
    measure_shortest_distance(pdb_file, residue1_type, residue1_number, residue2_type, residue2_number)
