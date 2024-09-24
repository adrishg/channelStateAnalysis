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

def get_atom_coordinates(lines, residue_type, residue_number, atom_name, chain=None):
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_residue_type = line[17:20].strip()
            atom_residue_number = int(line[22:26].strip())
            atom_chain = line[21].strip()  # Chain identifier
            atom_name_field = line[12:16].strip()

            # Ensure we're only working with the specified chain, if provided
            if chain and atom_chain != chain:
                continue
            
            # Check if the residue type, number, and atom name match
            if (atom_residue_type == residue_type and atom_residue_number == residue_number 
                and atom_name_field == atom_name):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    return (x, y, z)
                except ValueError:
                    print(f"Error parsing coordinates for {residue_type}{residue_number} {atom_name}")
                    return None
    print(f"Could not find {atom_name} atom for {residue_type} {residue_number}")
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

def calculate_area_from_residues_with_chain(pdb_file, residues, chain=None):
    """
    Calculates the area formed by connecting the C-alpha atoms of the provided residues, optionally by chain.
    
    Parameters:
    pdb_file (str): The path to the PDB file.
    residues (list of tuples): A list of tuples where each tuple contains the residue type and residue number.
                               Example: [('ARG', 10), ('GLU', 20), ('SER', 30)]
    chain (str): The chain identifier (optional). If not provided, it will assume all residues are in the same chain.
    
    Returns:
    tuple: A dictionary of residue pairs and distances, and the calculated area of the polygon.
    """

    lines = parse_pdb(pdb_file)
    ca_coords = []

    # Extract the CA coordinates for each residue
    for residue_type, residue_number in residues:
        coord = get_atom_coordinates(lines, residue_type, residue_number, 'CA', chain)
        if coord:
            ca_coords.append((residue_type, residue_number, coord))
        else:
            print(f"Could not find CA atom for {residue_type}{residue_number} in chain {chain if chain else 'any chain'}")
            return None

    if len(ca_coords) < 3:
        print("At least 3 residues are required to calculate an area.")
        return None

    # Dictionary to store residue pairs and their distances
    residue_pairs_distances = {}

    # Print the residue pairs and their coordinates
    print("\nResidue pairs and coordinates:")
    for i in range(len(ca_coords)):
        res1 = ca_coords[i]
        res2 = ca_coords[(i + 1) % len(ca_coords)]  # This will loop back to the first residue to close the polygon
        dist = calculate_distance(res1[2], res2[2])

        # Create pair key in the format resThreeLetterCodeResNum-resThreeLetterCodeResNum
        pair_key = f"{res1[0]}{res1[1]}-{res2[0]}{res2[1]}"
        residue_pairs_distances[pair_key] = dist

        print(f"{res1[0]}{res1[1]} (CA: {res1[2]}) to {res2[0]}{res2[1]} (CA: {res2[2]}): Distance = {dist:.2f} Å")

    # Project the 3D coordinates onto a 2D plane (ignoring z-coordinate)
    projected_coords = [(coord[0], coord[1]) for _, _, coord in ca_coords]

    # Calculate the area using the Shoelace formula
    n = len(projected_coords)
    area = 0

    for i in range(n):
        x1, y1 = projected_coords[i]
        x2, y2 = projected_coords[(i + 1) % n]  # Connects the last vertex back to the first
        area += x1 * y2 - x2 * y1

    area = abs(area) / 2.0

    print(f"\nArea formed by the residues: {area:.2f} square Å")
    
    return residue_pairs_distances, area
