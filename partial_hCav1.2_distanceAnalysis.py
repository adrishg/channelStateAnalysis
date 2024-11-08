import os
import pandas as pd
from stateAnalysis_tools import measure_shortest_distance, measure_ca_distances

# Folder containing the PDB files

# Subsampling set 8hlp
#folder_path = '/Users/adrianahernandezgonzalez/LabNotebook/10-24/states/partialAlphaCaV12HS8HLPlocalrun_b3702_256_512_10/pdb/'
folder_path = '/Users/adrianahernandezgonzalez/LabNotebook/10-24/states/partialAlphaCaV12HS8HLPlocalrun_b3702_8_16_10/pdb/'
#folder_path = '/Users/adrianahernandezgonzalez/LabNotebook/10-24/states/partialAlphaCaV12HS8HLPlocalrun_b3702_32_64_10/pdb/'

# List of residue pairs for distance measurement
#VSD1: S2: E52, F55?, F59, E62 S4:R126,R129,R132,R135
#VSD2: S2: E443, F456, E459 S3: E491 S4:R509,R512,R515,K518,521
#VSD3: S2: (N819, F823,) D826, T830, F833, E836 S3: D883, S4:K904,R907,R910,R913,R916,R920
#VSD4: S2: N1161, N1164, F1169, F1171, E1174, T1178, F1171, M1175 ,K1178 S3: D1196 , D1206 S4:R1259,R1266,R1269,R1272,R1275,R1279

residue_pairs_shortest = [
    ('ARG', 126, 'GLU', 52),
    ('ARG', 126, 'PHE', 59),
    ('ARG', 126, 'GLU', 62),
    ('ARG', 129, 'GLU', 52),
    ('ARG', 129, 'PHE', 59),
    ('ARG', 129, 'GLU', 62),
    ('ARG', 132, 'GLU', 52),
    ('ARG', 132, 'PHE', 59),
    ('ARG', 132, 'GLU', 62),
    ('ARG', 135, 'GLU', 52),
    ('ARG', 135, 'PHE', 59),
    ('ARG', 135, 'GLU', 62),
    ('ARG', 509, 'GLU', 443),
    ('ARG', 509, 'PHE', 456),
    ('ARG', 509, 'GLU', 459),
    ('ARG', 512, 'GLU', 443),
    ('ARG', 512, 'PHE', 456),
    ('ARG', 512, 'GLU', 459),
    ('ARG', 515, 'GLU', 443),
    ('ARG', 515, 'PHE', 456),
    ('ARG', 515, 'GLU', 459),
    ('LYS', 518, 'GLU', 443),
    ('LYS', 518, 'PHE', 456),
    ('LYS', 518, 'GLU', 459),
    ('ARG', 521, 'GLU', 443),
    ('ARG', 521, 'PHE', 456),
    ('ARG', 521, 'GLU', 459),
    ('LYS', 904, 'ASP', 826),
    ('LYS', 904, 'PHE', 829),
    ('LYS', 904, 'PHE', 833),
    ('LYS', 904, 'GLU', 836),
    ('ARG', 907, 'ASP', 826),
    ('ARG', 907, 'PHE', 829),
    ('ARG', 907, 'PHE', 833),
    ('ARG', 907, 'GLU', 836),
    ('ARG', 910, 'ASP', 826),
    ('ARG', 910, 'PHE', 829),
    ('ARG', 910, 'PHE', 833),
    ('ARG', 910, 'GLU', 836),
    ('ARG', 913, 'ASP', 826),
    ('ARG', 913, 'PHE', 829),
    ('ARG', 913, 'PHE', 833),
    ('ARG', 913, 'GLU', 836),
    ('ARG', 916, 'ASP', 826),
    ('ARG', 916, 'PHE', 829),
    ('ARG', 916, 'PHE', 833),
    ('ARG', 916, 'GLU', 836),
    ('ARG', 920, 'ASP', 826),
    ('ARG', 920, 'PHE', 833),
    ('ARG', 920, 'GLU', 836),
    ('ARG', 1259, 'ASN', 1164),
    ('ARG', 1259, 'PHE', 1167),
    ('ARG', 1259, 'GLU', 1174),
    ('ARG', 1266, 'ASN', 1164),
    ('ARG', 1266, 'PHE', 1167),
    ('ARG', 1266, 'GLU', 1174),
    ('ARG', 1269, 'ASN', 1164),
    ('ARG', 1269, 'PHE', 1167),
    ('ARG', 1269, 'GLU', 1174),
    ('ARG', 1272, 'ASN', 1164),
    ('ARG', 1272, 'PHE', 1167),
    ('ARG', 1272, 'GLU', 1174),
    ('LYS', 1275, 'ASN', 1164),
    ('LYS', 1275, 'PHE', 1167),
    ('LYS', 1275, 'GLU', 1174),
    ('ARG', 1279, 'ASN', 1164),
    ('ARG', 1279, 'PHE', 1167),
    ('ARG', 1279, 'GLU', 1174)
]

residue_pairs_ca = [
    ('ARG', 126, 'THR', 44),
    ('ARG', 126, 'GLU', 52),
    ('ARG', 126, 'PHE', 59),
    ('ARG', 126, 'GLU', 62),
    ('ARG', 126, 'ALA', 69),
    ('ARG', 129, 'THR', 44),
    ('ARG', 129, 'GLU', 52),
    ('ARG', 129, 'PHE', 59),
    ('ARG', 129, 'GLU', 62),
    ('ARG', 129, 'ALA', 69),
    ('ARG', 132, 'THR', 44),
    ('ARG', 132, 'GLU', 52),
    ('ARG', 132, 'PHE', 59),
    ('ARG', 132, 'GLU', 62),
    ('ARG', 132, 'ALA', 69),
    ('ARG', 135, 'THR', 44),
    ('ARG', 135, 'GLU', 52),
    ('ARG', 135, 'PHE', 59),
    ('ARG', 135, 'GLU', 62),
    ('ARG', 135, 'ALA', 69),
    ('ARG', 509, 'ASN', 439),
    ('ARG', 509, 'GLU', 443),
    ('ARG', 509, 'PHE', 456),
    ('ARG', 509, 'GLU', 459),
    ('ARG', 509, 'MET', 464),
    ('ARG', 512, 'ASN', 439),
    ('ARG', 512, 'GLU', 443),
    ('ARG', 512, 'PHE', 456),
    ('ARG', 512, 'GLU', 459),
    ('ARG', 515, 'GLU', 443),
    ('ARG', 515, 'PHE', 456),
    ('ARG', 515, 'GLU', 459),
    ('ARG', 515, 'MET', 464),
    ('LYS', 518, 'ASN', 439),
    ('LYS', 518, 'GLU', 443),
    ('LYS', 518, 'ASN', 439),
    ('LYS', 518, 'PHE', 456),
    ('LYS', 518, 'GLU', 459),
    ('ARG', 521, 'GLU', 443),
    ('ARG', 521, 'PHE', 456),
    ('ARG', 521, 'GLU', 459),
    ('ARG', 521, 'MET', 464),
    ('LYS', 904, 'ASN', 819),
    ('LYS', 904, 'ASP', 826),
    ('LYS', 904, 'PHE', 829),
    ('LYS', 904, 'PHE', 833),
    ('LYS', 904, 'GLU', 836),
    ('LYS', 904, 'ILE', 837),
    ('ARG', 907, 'ASN', 819),
    ('ARG', 907, 'ASP', 826),
    ('ARG', 907, 'PHE', 829),
    ('ARG', 907, 'PHE', 833),
    ('ARG', 907, 'GLU', 836),
    ('ARG', 907, 'ILE', 837),
    ('ARG', 910, 'ASN', 819),
    ('ARG', 910, 'ASP', 826),
    ('ARG', 910, 'PHE', 829),
    ('ARG', 910, 'PHE', 833),
    ('ARG', 910, 'GLU', 836),
    ('ARG', 913, 'ASP', 826),
    ('ARG', 913, 'PHE', 829),
    ('ARG', 913, 'PHE', 833),
    ('ARG', 913, 'GLU', 836),
    ('ARG', 913, 'ILE', 837),
    ('ARG', 916, 'ASN', 819),
    ('ARG', 916, 'ASP', 826),
    ('ARG', 916, 'PHE', 829),
    ('ARG', 916, 'PHE', 833),
    ('ARG', 916, 'GLU', 836),
    ('ARG', 920, 'ASP', 826),
    ('ARG', 920, 'PHE', 829),
    ('ARG', 920, 'PHE', 833),
    ('ARG', 920, 'GLU', 836),
    ('ARG', 920, 'ILE', 837),
    ('ARG', 1259, 'LYS', 1157),
    ('ARG', 1259, 'ASN', 1164),
    ('ARG', 1259, 'PHE', 1167),
    ('ARG', 1259, 'GLU', 1174),
    ('ARG', 1259, 'LEU', 1179),
    ('ARG', 1266, 'LYS', 1157),
    ('ARG', 1266, 'ASN', 1164),
    ('ARG', 1266, 'PHE', 1167),
    ('ARG', 1266, 'GLU', 1174),
    ('ARG', 1266, 'LEU', 1179),
    ('ARG', 1269, 'LYS', 1157),
    ('ARG', 1269, 'ASN', 1164),
    ('ARG', 1269, 'PHE', 1167),
    ('ARG', 1269, 'GLU', 1174),
    ('ARG', 1269, 'LEU', 1179),
    ('ARG', 1272, 'LYS', 1157),
    ('ARG', 1272, 'ASN', 1164),
    ('ARG', 1272, 'PHE', 1167),
    ('ARG', 1272, 'GLU', 1174),
    ('ARG', 1272, 'LEU', 1179),
    ('LYS', 1275, 'LYS', 1157),
    ('LYS', 1275, 'ASN', 1164),
    ('LYS', 1275, 'PHE', 1167),
    ('LYS', 1275, 'LEU', 1179),
    ('LYS', 1275, 'GLU', 1174),
    ('ARG', 1279, 'LYS', 1157),
    ('ARG', 1279, 'ASN', 1164),
    ('ARG', 1279, 'PHE', 1167),
    ('ARG', 1279, 'GLU', 1174),
    ('GLU', 252, 'GLU', 595),
    ('GLU', 252, 'GLU', 1024),
    ('GLU', 252, 'GLU', 1353),
    ('GLU', 595, 'GLU', 1024),
    ('GLU', 595, 'GLU', 1353),
    ('GLU', 1024, 'GLU', 1353),
    ('LEU', 290, 'LEU', 638),
    ('LEU', 290, 'VAL', 1071),
    ('LEU', 290, 'ILE', 1402),
    ('LEU', 638, 'VAL', 1071),
    ('LEU', 638, 'ILE', 1402),
    ('ILE', 1402, 'VAL', 1071),
    ]

# S6_I: LEU290, SER 294, S6_II: LEU 638, V642, S6_III: V1071, ILE1075, S6_IV: I1401, F1405

# Initialize an empty list to store the data
data = []

# Loop through all PDB files in the folder
for pdb_file in os.listdir(folder_path):
    if pdb_file.endswith('.pdb'):
        pdb_path = os.path.join(folder_path, pdb_file)
        row = {'pdb_file': pdb_file}
        
        # Measure shortest distance for each residue pair in residue_pairs_shortest
        for residue1_type, residue1_number, residue2_type, residue2_number in residue_pairs_shortest:
            shortest_distance = measure_shortest_distance(pdb_path, residue1_type, residue1_number, residue2_type, residue2_number)
            if shortest_distance:
                shortest_distance_value = float(shortest_distance.split(": ")[1].replace(' Å', ''))
                row[f"shortest_{residue1_type}{residue1_number}-{residue2_type}{residue2_number}"] = round(shortest_distance_value, 2)
            else:
                row[f"shortest_{residue1_type}{residue1_number}-{residue2_type}{residue2_number}"] = None
        
        # Measure CA distances for each residue pair in residue_pairs_ca
        ca_distances = measure_ca_distances(pdb_path, residue_pairs_ca)
        if ca_distances:
            for pair_label_ca, ca_distance in ca_distances.items():
                row[f"CA_{pair_label_ca}"] = round(ca_distance, 2)
        
        # Add the row to the data list
        data.append(row)

# Convert the data list to a DataFrame
df = pd.DataFrame(data)

# Write the DataFrame to a CSV file
csv_file = '10-22-2024_combined_distances_subsampling_8_16.csv'
#csv_file = '10-17-2024_combined_distances_subsampling_32_64.csv'
#csv_file = '10-21-2024_combined_distances_subsampling_256_512.csv'
df.to_csv(csv_file, index=False)

print(f"Results saved to {csv_file}")