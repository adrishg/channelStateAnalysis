import os
import pandas as pd
from stateAnalysis_tools import measure_shortest_distance

# Folder containing the PDB files
folder_path = '/Users/adrianahernandezgonzalez/Documents/YarovLab/statesModels/saved_output_VSD1_v4-selected/VSD1CaV12HS8WEAlocalrunv2_08453_256_512_10/pdb/'

# List of residue pairs for distance measurement
residue_pairs = [
    ('GLU', 49, 'ARG', 99),
    ('GLU', 49, 'ARG', 102),
    ('GLU', 49, 'ARG', 105),
    ('GLU', 49, 'ARG', 108),
    ('GLU', 59, 'ARG', 99),
    ('GLU', 59, 'ARG', 102),
    ('GLU', 59, 'ARG', 105),
    ('GLU', 59, 'ARG', 108)
]

# Initialize an empty list to store the data
data = []

# Loop through all PDB files in the folder
for pdb_file in os.listdir(folder_path):
    if pdb_file.endswith('.pdb'):
        pdb_path = os.path.join(folder_path, pdb_file)
        row = {'pdb_file': pdb_file}
        
        # Measure shortest distance for each residue pair
        for residue1_type, residue1_number, residue2_type, residue2_number in residue_pairs:
            distance = measure_shortest_distance(pdb_path, residue1_type, residue1_number, residue2_type, residue2_number)
            pair_label = f"{residue1_type}{residue1_number}-{residue2_type}{residue2_number}"
            row[pair_label] = distance
        
        # Add the row to the data list
        data.append(row)

# Convert the data list to a DataFrame
df = pd.DataFrame(data)

# Write the DataFrame to a CSV file
csv_file = '08-16-2024_shortest_distances_VSD1_256_512.csv'
df.to_csv(csv_file, index=False)

print(f"Results saved to {csv_file}")
