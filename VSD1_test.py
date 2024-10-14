import os
import pandas as pd
from stateAnalysis_tools import measure_shortest_distance

# Folder containing the PDB files
folder_path = '/Users/adrianahernandezgonzalez/LabNotebook/10-24/QE/subSampled_VSD2_r0_8-16/'

#97,100,103,106,109
# List of residue pairs for distance measurement
residue_pairs = [
    ('ASP', 34, 'ARG', 97),
    ('ASP', 34, 'ARG', 100),
    ('ASP', 34, 'ARG', 103),
    ('ASP', 34, 'LYS', 106),
    ('ASP', 34, 'ARG', 109),
    ('GLU', 47, 'ARG', 97),
    ('GLU', 47, 'ARG', 100),
    ('GLU', 47, 'ARG', 103),
    ('GLU', 47, 'LYS', 106),
    ('GLU', 47, 'ARG', 109),
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
csv_file = '10-10-2024_shortest_distances_VSD2_8-16_r0.csv'
df.to_csv(csv_file, index=False)

print(f"Results saved to {csv_file}")
