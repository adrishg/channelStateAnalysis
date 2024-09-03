import os
import pandas as pd
from stateAnalysis_tools import extractBFactor, calculate_average_plddt

# Folder containing the PDB files
folder_path = '/Users/adrianahernandezgonzalez/Documents/YarovLab/statesModels/saved_output_VSD1_v4-selected/VSD1CaV12HS8WEAlocalrunv2_08453_256_512_10/pdb/'

def process_pdb_files(folder_path):
    data = []

    # Loop through all PDB files in the folder
    for pdb_file in os.listdir(folder_path):
        if pdb_file.endswith('.pdb'):
            pdb_path = os.path.join(folder_path, pdb_file)
            # Extract B-factors (pLDDT scores) from the PDB file
            bFactors = extractBFactor(pdb_path, 'A')  # Assuming chain A, adjust if necessary
            # Calculate the average pLDDT score
            avg_plddt = calculate_average_plddt(bFactors)
            # Append the filename and average pLDDT to the data list
            data.append({'filename': pdb_file, 'average_plddt': avg_plddt})

    # Convert the data to a DataFrame and save it to a CSV file
    df = pd.DataFrame(data)
    csv_file = '08-16-2024_average_plddt_VSD1_256_512.csv'
    df.to_csv(csv_file, index=False)

    print(f"Results saved to {csv_file}")

# Call the function to process the PDB files
process_pdb_files(folder_path)
