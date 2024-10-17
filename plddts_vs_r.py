import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from stateAnalysis_tools import extractBFactor, calculate_average_plddt


# Function to process all PDB files in a folder and calculate pLDDT values grouped by 'r' number
def process_pdb_files(folder_path):
    plddt_data = {f'r{i}': [] for i in range(13)}  # Create dictionary for r0 to r12

    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.pdb'):
            # Identify files containing '_rNUMBER_'
            match = re.search(r'_r(\d+)_', filename)
            if match:
                r_number = int(match.group(1))
                if 0 <= r_number <= 12:  # Ensure valid r number
                    pdb_file_path = os.path.join(folder_path, filename)
                    b_factors = extractBFactor(pdb_file_path, 'A')  # Assuming chain A for simplicity
                    avg_plddt = calculate_average_plddt(b_factors)
                    if avg_plddt is not None:
                        plddt_data[f'r{r_number}'].append(avg_plddt)

    return plddt_data

# Function to plot the violin plot for pLDDT values
def plot_plddt_violin(plddt_data):
    # Prepare data for plotting
    data = []
    labels = []
    for r_number, plddts in plddt_data.items():
        data.extend(plddts)
        labels.extend([r_number] * len(plddts))
    
    # Create violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(x=labels, y=data)
    plt.title('pLDDT Distribution by r Number')
    plt.xlabel('r Number')
    plt.ylabel('pLDDT')
    plt.show()

# Example usage
folder_path = '/Users/adrianahernandezgonzalez/LabNotebook/10-24/states/partialAlphaCaV12HS8HLPlocalrun_b3702_256_512_10/pdb/'
  # Update with the actual folder path
plddt_data = process_pdb_files(folder_path)
plot_plddt_violin(plddt_data)
