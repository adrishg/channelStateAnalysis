import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Folder path where csvs are at
folder_path = '/Users/adrianahernandezgonzalez/Documents/YarovLab/repositories/stateAnalysis/data_plddts/'
# Sample input: dictionary of aliases and CSV file paths
csv_files = {
    "8_16": folder_path+"08-16-2024_average_plddt_VSD1_8_16.csv",
    "32_64": folder_path+"08-16-2024_average_plddt_VSD1_32_64.csv",
    "256_512": folder_path+"08-16-2024_average_plddt_VSD1_256_512.csv"
    # Add more as needed
}

def plot_violin_distribution(csv_files):
    # Initialize an empty list to store the data for plotting
    plot_data = []

    # Loop through each alias and CSV file
    for alias, csv_path in csv_files.items():
        # Load the CSV file into a DataFrame
        df = pd.read_csv(csv_path)
        # Create a new DataFrame containing the alias and pLDDT values
        temp_df = pd.DataFrame({
            "Alias": [alias] * len(df),
            "pLDDT": df['average_plddt']
        })
        plot_data.append(temp_df)

    # Concatenate all the DataFrames into a single DataFrame for plotting
    plot_df = pd.concat(plot_data)

    # Create a violin plot
    sns.violinplot(x="Alias", y="pLDDT", data=plot_df, palette="Set3")

    # Customize the plot
    plt.title("Distribution of pLDDT Scores")
    plt.xlabel("Model Alias")
    plt.ylabel("pLDDT Score")
    plt.xticks(rotation=45)
    
    # Show the plot
    plt.tight_layout()
    plt.show()

# Call the function to create the violin plot
plot_violin_distribution(csv_files)
