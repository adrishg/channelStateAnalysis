import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data_folder = '/Users/adrianahernandezgonzalez/Documents/YarovLab/repositories/stateAnalysis/'

# Sample input: dictionary of aliases and CSV file paths
csv_files = {
    "8_16": data_folder + "10-10-2024_shortest_distances_VSD2_8-16_r0.csv",
    # Add more as needed
}

# List of column names to compare (filtered for specific pairs: GLU47-ARG103 and GLU47-LYS106)
columns_to_compare = ['GLU47-ARG103', 'GLU47-LYS106']

def plot_split_violin_distribution(csv_files, columns_to_compare):
    # Initialize an empty list to store the data for plotting
    plot_data = []

    # Loop through each alias and CSV file
    for alias, csv_path in csv_files.items():
        # Load the CSV file into a DataFrame
        df = pd.read_csv(csv_path)
        
        # Loop through each column to compare
        for column in columns_to_compare:
            if column in df.columns:
                # Create a new DataFrame containing the alias and column values
                temp_df = pd.DataFrame({
                    "Alias": [alias] * len(df),
                    "Value": df[column],
                    "Column": [column] * len(df)
                })
                plot_data.append(temp_df)
            else:
                print(f"Warning: Column {column} not found in {csv_path}")

    # Concatenate all the DataFrames into a single DataFrame for plotting
    plot_df = pd.concat(plot_data)

    # Filter only the relevant residue pairs for the split plot
    plot_df = plot_df[plot_df['Column'].isin(columns_to_compare)]

    # Create a split violin plot to contrast GLU47 with ARG103 and LYS106
    sns.violinplot(x="Alias", y="Value", hue="Column", data=plot_df, split=True, inner="quart", palette="Set3")

    # Customize the plot
    plt.title("GLU47-ARG103 vs. GLU47-LYS106: Split Violin Plot")
    plt.xlabel("Alias")
    plt.ylabel("Value")
    plt.legend(title="Column")
    plt.xticks(rotation=45)
    
    # Show the plot
    plt.tight_layout()
    plt.show()

# Call the function to create the split violin plot
plot_split_violin_distribution(csv_files, columns_to_compare)
