import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Data taken from supplementary file:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6027877/

def calculate_cm_per_mb(df):
    """
    Calculate centimorgans per megabase (cM/Mb) for male, female, and sex-averaged values.
    Assumes the 'pos' column is in base pairs (bp) and converts to megabases (Mb).
    :param df: Pandas DataFrame with columns ['chr', 'pos', 'male', 'female', 'Sex_avg']
    :return: DataFrame with additional columns ['male_cM_Mb', 'female_cM_Mb', 'Sex_avg_cM_Mb']
    """
    df = df.copy()
    df['pos_Mb'] = df['pos'] / 1e6  # Convert position from bp to Mb
    df['male_cM_Mb'] = df['male'].diff() / df['pos_Mb'].diff()
    df['female_cM_Mb'] = df['female'].diff() / df['pos_Mb'].diff()
    df['Sex_avg_cM_Mb'] = df['Sex_avg'].diff() / df['pos_Mb'].diff()
    # First row will have NaN due to diff(), so fill with the next valid value
    df.fillna(method='bfill', inplace=True)
    return df


genetic_map = pd.read_csv("FileS2", sep="\t")
df = calculate_cm_per_mb(genetic_map)

df = df[['chr', 'pos', 'Sex_avg_cM_Mb', 'Sex_avg']]
df.columns = ['Chromosome', 'Position(bp)', 'Rate(cM/Mb)', 'Map(cM)']
df['Chromosome'] = df['Chromosome'].astype(str)

# Function to sort chromosomes correctly (numeric first, then X, Y, MT)
def sort_chromosomes(chrom_list):
    numeric_chroms = sorted([int(chrom) for chrom in chrom_list if chrom.isdigit()])
    special_chroms = sorted([chrom for chrom in chrom_list if not chrom.isdigit()], key=lambda x: (x != "X", x != "Y", x != "MT"))
    return [str(chrom) for chrom in numeric_chroms] + special_chroms

sorted_chromosomes = sort_chromosomes(df['Chromosome'].unique())
num_chromosomes = len(sorted_chromosomes)

cols = 5 
rows = math.ceil(num_chromosomes / cols)   
fig, axes = plt.subplots(rows, cols, figsize=(15, rows * 3), constrained_layout=True)
axes = axes.flatten()

for i, chrom in enumerate(sorted_chromosomes):
    ax = axes[i]
    chrom_df = df[df['Chromosome'] == chrom]

    ax.plot(chrom_df['Position(bp)'] / 1e6, chrom_df['Rate(cM/Mb)'], color='blue', alpha=0.7, linewidth=0.8)
    ax.set_title(f"Chr {chrom}", fontsize=10, fontweight='bold')
    ax.set_xlabel("Positions in Mb", fontsize=8)
    ax.set_ylabel("cM/Mb", fontsize=8)
    ax.set_ylim(0, df['Rate(cM/Mb)'].max())  
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=8)

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.savefig("rn6_map.png", dpi=300)
new_entries = pd.DataFrame([
	{'Chromosome': 'Y', 'Position(bp)': 0, 'Rate(cM/Mb)': 0, 'Map(cM)': 0},
	{'Chromosome': 'MT', 'Position(bp)': 0, 'Rate(cM/Mb)': 0, 'Map(cM)': 0}
])
df = pd.concat([df, new_entries], ignore_index=True)

print(df.tail())
output_dir = "Littrell2018_rn6"
os.makedirs(output_dir, exist_ok=True)

for chrom in df['Chromosome'].unique():
    chrom_df = df[df['Chromosome'] == chrom].reset_index(drop=True)
    output_path = os.path.join(output_dir, f"{chrom}.txt")
    chrom_df.to_csv(output_path, sep='\t', index=False)

plt.show()
