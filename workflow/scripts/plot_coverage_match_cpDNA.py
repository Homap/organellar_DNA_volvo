import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# Use glob to find all files matching the pattern
file_paths = glob.glob("results/cpDNA_genbank_info/*.main_output.tsv")

# Define gene group mapping and order
gene_groups = {
    "ATP synthase": ["atpA", "atpB", "atpE", "atpF", "atpH", "atpI"],
    "Photosystem I": ["psaA", "psaB", "psaC", "psaJ"],
    "Photosystem II": ["psb30", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH",
                       "psbI", "psbJ", "psbK", "psbL", "psbM", "psbT", "psbZ"],
    "Rubisco": ["rbcL"],
    "Cytochrome b6/f complex": ["petA", "petB", "petD", "petG", "petL"],
    "Photosystem Biogenesis factor 1": ["pbf1"],
    "Translation": ["tufA"],
    "Chlorophyll biosynthesis": ["chlB", "chlL", "chlN"],
    "Ribosomal Proteins-Large subunit": ["rpl14", "rpl16", "rpl2", "rpl20", "rpl23", "rpl36", "rpl5"],
    "Ribosomal Proteins-Small subunit": ["rps11", "rps12", "rps14", "rps18", "rps19", "rps2",
                                         "rps2a", "rps2b", "rps3", "rps4", "rps7", "rps8", "rps9"],
    "RNA polymerase": ["rpoA", "rpoB-1", "rpoB-2", "rpoBa", "rpoBb", "rpoC1", "rpoC1a", "rpoC1b", "rpoC2"],
    "Ribosomal RNAs": ["rrn16S", "rrn23S", "rrn5S"],
    "Transfer RNAs": ["trnA-UGC", "trnMf-CAU", "trnN-GUC", "trnW-CCA", "trnfM-CAU",
                      "trnS-GCU", "trnE-UUC", "trnG-UCC", "trnH-GUG", "trnL-UAG",
                      "trnR-ACG", "trnD-GUC", "trnR-UCU", "trnQ-UUG", "trnT-UGU", "trnL-UAA"],
    "Miscellaneous proteins": ["ccsA", "cemA", "clpP1", "ftsH", "pafI", "pafII", "ycf1"]
}
group_order = list(gene_groups.keys())

# Invert the mapping for easy lookup
gene_to_group = {gene: group for group, genes in gene_groups.items() for gene in genes}

# Read all files into a single DataFrame
dataframes = []
for file in file_paths:
    df = pd.read_csv(file, sep="\t")
    
    # Exclude GeneName starting with "orf" or containing "ORF" (case-insensitive)
    df = df[~df['GeneName'].str.contains('orf', case=False)]
    
    # Add a group column using the mapping
    df['Group'] = df['GeneName'].map(gene_to_group).fillna("NA")
    
    dataframes.append(df)

combined_df = pd.concat(dataframes, ignore_index=True)

# Ensure numeric columns are properly converted
combined_df['Coverage'] = pd.to_numeric(combined_df['Coverage'], errors='coerce')
combined_df['Match'] = pd.to_numeric(combined_df['Match'], errors='coerce')

# Melt the DataFrame for paired boxplots (Coverage and Match)
melted_df = combined_df.melt(id_vars=['Group', 'GeneName'], value_vars=['Coverage', 'Match'], 
                             var_name='Metric', value_name='Value')

# === Figure 1: Grouped by Gene Category (Exclude NA) ===
grouped_df = melted_df[melted_df['Group'] != "NA"]  # Exclude NA group

plt.figure(figsize=(18, 12))
sns.boxplot(data=grouped_df, x='Group', y='Value', hue='Metric', order=group_order, showfliers=False)

# Add the abline at 50
plt.axhline(50, color='red', linestyle='--', linewidth=1.5, label='Threshold (50%)')

# Customize the plot
plt.title("Paired Boxplot for Coverage and Match Across Gene Groups (Excluding NA)", fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.xlabel("Gene Group", fontsize=14)
plt.ylabel("Percentage (%)", fontsize=14)
plt.legend(loc='upper right', fontsize=12)
plt.tight_layout()

# Save the plot
plt.savefig("results/figures/cpDNA_paired_boxplot_gene_groups.png")
plt.show()

# === Figure 2: Individual Genes with Group Names (Including NA) ===
# Sort the genes by group, with "NA" at the end
group_order_with_na = group_order + ["NA"]
melted_df['Gene_with_Group'] = melted_df['GeneName'] + " (" + melted_df['Group'] + ")"
melted_df['Gene_with_Group'] = pd.Categorical(
    melted_df['Gene_with_Group'],
    categories=sorted(
        melted_df['Gene_with_Group'].unique(),
        key=lambda x: (group_order_with_na.index(x.split(' (')[1][:-1]), x)
    ),
    ordered=True
)

plt.figure(figsize=(25, 12))
sns.boxplot(data=melted_df, x='Gene_with_Group', y='Value', hue='Metric', showfliers=False)

# Add the abline at 50
plt.axhline(50, color='red', linestyle='--', linewidth=1.5, label='Threshold (50%)')

# Customize the plot
plt.title("Paired Boxplot for Coverage and Match Across Individual Genes with Group Labels", fontsize=16)
plt.xticks(rotation=90)
plt.xlabel("Gene (Group)", fontsize=14)
plt.ylabel("Percentage (%)", fontsize=14)
plt.legend(loc='upper right', fontsize=12)
plt.tight_layout()

# Save the plot
plt.savefig("results/figures/cpDNA_paired_boxplot_individual_genes_with_groups_sorted.png")
plt.show()

# === Calculate Medians and Save as TSV ===
medians = combined_df.groupby(['GeneName', 'Group']).median(numeric_only=True).reset_index()
medians.rename(columns={'Coverage': 'MedianCoverage', 'Match': 'MedianMatch'}, inplace=True)
medians = medians[['GeneName', 'Group', 'MedianCoverage', 'MedianMatch']]
medians.to_csv("results/cpDNA_genbank_info/cpDNA_gene_medians.tsv", sep="\t", index=False)

# Display the first few rows of the median table for reference
print(medians.head())
