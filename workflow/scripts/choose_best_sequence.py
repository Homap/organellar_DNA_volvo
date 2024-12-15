import os
from Bio import SeqIO
import pandas as pd
import sys

# Paths
fasta_dir = sys.argv[1]  # Update with your directory containing fasta files
table_path = sys.argv[2]  # Update with the path to your table
output_dir = sys.argv[3]  # Directory where filtered fasta files will be written
os.makedirs(output_dir, exist_ok=True)

# Load the table into a pandas DataFrame
table = pd.read_csv(table_path, sep="\t")

# Create a dictionary for quick lookup
locus_gene_map = table.set_index("GeneName")["LOCUS"].to_dict()

# Process each fasta file
for fasta_file in os.listdir(fasta_dir):
    if fasta_file.endswith(".fasta"):
        gene_name = os.path.splitext(fasta_file)[0]  # Extract gene name (e.g., atpA from atpA.fasta)
        
        # Check if the gene name exists in the table
        if gene_name in locus_gene_map:
            fasta_path = os.path.join(fasta_dir, fasta_file)
            output_path = os.path.join(output_dir, f"{gene_name}.filtered.fasta")
            
            # Read the fasta file
            sequences = list(SeqIO.parse(fasta_path, "fasta"))
            
            if len(sequences) == 1:
                # If only one sequence, convert to multiline fasta
                with open(output_path, "w") as output_handle:
                    SeqIO.write(sequences[0], output_handle, "fasta")
            else:
                # If more than one sequence, filter based on LOCUS
                locus = locus_gene_map[gene_name]
                filtered_sequences = [seq for seq in sequences if locus in seq.id]
                
                if filtered_sequences:
                    # Write the filtered sequence to the output file
                    with open(output_path, "w") as output_handle:
                        SeqIO.write(filtered_sequences[0], output_handle, "fasta")
                else:
                    print(f"No matching sequence found for {gene_name} with LOCUS {locus}. Skipping.")
        else:
            print(f"{gene_name} not found in the table. Skipping.")
