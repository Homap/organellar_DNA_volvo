rule find_gene_sets_cpDNA:
    output:
        "config/cpDNA_unique_gene_names.txt"
    run:
        import os
        import re

        # Directory containing the sample directories
        base_directory = "results/cpDNA_CDS"

        # Set to store unique GeneNames
        unique_gene_names = set()

        # Walk through all sample directories
        for root, dirs, files in os.walk(base_directory):
            for file in files:
                # Check if the file ends with ".filtered.fasta"
                if file.endswith(".filtered.fasta"):
                    # Extract the GeneName from the file name
                    match = re.match(r"(.+)\.filtered\.fasta", file)
                    if match:
                        unique_gene_names.add(match.group(1))

        # Convert the set to a sorted list
        unique_gene_names_sorted = sorted(unique_gene_names)

        # Write the unique GeneNames to a text file
        output_file = "config/cpDNA_unique_gene_names.txt"
        with open(output_file, "w") as f:
            for gene_name in unique_gene_names_sorted:
                f.write(f"{gene_name}\n")

        print(f"Unique GeneNames have been saved to {output_file}")
