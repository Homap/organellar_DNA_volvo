import sys
# Define the file path
input_file = sys.argv[1]

# Reading the file and initializing variables
blocks = []
current_block = []
locus_name = None

# Reading and grouping LOCUS blocks
with open(input_file, 'r') as file:
    for line in file:
        if line.startswith("LOCUS"):
            if current_block:
                blocks.append((locus_name, current_block))
            locus_name = line.split()[1]
            current_block = [line.strip()]
        elif current_block is not None:
            current_block.append(line.strip())
    if current_block:
        blocks.append((locus_name, current_block))

# Refined logic to handle multi-line qualifiers in the FEATURES section
results = []

for locus, block in blocks:
    in_features = False
    gene_name = None
    coverage = ''
    match = ''
    info_line = ""  # Accumulate multiline qualifiers
    
    for line in block:
        # Identify the FEATURES section
        if line.startswith("FEATURES"):
            in_features = True
            continue

        # Stop processing after FEATURES if another section starts
        if in_features and (line.startswith("ORIGIN") or line.startswith("LOCUS")):
            in_features = False
            continue

        # Process lines within FEATURES
        if in_features:
            line = line.strip()
            if line.startswith("/gene="):
                gene_name = line.split('=')[1].strip('"')
            elif line.startswith("/info="):
                info_line = line.split('=')[1].strip('"')  # Start capturing multiline info
            elif info_line:  # Handle continuation lines for /info
                if line.startswith("/") or line.startswith("gene") or line == "":
                    # End of multiline /info
                    info_parts = info_line.split(',')
                    for part in info_parts:
                        if 'coverage' in part.lower():
                            coverage = part.split('=')[-1].strip().replace("coverage", "").replace('"', '').strip()
                        elif 'match' in part.lower():
                            match = part.split('=')[-1].strip().replace("match", "").replace('"', '').strip()
                    info_line = ""
                else:
                    info_line += " " + line.strip()

            if gene_name and (coverage or match):
                try:
                    results.append([locus, gene_name, float(coverage.strip('%')), float(match.strip('%'))])
                except ValueError:
                    print(f"Skipping invalid data for {gene_name} in {locus}: coverage={coverage}, match={match}")
                gene_name = None  # Reset for the next gene
                coverage, match = '', ''

# Create a DataFrame for the results
import pandas as pd
columns = ['LOCUS', 'GeneName', 'Coverage', 'Match']
df = pd.DataFrame(results, columns=columns)

# Group by GeneName to keep the one with higher coverage and match
sorted_df = df.sort_values(['GeneName', 'Coverage', 'Match'], ascending=[True, False, False])
main_output = sorted_df.drop_duplicates(subset=['GeneName'], keep='first')
secondary_output = sorted_df[~sorted_df.index.isin(main_output.index)]

# Save or display the results
main_output.to_csv(open(sys.argv[2], "w"), sep='\t', index=False)
secondary_output.to_csv(open(sys.argv[3], "w"), sep='\t', index=False)

print("Main output (higher coverage and match):")
print(main_output)
print("\nSecondary output (lower coverage and match):")
print(secondary_output)
