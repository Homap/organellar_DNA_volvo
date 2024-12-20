from Bio import SeqIO
import sys

# Define the input and output file paths
input_fasta = sys.argv[1] #"input.fasta"  # Replace with your input FASTA file
output_fasta = sys.argv[2] #"output_trimmed.fasta"  # Replace with your desired output FASTA file

# Define the trimming range
trim_start = int(sys.argv[3]) #297  # Remove bases from 1 to 297 (0-based index)

# Open the output file for writing
with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Trim the sequence by removing the first 297 bases
        record.seq = record.seq[trim_start:]
        # Write the trimmed sequence to the output file
        SeqIO.write(record, output_handle, "fasta")

print(f"Trimmed FASTA file saved to {output_fasta}")
