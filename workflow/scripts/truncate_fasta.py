from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Define the input and output file paths
input_fasta = sys.argv[1]  # Input FASTA file
output_fasta = sys.argv[2]  # Output FASTA file

# Define the sequence ID and cutoff positions
sequence_id = sys.argv[3]  # ID of the sequence to modify
start_position = int(sys.argv[4])  # Start position for the region to retain (0-based index)
cutoff_position = int(sys.argv[5])  # End position for the region to retain (0-based index)

# Open the output file for writing
with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id == sequence_id:
            # Calculate the length of the replaced region
            removed_length = len(record.seq) - (cutoff_position - start_position)
            # Replace the removed bases with '-' and keep the retained region
            retained_seq = record.seq[start_position:cutoff_position]
            padded_seq = "-" * removed_length + str(retained_seq)
            record.seq = Seq(padded_seq)
        # Write the (modified or unmodified) sequence to the output file
        SeqIO.write(record, output_handle, "fasta")

print(f"Modified FASTA file saved to {output_fasta}")
