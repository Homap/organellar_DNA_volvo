from Bio import SeqIO
import sys

# Input FASTA file and sequence ID to remove
input_fasta = sys.argv[1] #"input.fasta"  # Replace with your input file
sequence_to_remove = sys.argv[2] #"seq1"  # Replace with the sequence ID to remove

# Read the FASTA file and filter out the sequence to remove
remaining_sequences = []
for record in SeqIO.parse(input_fasta, "fasta"):
    if record.id != sequence_to_remove:
        remaining_sequences.append(record)

# Print the remaining sequences in FASTA format
for record in remaining_sequences:
    print(f">{record.id}\n{record.seq}")
