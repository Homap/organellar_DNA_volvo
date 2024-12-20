from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# Input and output file paths
input_fasta = sys.argv[1]  # Input FASTA file
output_fasta = sys.argv[2]  # Output FASTA file
sample_list = sys.argv[3]  # File containing the list of all expected sample IDs

# Read the input FASTA and find the sequence length
records = list(SeqIO.parse(input_fasta, "fasta"))

# Assume all sequences have the same length; get the length from the first sequence
if records:
    sequence_length = len(records[0].seq)
else:
    raise ValueError("Input FASTA file is empty.")

# Collect existing sample IDs (use only the first part of the header)
existing_samples = {}
for record in records:
    sample_id = record.id.split()[0]  # Take the first part of the header
    record.id = sample_id  # Update record ID to be the simplified sample ID
    record.description = ""  # Clear the description to avoid duplication
    existing_samples[sample_id] = record

# Ensure sequences are uppercase
for record in existing_samples.values():
    record.seq = record.seq.upper()

# Load the sample list and ensure order
with open(sample_list, "r") as sample_file:
    sample_order = [line.strip() for line in sample_file]

# Prepare the final list of records in the order of the sample list
final_records = []
for sample in sample_order:
    if sample in existing_samples:
        final_records.append(existing_samples[sample])  # Use existing sequence
    else:
        # Create a sequence of '-' with the same length as existing sequences
        new_record = SeqRecord(Seq("-" * sequence_length), id=sample, description="")
        final_records.append(new_record)

# Write the final records to the output file
with open(output_fasta, "w") as output_handle:
    SeqIO.write(final_records, output_handle, "fasta")

print(f"Updated FASTA file saved to {output_fasta}")
