from collections import defaultdict
import sys


def read_fasta(file_path):
    """
    Reads a FASTA file and stores sequences in a dictionary.
    If multiple sequences have the same key, they are stored as a list.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with headers as keys and sequences as values (list of sequences).
    """
    fasta_dict = defaultdict(list)
    current_key = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_key = line[1:]  # Remove '>' from header
            else:
                fasta_dict[current_key].append(line)

    # Combine sequences split across multiple lines
    for key in fasta_dict:
        fasta_dict[key] = [''.join(seq) for seq in fasta_dict[key]]

    return fasta_dict

def print_sequences(fasta_dict):
    """
    Prints the longest sequence for keys with multiple sequences,
    and prints the sequence for keys with a single sequence.

    Args:
        fasta_dict (dict): A dictionary with headers as keys and sequences as values (list of sequences).
    """
    for key, sequences in fasta_dict.items():
        if len(sequences) > 1:
            longest_sequence = max(sequences, key=len)
            print(f">{key}\n{longest_sequence}")
        else:
            print(f">{key}\n{sequences[0]}")

if __name__ == "__main__":
    fasta_file = sys.argv[1]  # Replace with the path to your FASTA file
    fasta_dict = read_fasta(fasta_file)
    print_sequences(fasta_dict)