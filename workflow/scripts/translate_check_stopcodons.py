from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import random


def pad_to_multiple_of_three(sequence):
    remainder = len(sequence) % 3
    if remainder != 0:
        sequence += 'N' * (3 - remainder)
    return sequence


def translate_sequence(sequence, frame):
    sliced_seq = sequence[frame:]
    padded_seq = pad_to_multiple_of_three(sliced_seq)
    return padded_seq.translate()


def count_internal_stop_codons(translated_seq):
    return translated_seq[:-1].count('*')  # Exclude last codon


def first_stop_codon_position(translated_seq):
    stops = [i for i, codon in enumerate(translated_seq[:-1]) if codon == '*']
    if stops:
        relative_pos = round(stops[0] / len(translated_seq), 4)
        absolute_pos = stops[0] * 3
        return relative_pos, absolute_pos
    return None, None


def replace_internal_stops(sequence):
    stop_codons = {"TAA", "TAG", "TGA"}
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    result = []
    modified = False

    for codon in codons[:-1]:  # Exclude last codon
        if codon in stop_codons:
            result.append("NNN")
            modified = True
        else:
            result.append(str(codon))
    result.append(str(codons[-1]))  # Add last codon
    return Seq("".join(result)), modified


def trim_from_first_stop(sequence, frame, stop_position):
    return sequence[:frame + stop_position]


def process_fasta(input_fasta, output_report, output_trimmed, output_stop_codon_seqs):
    trimmed_records = []
    stop_codon_entries = []

    with open("warnings.log", "w") as warning_log, open(output_report, 'w') as report_file:
        # Write headers
        report_file.write("Sequence_Name\tFrame\tInternal_Stop_Codons\n")

        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_name = record.id
            sequence = record.seq
            sequence_length = len(sequence)

            stop_codons = []
            stop_positions = []
            for frame in range(3):
                translated_seq = translate_sequence(sequence, frame)
                stop_codons.append(count_internal_stop_codons(translated_seq))
                rel_pos, abs_pos = first_stop_codon_position(translated_seq)
                stop_positions.append((rel_pos, abs_pos))
                report_file.write(f"{seq_name}\tFrame {frame + 1}\t{stop_codons[-1]}\n")

            # Populate stop_codon_entries for all frames
            stop_codon_entries.append(
                f"{seq_name}\t{stop_codons[0]}\t{stop_positions[0][0] or 'NA'}\t{stop_positions[0][1] or 'NA'}\t"
                f"{stop_codons[1]}\t{stop_positions[1][0] or 'NA'}\t{stop_positions[1][1] or 'NA'}\t"
                f"{stop_codons[2]}\t{stop_positions[2][0] or 'NA'}\t{stop_positions[2][1] or 'NA'}"
            )

            # Handle case with 0 stop codons
            frames_with_no_stop = [i for i, count in enumerate(stop_codons) if count == 0]
            if frames_with_no_stop:
                if len(frames_with_no_stop) > 1:
                    target_frame = 0 #random.choice(frames_with_no_stop)
                    warning_log.write(f"Warning: '{seq_name}' in '{input_fasta}' - Multiple frames with 0 stop codons. "
                                      f"Frame {target_frame + 1} chosen randomly.\n")
                else:
                    target_frame = frames_with_no_stop[0]
                adjusted_seq = sequence[target_frame:]
                padded_seq = pad_to_multiple_of_three(adjusted_seq)
                if padded_seq[-3:] in {"TAA", "TAG", "TGA"}:
                    padded_seq = padded_seq[:-3]  # Remove terminal stop codon
                header = f"{seq_name} Trimmed-Terminal-stop (Frame {target_frame + 1})"
                trimmed_records.append(SeqRecord(padded_seq, id=seq_name, description=header))
                continue

            # Handle case with 1 stop codon
            frames_with_one_stop = [i for i, count in enumerate(stop_codons) if count == 1]
            if frames_with_one_stop:
                target_frame = frames_with_one_stop[0]
                sliced_seq = sequence[target_frame:]  # Adjust sequence for target frame
                padded_seq = pad_to_multiple_of_three(sliced_seq)  # Pad to multiple of 3

                # Replace internal stop codons with NNN
                repaired_seq, modified = replace_internal_stops(sliced_seq)

                # Check if the sequence ends with incomplete codons (e.g., last 1-2 bases)
                stop_codons = {"TAA", "TAG", "TGA"}
                last_codon = repaired_seq[-3:]  # Extract the last codon
                remaining_bases = len(sliced_seq) % 3  # Check number of trailing bases

                # Handle case of internal stop codon with remaining bases
                if remaining_bases in {1, 2} and last_codon in stop_codons:
                    repaired_seq = repaired_seq[:-3] + "NNN" + sliced_seq[-remaining_bases:]  # Replace stop codon, keep trailing bases
                    had_terminal_stop = False  # No terminal stop codon
                elif last_codon in stop_codons:
                    repaired_seq = repaired_seq[:-3]  # Remove terminal stop codon
                    had_terminal_stop = True
                else:
                    had_terminal_stop = False

                # Set the appropriate header
                header_suffix = "OneInternal" if modified and not had_terminal_stop else "OneInternal-Terminal-stop"
                header = f"{seq_name} {header_suffix} (Frame {target_frame + 1})"

                trimmed_records.append(SeqRecord(repaired_seq, id=seq_name, description=header))
                continue


            # Handle case with >1 stop codons in all frames
            if all(count > 1 for count in stop_codons):
                best_frame, best_pos = None, None
                for frame, (rel_pos, abs_pos) in enumerate(stop_positions):
                    # Look for frames where rel_pos >= 0.45
                    if rel_pos is not None and rel_pos >= 0.45:
                        if best_pos is None or rel_pos < stop_positions[best_frame][0]:
                            best_frame, best_pos = frame, abs_pos

                if best_frame is not None and best_pos is not None:
                    # Trim the sequence from the first stop codon in the best frame
                    trimmed_seq = trim_from_first_stop(sequence, best_frame, best_pos)
                    header = f"{seq_name} Trimmed-Frame{best_frame + 1}"
                    trimmed_records.append(SeqRecord(trimmed_seq, id=seq_name, description=header))
                else:
                    # No rel_pos >= 0.45 found in any frame; keep original sequence
                    header = f"{seq_name} Untrimmed_3stopcodons"
                    trimmed_records.append(SeqRecord(sequence, id=seq_name, description=header))


    # Write outputs
    with open(output_trimmed, 'w') as trimmed_file:
        SeqIO.write(trimmed_records, trimmed_file, "fasta")

    with open(output_stop_codon_seqs, 'w') as stop_codon_file:
        stop_codon_file.write("Sequence_Name\tFrame1_Stop_Codons\tFrame1_Relative_Pos\tFrame1_Absolute_Pos\t"
                              "Frame2_Stop_Codons\tFrame2_Relative_Pos\tFrame2_Absolute_Pos\t"
                              "Frame3_Stop_Codons\tFrame3_Relative_Pos\tFrame3_Absolute_Pos\n")
        stop_codon_file.write("\n".join(stop_codon_entries) + "\n")

    print(f"Warnings logged to 'warnings.log'.")
    print(f"Translation report written to {output_report}")
    print(f"Trimmed sequences written to {output_trimmed}")
    print(f"Stop codon sequences report written to {output_stop_codon_seqs}")


if __name__ == "__main__":
    input_fasta = sys.argv[1]
    output_report = sys.argv[2]
    output_trimmed = sys.argv[3]
    output_stop_codon_seqs = sys.argv[4]
    process_fasta(input_fasta, output_report, output_trimmed, output_stop_codon_seqs)
