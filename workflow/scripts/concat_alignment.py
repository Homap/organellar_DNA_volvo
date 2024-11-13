from __future__ import annotations
import glob
import sys
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio.Nexus import Nexus

# Written by Homa Papoli Yazdi September 2022
# The concat function was adapted from script by Etka Yapar.
# I learnt adding usage() and try() except() by looking at Etka Yapar's script.

def usage():
    msg = '''concatenate genewise alignments in fasta format into single NEXUS or PHYLIP formatted alignments

    Environment:
    This script requires Biopython packages. Before running, activate the relevant conda environment as follows:
    module load conda/latest
    export CONDA_ENVS_PATH=/proj/naiss2023-23-145/algae_conda
    export CONDA_PKGS_DIRS=/proj/naiss2023-23-145/algae_conda/algae_phylo_pkgs
    conda activate algae_phylo_env

    Usage:
    concat_alignment.py <alignment-input-dir> <output-file-name> <minimum_num_taxa> <molecule_type>

    alignment-input-dir: Path to directory containing alignment files. It should contain either protein or DNA. 
    Header of fasta sequences in alignment file must only contain species name.
    output-file-name: Name of the Nexus or Phylip output
    minimum_num_taxa: Minimum number of taxa for any given gene, e.g. 3
    molecule_type: protein/DNA

    A note about the input:
    input directory must contain either only protein or DNA sequences.
    The header of each alignment sequence in the fasta format must only contain the species name. 
    If headers are ASTGUB|12455AT3041, do the following:
    sed -i 's/|.*//g' *_aa.aln
    sed -i 's/|.*//g' *_nt.aln
    '''
    print(msg)

def check_taxa_number(nexi, tax_num=2):
    """
    Checks that at least a given number of taxa are present for each gene.
    """
    taxa_length = [len(nexi[i][1].taxlabels) for i in range(len(nexi))]

    gene_to_exclude = [taxa_length.index(i) for i in taxa_length if i < tax_num]
    new_nexi = [j for i, j in enumerate(nexi) if i not in gene_to_exclude]
    return new_nexi

def concat(gene_list, tax_num=2, molecule_type="DNA", input_dir="."):
    """Combine multiple nexus data matrices in one partitioned file.
    
    By default this will only work if the same taxa are present in each file 
    use same_taxa=False if you are not concerened by this"""
#re.sub("/|.*", "", gene.replace(input_dir, "").replace("_aa.aln", ""))
    # Read genes in a dictionary
    if molecule_type == "protein":
        gene_dict = {gene.replace(input_dir+"/", "").replace("_aa.aln", ""):AlignIO.read(gene, "fasta") for gene in gene_list}
    elif molecule_type == "DNA":
        gene_dict = {gene.replace(input_dir+"/", "").replace("_nt.aln", ""):AlignIO.read(gene, "fasta") for gene in gene_list}
    # Add feature to each gene
    [seq.annotations.update({"molecule_type":molecule_type})for gene in gene_dict for seq in gene_dict[gene]]

    # Convert each gene format into nexus 
    nexi = [(gene, Nexus.Nexus(format(gene_dict[gene], "nexus"))) for gene in gene_dict]

    new_nexi = check_taxa_number(nexi, tax_num)
    combined = Nexus.combine(new_nexi)
    return combined

def main():

    try:
        input_dir = sys.argv[1]
        output_name = sys.argv[2]
        mininum_num = int(sys.argv[3])
        molecule_type = sys.argv[4]
    except:
        usage()
        sys.exit()
    if len(sys.argv) < 5 or input_dir in ('-h', '--help'):
        usage()
        sys.exit()

    gene_list = [gene for gene in glob.glob(input_dir+"/*") if gene.endswith(".aln") or gene.endswith(".fa") or gene.endswith(".fasta")]
    combined = concat(gene_list, mininum_num, molecule_type, input_dir)
    # Write out the nexus file into output
    combined.write_nexus_data(filename=output_name+".nex")
    # Write out the phylip file
    combined.export_phylip(filename=output_name+".phy")
    print(f"Aren't we brilliant? YES WE ARE! We just made the computer write concatenated alignment from directory {input_dir} to files {output_name}.nex and {output_name}.phy. Well done to us!")
    return 0

if __name__ == "__main__":
    main()