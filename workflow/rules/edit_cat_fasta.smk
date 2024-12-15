rule edit_cat_fasta:
    input:
        "results/cpDNA_concat/cpDNA.done.txt",
        "results/mtDNA_concat/mtDNA.done.txt"
    output:
        cpdna_cat="results/cpDNA_concat/cpDNA.nameedit.done.txt",
        mtdna_cat="results/mtDNA_concat/mtDNA.nameedit.done.txt"
    shell:
        """
        for gene in results/cpDNA_concat/*.temp.fasta
        do
        prefix=$(echo $gene | cut -f3 -d "/" | cut -f1 -d ".")
        python workflow/scripts/modify_fasta_names.py ${{gene}} > results/cpDNA_concat/${{prefix}}.edited.fasta
        done
        touch {output.cpdna_cat}
        for gene in results/mtDNA_concat/*.temp.fasta
        do
        prefix=$(echo $gene | cut -f3 -d "/" | cut -f1 -d ".")
        python workflow/scripts/modify_fasta_names.py ${{gene}} > results/mtDNA_concat/${{prefix}}.edited.fasta
        done
        touch {output.mtdna_cat}
        """