rule concatenate_CDS:
    output:
        cpdna_cat="results/cpDNA_concat/cpDNA.done.txt",
        mtdna_cat="results/mtDNA_concat/mtDNA.done.txt"
    shell:
        """
        # Concatenate per gene for the sequenced cpDNA genomes
        mkdir -p results/cpDNA_concat
        while read -r genename                                                 
        do  
            echo "Processing $genename..."
            # Check if any matching file exists
            if ls results/cpDNA_CDS/*cds/gene/filtered/${{genename}}.filtered.fasta 1> /dev/null 2>&1; then
                # If files exist, concatenate them
                cat results/cpDNA_CDS/*cds/gene/filtered/${{genename}}.filtered.fasta >> results/cpDNA_concat/${{genename}}.temp.fasta
                echo "$genename has been concatenated."
            else
                # If no matching files, print a warning message
                echo "Warning: No files found for $genename. Skipping..."
            fi
        done < config/cpDNA_unique_gene_names.txt

        # Concatenate per gene for the NCBI cpDNA genomes
        while read -r genename
        do  
            echo "Processing $genename..."
            # Check if any matching file exists
            if ls data/cpDNA_genome/cpDNA_ncbi/*_${{genename}}.fasta 1> /dev/null 2>&1; then
                # If files exist, concatenate them
                cat data/cpDNA_genome/cpDNA_ncbi/*_${{genename}}.fasta >> results/cpDNA_concat/${{genename}}.temp.fasta
                echo "$genename has been concatenated."
            else
                # If no matching files, print a warning message
                echo "Warning: No files found for $genename. Skipping..."
            fi
        done < config/cpDNA_unique_gene_names.txt
        touch {output.cpdna_cat}

        # Concatenate per gene for the sequenced mtDNA genomes
        mkdir -p results/mtDNA_concat
        while read -r genename                                                 
        do  
            echo "Processing $genename..."
            # Check if any matching file exists
            if ls results/mtDNA_CDS/*cds/gene/filtered/${{genename}}.filtered.fasta 1> /dev/null 2>&1; then
                # If files exist, concatenate them
                cat results/mtDNA_CDS/*cds/gene/filtered/${{genename}}.filtered.fasta >> results/mtDNA_concat/${{genename}}.temp.fasta
                echo "$genename has been concatenated."
            else
                # If no matching files, print a warning message
                echo "Warning: No files found for $genename. Skipping..."
            fi
        done < config/mtDNA_unique_gene_names.txt

        # Concatenate per gene for the NCBI mtDNA genomes
        while read -r genename
        do  
            echo "Processing $genename..."
            # Check if any matching file exists
            if ls data/mtDNA_genome/mtDNA_ncbi/*_${{genename}}.fasta 1> /dev/null 2>&1; then
                # If files exist, concatenate them
                cat data/mtDNA_genome/mtDNA_ncbi/*_${{genename}}.fasta >> results/mtDNA_concat/${{genename}}.temp.fasta
                echo "$genename has been concatenated."
            else
                # If no matching files, print a warning message
                echo "Warning: No files found for $genename. Skipping..."
            fi
        done < config/mtDNA_unique_gene_names.txt
        touch {output.mtdna_cat}
        """
