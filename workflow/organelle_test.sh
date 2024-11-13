# Run Get organelle
#P26503
# Add the extra species as Maria has to the tree

# datadirP26503='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P26503_filtered_reads'

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP26503}/${sample_name}.ucseqs_1.fastq.gz
# fastq2=${datadirP26503}/${sample_name}.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_run.sh ../data/cpDNA_genome/cpDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_cp.label.fasta ${fastq1} ${fastq2} ../result/cpDNA_get_organelle/${sample_name}_cpDNA
# done < sample_names_P26503.txt

# datadirP26503='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P26503_filtered_reads'

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP26503}/${sample_name}.ucseqs_1.fastq.gz
# fastq2=${datadirP26503}/${sample_name}.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_mt_run.sh ../data/mtDNA_genome/mtDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_mt.label.fasta ${fastq1} ${fastq2} ../result/mtDNA_get_organelle/${sample_name}_mtDNA
# done < sample_names_P26503.txt

# datadirP28566='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P28566_filtered_reads'

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP28566}/${sample_name}.ucseqs_1.fastq.gz
# fastq2=${datadirP28566}/${sample_name}.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_run.sh ../data/cpDNA_genome/cpDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_cp.label.fasta ${fastq1} ${fastq2} ../result/cpDNA_get_organelle/${sample_name}_cpDNA
# done < sample_names_P28566.txt

# datadirP28566='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P28566_filtered_reads'

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP28566}/${sample_name}.ucseqs_1.fastq.gz
# fastq2=${datadirP28566}/${sample_name}.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_mt_run.sh ../data/mtDNA_genome/mtDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_mt.label.fasta ${fastq1} ${fastq2} ../result/mtDNA_get_organelle/${sample_name}_mtDNA
# done < sample_names_P28566.txt

# datadirmerged='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/merged_reads'
# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirmerged}/${sample_name}.R1.fastq.gz
# fastq2=${datadirmerged}/${sample_name}.R2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_run.sh ../data/cpDNA_genome/cpDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_cp.label.fasta ${fastq1} ${fastq2} ../result/cpDNA_get_organelle/${sample_name}_cpDNA
# done < merged.txt

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirmerged}/${sample_name}.R1.fastq.gz
# fastq2=${datadirmerged}/${sample_name}.R2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_mt_run.sh ../data/mtDNA_genome/mtDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_mt.label.fasta ${fastq1} ${fastq2} ../result/mtDNA_get_organelle/${sample_name}_mtDNA
# done < merged.txt

# # P29912
# datadirP29912='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P29912_filtered_reads'
# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP29912}/${sample_name}.merged.ucseqs_1.fastq.gz
# fastq2=${datadirP29912}/${sample_name}.merged.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_run.sh ../data/cpDNA_genome/cpDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_cp.label.fasta ${fastq1} ${fastq2} ../result/cpDNA_get_organelle/${sample_name}_cpDNA
# done < sample_names_P29912.txt

# while read -r sample_name                                                  
# do                                                                       
# echo ${sample_name}     
# fastq1=${datadirP29912}/${sample_name}.merged.ucseqs_1.fastq.gz
# fastq2=${datadirP29912}/${sample_name}.merged.ucseqs_2.fastq.gz
# echo ${fastq1} ${fastq2} 
# sbatch get_organelle_mt_run.sh ../data/mtDNA_genome/mtDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_mt.label.fasta ${fastq1} ${fastq2} ../result/mtDNA_get_organelle/${sample_name}_mtDNA
# done < sample_names_P29912.txt

# # 
# fastq1='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P26503_filtered_reads/P26503_122_S174.ucseqs_1.fastq.gz'
# fastq2='/proj/snic2022-23-81/algae_phylo_volv/data/soft_links/filtered_reads/P26503_filtered_reads/P26503_122_S174.ucseqs_2.fastq.gz'
# sbatch get_organelle_mt_run.sh ../data/mtDNA_genome/mtDNA_chlorophyceae_complete_genome.fasta ../data/get_organelle_label_dataset/chlorophyceae_mt.label.fasta ${fastq1} ${fastq2} ../result/mtDNA_get_organelle/P26503_122_S174_mtDNA

# Create contiguous scaffolds from contigs
for i in *
do
echo $i
mkdir ../cpDNA_output/$i
cp $i/*1.path_sequence.fasta ../cpDNA_output/$i
done

cd ../cpDNA_output/
for i in *
do
echo $i
mv $i/*1.path_sequence.fasta ${i}.fasta
done

rmdir *cpDNA

# For chloroplast DNA (do it based on order and orientation of C. reinhardtii)
# mkdir -p result/scaffolding/cpDNA result/scaffolding/mtDNA

# cd result/cpDNA_output
# for i in *
# do
# echo $i
# prefix=$(echo $i | cut -f1 -d ".")
# ragtag.py scaffold ../../data/cpDNA_genome/reinhardtii.fa $i -o ../scaffolding/cpDNA/$prefix
# cd ../scaffolding/cpDNA/$prefix 
# echo "nucmer"
# nucmer ../../../../data/cpDNA_genome/reinhardtii.fa ragtag.scaffold.fasta -p $prefix
# echo "mummer"
# mummerplot -l ${prefix}.delta --png -p $prefix
# dnadiff -d ${prefix}.delta -p $prefix
# mv ragtag.scaffold.fasta ${prefix}.ragtag.scaf.fasta
# cd ../../../cpDNA_output
# done

# Gather all cpDNA assemblies in one folder
# mkdir result/scaffolding/cpDNA/scaffolded_cpDNA
# cd result/scaffolding/cpDNA
# for i in *
# do
# echo $i
# cp $i/${i}.ragtag.scaf.fasta scaffolded_cpDNA
# done

# Transfer to computer under /Users/homapapoli/Dropbox/Algae/Projects/Algae_phylo_comparative/result/organelle_DNA
scp -r homap@rackham.uppmax.uu.se:/proj/snic2022-23-81/cpDNA_assembly_phylo_volv/result/scaffolding/cpDNA/scaffolded_cpDNA .

# Submit the concatenated fasta to GeSeq
# Used all available Chlorophyceae chloroplast for annotation
# Copied results to Uppmax
scp -r geseq_cpDNA_annotation homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/cpDNA_assembly/result

# Collect all annotations and change the fasta header to gene name and sample name
for i in *gb
do
echo $i
samplename=$(grep 'LOCUS' $i | awk '{print $2}')
echo $samplename
mv $i ${samplename}.gb
done
# Convert gene bank format to CDS
for i in *gb
do
echo $i
samplename=$(echo $i | cut -f1 -d ".")
get_annotated_regions_from_gb.py $i -o ${samplename}_cds -t CDS
done

# Concatenate all samples per gene
# /proj/snic2022-23-81/cpDNA_assembly_phylo_volv/result/geseq_cpDNA_annotation/gb
while read -r genename                                                 
do  
echo $genename
cat *cds/gene/${genename}.fasta > ../../concat_per_gene_cpDNA/${genename}.fasta
done < cpDNA_gene_set.txt

# module load bioinfo-tools cufflinks/2.2.1
# mkdir cpDNA_all_species
# mkdir mtDNA_all_species
# gffread -w Chlamydomonas_incerta_cpDNA_cds.fasta -g chlamydomonas_incerta_cpDNA.fasta  GeSeqJob-20240601-103750_MW465979.1_GFF3.gff3
# tr -d '\n' < Chlamydomonas_schloesseri_cpDNA.fasta > tt.fasta 
# mv tt.fasta Chlamydomonas_schloesseri_cpDNA.fasta
# gffread -w Chlamydomonas_schloesseri_cpDNA_cds.fasta -g Chlamydomonas_schloesseri_cpDNA.fasta GeSeqJob-20240601-104504_MW465980.1_GFF3.gff3
# Names in Haematococcus lacustris, re-annotated using GeSeq
# gffread -w Haematococcus_lacustris_cpDNA_cds.fasta -g Haematococcus_lacustris.fasta GeSeqJob-20240604-100807_NC_037007.1_GFF3.gff3

cd /proj/snic2022-23-81/cpDNA_assembly_phylo_volv/result/cpDNA_all_species/cpDNA_ncbi
mkdir concantenated_cpDNA
# Concatenate all genes from the NCBI
while read -r genename                                                 
do  
echo $genename
cat *_${genename}.fasta >> concantenated_cpDNA/${genename}.fasta
done < ../cpDNA_gene_set.txt

# Concatenate all NCBI genes with the concatenated sample genes
for i in ../../../concat_per_gene_cpDNA/*fasta
do
prefix=$(echo $i | cut -f5 -d "/" | cut -f1 -d ".")
echo $prefix
cat $i ${prefix}.fasta > ${prefix}.all.fasta
done

# Edit sequence names. Keep original sample IDs since it is easlier to track.
cd /proj/snic2022-23-81/cpDNA_assembly_phylo_volv/result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA
for gene in *.all.fasta
do
prefix=$(echo $gene | cut -f1 -d ".")
echo $prefix
python ../../../../code/modify_fasta_names.py ${prefix}.all.fasta > ${prefix}.all.edited.fasta
done

rm -f *all.fasta

# Do one alignment per gene
module load bioinfo-tools MAFFT/7.407
../result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA/aligned
for input_fasta in ../result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA/*all.edited.fasta
do
genename=$(echo $input_fasta | cut -f6 -d "/" | cut -f1 -d ".")
echo $genename
sbatch mafft_alignment.sh ${input_fasta} ../result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA/aligned/${genename}.aligned.fasta
done

# Cleaning the alignment by removing common issues such as 
# gaps, divergent sequences, large insertions and deletions and poorly aligned sequence ends can substantially improve analyses. 
# Manual editing of MSAs is very widespread but is time-consuming and difficult to reproduce.

# Edit each gene alignment separately before concatenation. Then concatenate after alignment. Then there is no need to edit because each edited 
# alignment is now one block of the concatenated alignment. View each gene in alignment viewer and find the problematic genes.
# Plus, these are protein coding genes, check that the translation of the codons gives meaningful results and if you are adding non-coding DNA sequences, check those as well. 
### Remove excessive gaps with CIAlign
# module load bioinfo-tools biopython/1.80-py3.10.8

# module load conda/latest biopython/1.80-py3.10.8
# export CONDA_ENVS_PATH=/proj/naiss2023-23-145/algae_conda
# export CONDA_PKGS_DIRS=/proj/naiss2023-23-145/algae_conda/algae_phylo_pkgs
# conda activate algae_phylo_env

# CIAlign --infile cpDNA_concatented.fasta.fasta --outfile_stem cpDNA_concatented.fasta.cialign.fasta --insertion_min_perc 0.9 --remove_insertions --crop_ends_redefine_perc 0.5 --crop_ends --remove_divergent_minperc 0.7 --remove_divergent

# ### Remove poorly fit sequences with SequenceBouncer
# python SequenceBouncer/SequenceBouncer.py -i cpDNA_concatented.fasta.cialign.fasta_cleaned.fasta -o cpDNA_concatented.fasta.cialign.bouncer.fasta

# CIAlign --infile cpDNA_concatented.fasta.cialign.bouncer.fasta_output_clean.fasta --outfile_stem cpDNA_concatented.fasta.cialign.bouncer.cialign.fasta --crop_ends_redefine_perc 0.5 --crop_ends


# I manually edited the alignments.

# Run ModelFinder in iqtree for each gene to find the best model for each. 
module load bioinfo-tools iqtree/2.2.2.6-omp-mpi gcc/9.3.0 openmpi/3.1.5
for i in *aligned.fasta
do
iqtree2 -s cpDNA_concatented.fasta.cialign.fasta_cleaned.fasta -m MFP
done
# Partition model
# When you concatenate the genes, each gene might have a different rate of evolution. And it is not clear what the partitions are, for example, two genes concatenated could form 1 parition since they
# have similar evolutionary rates, etc. You need to run a partition finder in order to find the right paritions. 
# nexus: This is the nexus file. You can run modelfinder on each gene to find the best model and then based on that, make the partition file. 
begin sets;
    charset part1 = 1-100;
    charset part2 = 101-384;
    charpartition mine = HKY+G:part1, GTR+I+G:part2;
end;

iqtree -s example.phy -p example.nex

# Also a program to find paritions

# Bootstrapping : Ultrafast bootstrap approximation
iqtree -s example.phy -m TIM2+I+G -B 1000
# Standard non-parametric bootstrap
iqtree -s example.phy -m TIM2+I+G -b 100




# Check alignments by eye and discard those with obvious problems

# Concatenate alignments
sed -i 's/-.*//g' *.aligned.fasta | sed -i 's/_//g'
# Run from /proj/snic2020-16-269/private/homap/cpDNA_assembly/result/aligned_gene_cpDNA
python ../../code/concat_alignment.py ./ cpDNA_concatented.fasta 3 DNA

# Clean fasta and nexus files
sed -i 's/n/-/g' *fasta
sed -i 's/?/-/g' cpDNA_concatenated.fasta.nex
sed -i 's/?/-/g' cpDNA_concatenated.fasta.phy

# Run IQtree
module load bioinfo-tools iqtree/2.2.2.6-omp-mpi gcc/9.3.0 openmpi/3.1.5
iqtree2 -nt 4 -s ../aligned_gene_cpDNA/cpDNA_concatented.fasta.nex
sed 's/.aligned.fasta//g' cpDNA_concatenated.fasta.nex > cpDNA_concatenated.fasta.edited.nex
# Find best partition scheme followed by tree inference and bootstrap:
iqtree2 -s cpDNA_concatenated.fasta.edited.nex -p partition_file.nex -m MFP --merge -T 8 -b 1000

# Run iqtree2 to find the best partitioning scheme
# sbatch iqtree2.sh 
# The best partition scheme is under partition_file.nex.best_scheme.nex
# Now Running IQtree with the best partition scheme and bootstrap
iqtree2 -s cpDNA_concatenated.fasta.edited.phy -p partition_file.nex.best_scheme.nex -b 1000


#!/usr/bin/bash

#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 04:00:00
#SBATCH -J iqtree2-bootstrap
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools iqtree/2.2.2.6-omp-mpi gcc/9.3.0 openmpi/3.1.5

iqtree2 -s $1 -p $2 -T 16 -b 20 --prefix $3

for i in {1..50}
do
echo $i
sbatch ../../../../../code/iqtree2_bootstrap.sh cpDNA_concatenated.fasta.edited.phy partition_file.nex.best_scheme.nex cpDNA_boot.${i}
done

iqtree2 -s cpDNA_concatenated.fasta.edited.phy -p partition_file.nex.best_scheme.nex -T 16 -b 1000

iqtree2 -s ../result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA/aligned_edited/cpDNA_concatenated.fasta.edited.nex -p ../result/cpDNA_all_species/cpDNA_ncbi/concantenated_cpDNA/aligned_edited/partition_file.nex -m MFP --merge -T 16



# Open tree output in FigTree
/proj/snic2020-16-269/private/homap/cpDNA_assembly/result/aligned_gene_cpDNA/cpDNA_concatented.fasta.nex.treefile

# For mitohondrion DNA
# Concatenate scaffolds with 20 N in between if the assembly is not complete
for i in *
do
echo $i
mkdir ../mtDNA_output/$i
cp $i/*1.path_sequence.fasta ../mtDNA_output/$i
done

cd ../mtDNA_output/
for i in *
do
echo $i
mv $i/*1.path_sequence.fasta ${i}.fasta
done

rmdir *mtDNA

# For mtDNA (do it based on order and orientation of C. reinhardtii)
cd result/mtDNA_output
module load bioinfo-tools MUMmer
for i in *
do
echo $i
prefix=$(echo $i | cut -f1 -d ".")
ragtag.py scaffold ../../data/mtDNA_genome/reinhardtii.fa $i -o ../scaffolding/mtDNA/$prefix
cd ../scaffolding/mtDNA/$prefix 
echo "nucmer"
nucmer ../../../../data/mtDNA_genome/reinhardtii.fa ragtag.scaffold.fasta -p $prefix
echo "mummer"
mummerplot -l ${prefix}.delta --png -p $prefix
dnadiff -d ${prefix}.delta -p $prefix
mv ragtag.scaffold.fasta ${prefix}.ragtag.scaf.fasta
cd ../../../mtDNA_output
done

# Gather all cpDNA assemblies in one folder
mkdir result/scaffolding/mtDNA/scaffolded_mtDNA
cd result/scaffolding/mtDNA
for i in *
do
echo $i
cp $i/${i}.ragtag.scaf.fasta scaffolded_mtDNA
done

# Transfer to computer under /Users/homapapoli/Dropbox/Algae/Projects/Algae_phylo_comparative/result/organelle_DNA
scp -r homap@rackham.uppmax.uu.se:/proj/snic2022-23-81/cpDNA_assembly_phylo_volv/result/scaffolding/cpDNA/scaffolded_cpDNA .





for i in ../result/mtDNA_output/*
do
samplename=$(echo $i | cut -f4 -d "/" | cut -f1 -d ".")
echo $samplename
python stitch_scaffolds.py $i ${samplename} ../result/mtDNA_output/stitched/${samplename}_stitched.fasta
done

# Transfer to computer under /Users/homapapoli/Dropbox/Algae/Projects/Algae_phylo_comparative/result/organelle_DNA/stitched
scp -r homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/cpDNA_assembly/result/mtDNA_output/stitched .

# Submit the concatenated fasta to GeSeq
# Used all available Chlorophyceae chloroplast for annotation
# Copied results to Uppmax
scp -r geseq_cpDNA_annotation homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/cpDNA_assembly/result

# Collect all annotations and change the fasta header to gene name and sample name
for i in *gb
do
echo $i
samplename=$(grep 'LOCUS' $i | awk '{print $2}')
echo $samplename
mv $i ${samplename}.gb
done
# Convert gene bank format to CDS
for i in *gb
do
echo $i
samplename=$(echo $i | cut -f1 -d ".")
get_annotated_regions_from_gb.py $i -o ${samplename}_cds -t CDS
done

# Concatenate all samples per gene
while read -r genename                                                 
do  
echo $genename
cat *cds/gene/${genename}.fasta > ../concat_per_gene_mtDNA/${genename}.fasta
done < mtDNA_gene_set.txt

# Do one alignment per gene
module load bioinfo-tools MAFFT/7.407

for input_fasta in *fasta
do
genename=$(echo $input_fasta | cut -f1 -d ".")
echo $genename
mafft-linsi ${input_fasta} > ${genename}.aligned.fasta
done

# Check alignments by eye and discard those with obvious problems

# Concatenate alignments
# Run from /proj/snic2020-16-269/private/homap/cpDNA_assembly/result/aligned_gene_cpDNA
python ../../code/concat_alignment.py ./ mtDNA_concatented.fasta 3 DNA

# Clean fasta and nexus files
sed -i 's/n/-/g' *nex
sed -i 's/?/-/g' *nex

# Run IQtree
module load bioinfo-tools iqtree/2.2.2.6-omp-mpi gcc/9.3.0 openmpi/3.1.5
iqtree2 -nt 4 -s ../aligned_gene_mtDNA/mtDNA_concatented.fasta.nex

# Open tree output in FigTree
/proj/snic2020-16-269/private/homap/cpDNA_assembly/result/aligned_gene_mtDNA/mtDNA_concatented.fasta.nex.treefile

for i in {1..50}
do
echo $i
sbatch ../../../../../code/iqtree2_bootstrap.sh mtDNA_concantenated.fasta.phy partition.nex.best_scheme.nex mtDNA_boot.${i}
done

cat cpDNA*.boottrees > alltrees
iqtree2 -con -t alltrees
iqtree2 -T 8 -s cpDNA_concatenated.fasta.edited.phy -p partition_file.nex.best_scheme.nex -te alltrees.contree -pre alltrees.contree
mv alltrees.contree.treefile cpDNA.alltrees.contree.treefile

cat mtDNA*.boottrees > alltrees
iqtree2 -con -t alltrees
iqtree2 -T 8 -s mtDNA_concantenated.fasta.phy -p partition.nex.best_scheme.nex -te alltrees.contree -pre alltrees.contree
mv alltrees.contree.treefile mtDNA.alltrees.contree.treefile

# NOTE: P28566_112 is identical to P28566_106 but kept for subsequent analysis
# NOTE: P28566_114 is identical to P28566_110 but kept for subsequent analysis

# NOTE: P26503_118 is identical to P26503_117 but kept for subsequent analysis
# NOTE: P28566_107 is identical to P28566_103 but kept for subsequent analysis
# NOTE: P28566_112 is identical to P28566_106 but kept for subsequent analysis
# NOTE: P28566_108 (identical to P28566_103) is ignored but added at the end