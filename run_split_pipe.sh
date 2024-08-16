#script saved on /cephfs2/hannas/tetraploids/analysis/run_split_pipe.sh
#make executable: chmod +x run_split_pipe.sh


#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

# Run split-pipe commands
split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-1_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-1_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_1

split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-2_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-2_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_2

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-3_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-3_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_3 &

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-4_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-4_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_4 &

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-5_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-5_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_5 &

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-6_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-6_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_6 &

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-7_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-7_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_7 &

nohup split-pipe --mode dge --kit WT --chemistry v2 \
    --genome_dir /public/genomics/species_references/nextflow/Genome_References/Ensembl/Ensembl_Parse/homo_sapiens/GRCh38/Release_102/Parse_index/homo_sapiens.GRCh38.dna.102.Parse_index \
    --fq1 $path2data/expdata/UDI-WT-8_1.fastq.gz \
    --fq2 $path2data/expdata/UDI-WT-8_2.fastq.gz \
    --output_dir $path2data/analysis/UDI_WT_8 &

split-pipe --mode comb \
    --sublibraries $path2data/analysis/UDI_WT_1 $path2data/analysis/UDI_WT_2 $path2data/analysis/UDI_WT_3 $path2data/analysis/UDI_WT_4 $path2data/analysis/UDI_WT_5 $path2data/analysis/UDI_WT_6 $path2data/analysis/UDI_WT_7 $path2data/analysis/UDI_WT_8 \
    --output_dir $path2data/analysis/combined
