#open the terminal
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe 

path2data="/cephfs2/hannas/tetraploids"
nohup split-pipe \
--mode mkref \
--genome_name GRCh38 \
--fasta $path2data/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes $path2data/genomes/Homo_sapiens.GRCh38.112.gtf.gz \
--output_dir $path2data/genomes/GRCh38 &

#making bash script: 

path2data="/cephfs2/hannas/tetraploids"
split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-1_1.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-1_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_1

split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-2_1.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-2_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_2

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-3_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-3_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_3 &

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-4_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-4_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_4 &

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-5_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-5_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_5 &

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-6_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-6_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_6 &

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-7_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-7_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_7 &

nohup  split-pipe --mode dge --kit WT --chemistry v2 --genome_dir $path2data/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-8_1.fastq.gz  \
--fq2 $path2data/expdata/UDI-WT-8_2.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_8 &

split-pipe \
    --mode comb \
    --sublibraries $path2data/analysis/UDI_WT_1 $path2data/analysis/UDI_WT_2 $path2data/analysis/UDI_WT_3 $path2data/analysis/UDI_WT_4 $path2data/analysis/UDI_WT_5 $path2data/analysis/UDI_WT_6 $path2data/analysis/UDI_WT_7 $path2data/analysis/UDI_WT_8 \
    --output_dir $path2data/analysis/combined

    
