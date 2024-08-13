conda activate spipe 

#making human reference genome with spipe v2.3.1
nohup split-pipe \
--mode mkref \
--genome_name GRCh38 \
--fasta /data1/hania/MultiSpeciesComp/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes /data1/hania/MultiSpeciesComp/newvolume/genomes/Homo_sapiens.GRCh38.112.gtf.gz \
--output_dir $path2data/genomes/GRCh38 &

path2data="/data1/hania/Tetraploids/newvolume"
split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-1_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-1_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_1

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-2_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-2_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_2

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-3_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-3_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_3

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-4_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-4_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_4

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-5_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-5_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_5

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-6_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-6_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_6

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-7_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-7_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_7

split-pipe --mode all --kit WT --chemistry v2 --genome_dir /data1/hania/MultiSpeciesComp/newvolume/genomes/GRCh38/ \
--fq1 $path2data/expdata/UDI-WT-8_R1_001.fastq.gz \
--fq2 $path2data/expdata/UDI-WT-8_R2_001.fastq.gz \
--output_dir $path2data/analysis/UDI_WT_8
