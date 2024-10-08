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
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-1_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-1_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_1 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

#each of the following individually: 
split-pipe --mode all --kit WT --chemistry v2 \
     --genome_dir $path2data/genomes/GRCh38/ \
     --fq1 $path2data/expdata/raw_data/UDI-WT-2_1.fastq.gz \
     --fq2 $path2data/expdata/raw_data/UDI-WT-2_2.fastq.gz \
     --output_dir $path2data/expdata/data_samptab/UDI_WT_2 \
     --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-3_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-3_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_3 \ 
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-4_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-4_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_4 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-5_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-5_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_5 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-6_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-6_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_6 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-7_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-7_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_7 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids"

split-pipe --mode all --kit WT --chemistry v2 \
    --genome_dir $path2data/genomes/GRCh38/ \
    --fq1 $path2data/expdata/raw_data/UDI-WT-8_1.fastq.gz \
    --fq2 $path2data/expdata/raw_data/UDI-WT-8_2.fastq.gz \
    --output_dir $path2data/expdata/data_samptab/UDI_WT_8 \
    --samp_sltab $path2data/expdata/raw_data/Parse_WT_100K_Sample_Loading_Table.xlsm

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sbatch -J UDI_WT_1 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_1.sh
sbatch -J UDI_WT_2 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_2.sh
sbatch -J UDI_WT_3 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_3.sh
sbatch -J UDI_WT_4 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_4.sh
sbatch -J UDI_WT_5 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_5.sh
sbatch -J UDI_WT_6 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_6.sh
sbatch -J UDI_WT_7 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_7.sh
sbatch -J UDI_WT_8 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_8.sh


#each fo these are in seperate files. 

sbatch -J UDI_WT_1 -c 1 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_1.sh
sbatch -J UDI_WT_2 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_2.sh
sbatch -J UDI_WT_3 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_3.sh
sbatch -J UDI_WT_4 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_4.sh
sbatch -J UDI_WT_5 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_5.sh
sbatch -J UDI_WT_6 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_6.sh
sbatch -J UDI_WT_7 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_7.sh
sbatch -J UDI_WT_8 -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G UDI_WT_8.sh


#!/bin/bash

# Activate the conda environment
source /lmb/home/hannas/miniconda3/etc/profile.d/conda.sh
conda activate spipe

# Set path to data
path2data="/cephfs2/hannas/tetraploids/"

# Run split-pipe commands
split-pipe --mode comb \
    --sublibraries $path2data/expdata/data_samptab/UDI_WT_1 $path2data/expdata/data_samptab/UDI_WT_2 $path2data/expdata/data_samptab/UDI_WT_3 $path2data/expdata/data_samptab/UDI_WT_4 $path2data/expdata/data_samptab/UDI_WT_5 $path2data/expdata/data_samptab/UDI_WT_6 $path2data/expdata/data_samptab/UDI_WT_7 $path2data/expdata/data_samptab/UDI_WT_8 \
    --output_dir $path2data/expdata/data_samptab/combined

sbatch -J combine -c 2 --mail-type=ALL --mail-user=hannas@mrc-lmb.cam.ac.uk --mem=300G combine.sh


