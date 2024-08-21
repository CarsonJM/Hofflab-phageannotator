#!/bin/bash -l

#SBATCH --job-name=parallel-nf-test
##SBATCH --mail-type=<status>
##SBATCH --mail-user=<email>

#SBATCH --account=pedslabs
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=15GB
##SBATCH --gpus=<type:quantity>
#SBATCH --time=02:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator
#SBATCH --export=all
#SBATCH --output=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/%A_%a.outerr # where STDOUT goes
#SBATCH --error=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/%A_%a.outerr # where STDERR goes

#SBATCH --array=1-5%5

# activate nf-core environment
micromamba activate nf-core

# export tmp workdirectory for nf-test
# mkdir -p /tmp/.nf-test_${SLURM_ARRAY_TASK_ID}/
# NFT_WORKDIR=/tmp/.nf-test_${SLURM_ARRAY_TASK_ID}/
# export NFT_WORKDIR

# run nf-test on each file in the array
nf-test \
test \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/local/viromeqc/install/main.nf \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/nf-core/fastp/main.nf \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/nf-core/fastqc/main.nf \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/nf-core/bowtie2/build/main.nf \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/modules/nf-core/cat/fastq/main.nf \
--related-tests \
--follow-dependencies \
--tap /mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/${SLURM_ARRAY_TASK_ID}.tap \
--profile +nfcore_hyak \
--silent \
--shard ${SLURM_ARRAY_TASK_ID}/5 \
# --updateSnapshot \
# --cleanSnapshot \
--log /mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/${SLURM_ARRAY_TASK_ID}.log
