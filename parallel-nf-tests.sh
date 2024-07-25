#!/bin/bash -l

#SBATCH --job-name=hofflab-phageannotator
##SBATCH --mail-type=<status>
##SBATCH --mail-user=<email>

#SBATCH --account=pedslabs
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15GB
##SBATCH --gpus=<type:quantity>
#SBATCH --time=04:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator
#SBATCH --export=all
#SBATCH --output=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/%A_%a.outerr # where STDOUT goes
#SBATCH --error=/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/HoffLab-phageannotator/nf-test-parallel/%A_%a.outerr # where STDERR goes

#SBATCH --array=0-20%21

# activate nf-core environment
micromamba activate nf-core

# export tmp workdirectory for nf-test
NFT_WORKDIR=/tmp/.nf-test_${SLURM_ARRAY_TASK_ID}/
export NFT_WORKDIR

# make an array of nf-test files
nf_tests=(
# CSV_SRADOWNLOAD_FETCHNGS
    "/modules/local/fetchngs/sra_ids_to_runinfo/tests/main.nf.test" \
    "/modules/local/fetchngs/sra_runinfo_to_ftp/tests/main.nf.test" \
    "/modules/local/fetchngs/sra_fastq_ftp/tests/main.nf.test" \
    # FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS
        "/modules/nf-core/custom/sratoolsncbisettings/tests/main.nf.test" \
        "/modules/nf-core/sratools/prefetch/tests/main.nf.test" \
        "/modules/nf-core/sratools/fasterqdump/tests/main.nf.test" \
        "/subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/tests/main.nf.test" \
    "/modules/local/fetchngs/asperacli/tests/main.nf.test" \
"/subworkflows/local/csv_sradownload_fetchngs/tests/main.nf.test" \
# READ_PREPROCESS


)

echo $(pwd)
echo /mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/VTDB-createdb${nf_tests[$SLURM_ARRAY_TASK_ID]}

# run nf-test on each file in the array
nf-test \
test \
/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/VTDB-createdb${nf_tests[$SLURM_ARRAY_TASK_ID]} \
--verbose \
--tap /mmfs1/gscratch/pedslabs_hoffman/carsonjm/CarsonJM/VTDB-createdb/nf-test-parallel/${SLURM_ARRAY_TASK_ID}.tap \
--profile +singularity \
# --updateSnapshot \
# --cleanSnapshot
