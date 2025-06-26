#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --cpus-per-task 12 
#SBATCH --mem=24GB
#SBATCH --output=/scratch/mwasci/duchesst/processing/logs/%x.out
#SBATCH --export=ALL
#SBATCH --time=02:00:00
#SBATCH --account=pawsey1160

module load singularity/4.1.0-slurm

# A simple helper script to string together the completeness tasks that need to be run.

SOURCE_DENSITY=10
flux=-2.4,0.0,0.1
nfiles=1
sep_min=4
imageset_dir=./


export MYCODE=/software/projects/mwasci/duchesst/Completeness-GLEAMX_DRII/
export NCPUS=12
export CONTAINER=/software/projects/mwasci/duchesst/piip.img

set -x

TILE=$1
TILE_FILE=$2
params=$(more ${TILE_FILE} | grep "${TILE}")
ra=$(echo "$params" | awk -F',' '{print $2}')
dec=$(echo "$params" | awk -F',' '{print $3}')
ra_size=$(echo "$params" | awk -F',' '{print $4}')
dec_size=$(echo "$params" | awk -F',' '{print $5}')

nsrc=$(echo "${SOURCE_DENSITY}*${ra_size}*${dec_size}" | bc -l)
imageset="GLEAM300_${TILE}"
outdir="/scratch/mwasci/duchesst/processing/${TILE}/completeness"
export GLEAMX="${outdir}"
input="/scratch/mwasci/duchesst/processing/${TILE}/"
imageset_dir="./"

cd ${input} || exit 1

ra_min=$(echo "${ra}-(0.5*${ra_size})" | bc -l)
ra_max=$(echo "${ra}+(0.5*${ra_size})" | bc -l)
dec_min=$(echo "${dec}-(0.5*${dec_size})" | bc -l)
dec_max=$(echo "${dec}+(0.5*${dec_size})" | bc -l)

region="${ra_min},${ra_max},${dec_min},${dec_max}"

"$MYCODE/generate_fluxes.sh" \
$nsrc \
$region \
$sep_min \
$flux \
$nfiles \
"$outdir/"

"$MYCODE/inject_sources.sh" \
"${input}" \
"${GLEAMX}/source_pos/source_pos.txt" \
"${GLEAMX}/fluxes" \
4.0 \
"${GLEAMX}/inject" \
"${imageset}"

"$MYCODE"/make_cmp_map.sh \
"${GLEAMX}/source_pos/source_pos.txt" \
"${GLEAMX}/inject" \
"$flux" \
"${input}/${imageset}_psfmap.fits" \
"${region}" \
6 \
"${GLEAMX}/results"
