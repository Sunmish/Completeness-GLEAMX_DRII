#!/bin/bash -l

# Generate list of random RA and Dec positions in specified region of sky

#SBATCH --account=mwasci
#SBATCH --partition=workq
#SBATCH --clusters=garrawarla
#SBATCH --nodes=1
#SBATCH --output=/astro/mwasci/software/kross/GLEAM-X-pipeline/log_garrawarla/generate_pos.o%A
#SBATCH --error=/astro/mwasci/software/kross/GLEAM-X-pipeline/log_garrawarla/generate_pos.e%A
#SBATCH --export=all

# echo "Reminder on this branch this is not needed. Exiting. "
# exit 1
set -x 
module load singularity
source /astro/mwasci/software/kross/GLEAM-X-pipeline/GLEAM-X-pipeline-garrawarla.profile
# echo $SINGULARITY_BINDPATH
export containerImage=/astro/mwasci/kross/GLEAM-X-pipeline/gleamx_container.img

start_time=$(date +%s)

# Read input parameters
if [[ $1 ]] && [[ $2 ]] && [[ $3 ]] && [[ $4 ]]; then
    # Number of random source positions to simulate
    nsrc=$(echo $1 | awk -F"=" '{print $NF}')
    # Region of sky - enter ra_min, ra_max, dec_min, dec_max, in deg.
    # For example, for 60 < RA < 300 enter: 60,300,-90,90
    # For RA < 60 or RA > 300, and -40 < Dec < -10, enter: 300,60,-40,-10
    region=$(echo $2 | awk -F"=" '{print $NF}')
    # Minimum separation between simulated sources, in arcmin
    sep_min=$(echo $3 | awk -F"=" '{print $NF}')
    # Output directory
    output_dir=$(echo $4 | awk -F"=" '{print $NF}')
else
    echo "Give me: nsrc region sep_min output_dir"
    exit 1
fi

# Create output directory
if [ -e $output_dir ]; then
    echo "Error: Output directory $output_dir already exists. Aborting."
    exit 1
else
    mkdir $output_dir
    cd $output_dir || exit 1
fi

# Write input parameters to file for record
cat >> input_parameters_generate_pos.txt <<EOPAR
nsrc = $nsrc
region = $region
sep_min = $sep_min
output_dir = $output_dir
EOPAR

# Run Python script to generate RA and Dec positions
srun singularity exec -B "/astro/mwasci/software/kross/Completeness-GLEAMX_DRII/" $containerImage /astro/mwasci/software/kross/Completeness-GLEAMX_DRII/generate_pos.py --nsrc=$nsrc --region=$region --sep_min=$sep_min source_pos.txt

end_time=$(date +%s)
duration=$(echo "$end_time-$start_time" | bc -l)
echo "Total runtime = $duration sec"

# Move output and error files to output directory
# root=/astro/mwasci/$USER/generate_pos
# mv $root.o${SLURM_JOB_ID} $root.e${SLURM_JOB_ID} .

exit 0
