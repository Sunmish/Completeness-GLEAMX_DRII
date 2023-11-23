#!/usr/bin/env bash

# A simple helper script to string together the completeness tasks that need to be run.

cluster=garrawarla
nsrc=90000
region=310,90,-90,30
flux=-3,-0.5,0.1
nfiles=6
sep_min=5
outdir=/astro/mwasci/kross/gleamx/GLEAMX_DRII/completeness_ims/
imageset='GLEAMX_DRII_170-231MHz'

imageset_dir=/astro/mwasci/kross/gleamx/GLEAMX_DRII/completeness_ims/

export GLEAMX="${outdir}"
export MYCODE=/astro/mwasci/software/kross/Completeness-GLEAMX_DRII/
export NCPUS=10
export CONTAINER=/astro/mwasci/kross/GLEAM-X-pipeline/gleamx_container.img

set -x

if [[ -z ${MYCODE} ]]
then
    echo "Error. Completeness code directory not set. Exiting. "
    return 1
fi

if [[ ! -d $outdir ]]
then
    echo "Making directory ${outdir}"
    mkdir -p "${outdir}"
fi

# mkdir "${GLEAMX}/input_images"

#TODO: See how well this works with symlinks. Need to be sure the container can follow them.
# for suffix in "" "_bkg" "_rms" "_projpsf_psf"
# do
#     if [[ -e "${imageset_dir}/${imageset}${suffix}.fits" ]]
#     then
#         cp -v "${imageset_dir}/${imageset}${suffix}.fits" "${GLEAMX}/input_images"
#     else
#         echo "Could not find ${imageset_dir}/${imageset}${suffix}.fits. Exiting. "
#         exit 1
#     fi
# done

# sbatch --time=06:00:00 --ntasks-per-node=1 $MYCODE/generate_pos.sh ${nsrc} ${region} 5 $GLEAMX/source_pos


"$MYCODE"/generate_fluxes.sh \
$nsrc \
$region \
$sep_min \
$flux \
$nfiles \
$outdir


if [[ $? -ne 0 ]]
then
    echo "Completeness simulation set up failed. Aborting."
    exit 1
fi

We will be blocking until we are finished
msg=($(sbatch \
    --array 1-$nfiles \
    --time 4:00:00 \
    --ntasks-per-node $NCPUS \
    --export ALL \
    -o "${outdir}/inject_source.o%A_%a" \
    -e "${outdir}/inject_source.e%A_%a" \
    "$MYCODE/inject_sources.sh" \
    "${GLEAMX}/input_images" \
    "${GLEAMX}/source_pos/source_pos.txt" \
    "${GLEAMX}/fluxes" \
    4.0 \
    "${GLEAMX}/inject" \
"${imageset}"))

jobid=${msg[3]}
# echo "$msg"
id=$(echo "$msg" | cut -d ' ' -f3)

msg=$(sbatch \
    --time 1:00:00 \
    --ntasks-per-node $NCPUS \
    --export ALL \
    -o "${outdir}/cmp_map.o%A" \
    -e "${outdir}/cmp_map.e%A" \
    "$MYCODE"/make_cmp_map.sh \
    "${GLEAMX}/source_pos/source_pos.txt" \
    "${GLEAMX}/inject" \
    "$flux" \
    "${GLEAMX}/input_images/${imageset}_projpsf_psf.fits" \
    "${region}" \
    6 \
"${GLEAMX}/results")

echo "$msg"
