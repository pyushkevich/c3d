#!/bin/bash
# Make the Convert3D logo (convert3d.jpg) and a smaller version without captions
# (convert3d-card.jpg) from the outputs of real c3d commands. See README.md.
set -e
if [[ $# -ne 3 ]]; then
  echo "usage: $0 MRIcrop-orig.gipl.gz work_dir Barlow-SemiBold.ttf"
  exit 1
fi
T=$1; P=$2; FONT=$3
D=$(cd $(dirname $0) && pwd)
mkdir -p $P

# the input and the outputs of the c3d commands shown in the logo
c3d $T -o $P/original.nii.gz
c3d $T $T -resample 20% -interpolation NearestNeighbor -reslice-identity -o $P/resample.nii.gz
c3d $T -smooth 2mm -o $P/smooth.nii.gz
c3d $T -smooth 1vox -grad -foreach -dup -times -endfor -accum -add -endaccum -sqrt -o $P/grad.nii.gz
c3d $T -smooth 1vox -otsu 3 -o $P/otsu.nii.gz
c3d $T -smooth 1vox -otsu 3 -thresh 0 0 1 0 -sdt -o $P/sdt.nii.gz

python3 $D/compose_logo.py $P $FONT $P/convert3d.jpg
python3 $D/compose_logo.py $P $FONT $P/convert3d-card.jpg thumb
