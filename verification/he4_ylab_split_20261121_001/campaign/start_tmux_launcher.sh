#!/bin/bash
set -eu
cd /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer
exec /cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/python-3.8.13-myixjs7jesskpbhwiowibizezfh52rky/bin/python3 -u /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/lambda_he4_ylab_split_20261121_001/launch_once.py --config /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/lambda_he4_ylab_split_20261121_001/launch_config.json --submit > /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/lambda_he4_ylab_split_20261121_001/launch.log 2>&1
