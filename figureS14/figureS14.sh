#!/usr/bin/env bash
if [[ $# -ne 2 ]]; then
    echo "2 Parameters must be provided."
    echo "Usage: figureS14.sh path_to_repository path_to_miptools_container"
    exit
fi

repo=$1
container=$2

analysis=$repo'/figureS14'
base=$repo'/source_data/miptools_data/base_resources'
proj=$repo'/source_data/miptools_data/ideel-barcode/'
work=$repo'/source_data/'
singularity exec -B $analysis:/opt/analysis -B $proj/:/opt/project_resources -B $base/:/opt/resources -B $work:/opt/work $container python /opt/analysis/figureS14.py
