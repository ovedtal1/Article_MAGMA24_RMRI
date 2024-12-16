#!/bin/bash
#Copyright 2024. TU Graz. Institute of Biomedical Imaging.
#Author: Moritz Blumenthal

set -eu
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )
cd $SCRIPT_DIR

wget https://zenodo.org/records/14497769/files/ksp_fully.hdr
wget https://zenodo.org/records/14497769/files/ksp_fully.cfl