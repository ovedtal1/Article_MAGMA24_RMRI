#!/bin/bash

bash Data/download.sh

(
    cd Fig3_TSE_seq_def_pulseq
    python3 ./run.py
)

(
    cd Fig4_TSE_2Dre-implementation
    python3 ./run.py
)

(
    cd Fig6_reproducible_recon
    bash ./run.sh
)

(
    cd Fig7_masking
    python3 ./run.py
)
