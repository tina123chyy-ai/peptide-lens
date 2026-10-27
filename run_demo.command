#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
bash build.sh
mkdir -p results
./build/peptide-lens --normal data/normal.fasta --mutant data/mutant.fasta --out results/demo --sample demo_single
printf '\n完成 / Done. Open the results folder for CSV tables.\n'
if [[ "$(uname -s)" == Darwin ]]; then open results; fi
