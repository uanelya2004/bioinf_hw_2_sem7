#!/usr/bin/env bash
set -euo pipefail
export PATH=/home/STUDY/FBMF/bioinformatics/chipseq/software/seminar1-chipseq/bin:$PATH
cd ~/chipseq_seminar
for s in 14ss 75ep 1-2ss; do bash scripts/run_stage.sh $s; done
bash scripts/compare_stages.sh
echo "ALL DONE"
