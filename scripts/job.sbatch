#!/bin/bash
#SBATCH -p E5-2630-TITAN_X
#SBATCH -J foxd3_chip
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH -o logs/slurm_%j.log
cd ~/chipseq_seminar
bash scripts/run_all.sh
