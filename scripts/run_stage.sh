#!/usr/bin/env bash
# Полный ChIP-seq workflow для одной стадии.
# Запуск: bash scripts/run_stage.sh 14ss   (или 75ep, 1-2ss)
set -euo pipefail

STAGE=$1
THREADS=4
GSIZE=1360000000          # effective genome size как на семинаре

cd ~/chipseq_seminar
CHIP=data/$STAGE/foxd3_${STAGE}_chip.nuclear.bam
INPUT=data/$STAGE/foxd3_${STAGE}_input.nuclear.bam
QC=results/qc; SIG=results/signal; PK=results/peaks; FIG=results/figures
mkdir -p $QC $SIG $PK $FIG logs
log(){ echo "[$(date +%H:%M:%S)] $STAGE: $*"; }

# 1. BAM QC
log "1/7 samtools flagstat"
samtools flagstat -@ $THREADS $CHIP  > $QC/foxd3_${STAGE}_chip.flagstat.txt
samtools flagstat -@ $THREADS $INPUT > $QC/foxd3_${STAGE}_input.flagstat.txt

# 2. CPM BigWig (пропускаем, если уже посчитано)
for S in chip input; do
  BAM=data/$STAGE/foxd3_${STAGE}_${S}.nuclear.bam
  BW=$SIG/foxd3_${STAGE}_${S}.CPM.bw
  if [ ! -s $BW ]; then
    log "2/7 bamCoverage $S"
    bamCoverage -b $BAM -o $BW --normalizeUsing CPM --binSize 10 \
      --numberOfProcessors $THREADS 2> logs/bamCoverage_${STAGE}_${S}.log
  fi
done

# 3. log2(ChIP/Input)
if [ ! -s $SIG/foxd3_${STAGE}.log2_IP_Input.bw ]; then
  log "3/7 bigwigCompare"
  bigwigCompare -b1 $SIG/foxd3_${STAGE}_chip.CPM.bw -b2 $SIG/foxd3_${STAGE}_input.CPM.bw \
    --operation log2 --pseudocount 1 --binSize 10 --numberOfProcessors $THREADS \
    -o $SIG/foxd3_${STAGE}.log2_IP_Input.bw 2> logs/bigwigCompare_${STAGE}.log
fi

# 4. MACS3 peak calling
NP=$PK/foxd3_${STAGE}_peaks.narrowPeak
if [ ! -s $NP ]; then
  log "4/7 MACS3 callpeak"
  macs3 callpeak -t $CHIP -c $INPUT -f BAMPE -g $GSIZE -n foxd3_${STAGE} \
    -q 0.01 --outdir $PK 2> logs/macs3_${STAGE}.log
fi
NPEAKS=$(wc -l < $NP)
log "   peaks: $NPEAKS"

# 5. Самый значимый peak (сортировка по колонке 9 = -log10 q)
log "5/7 top peaks"
sort -k9,9gr $NP | awk 'NR<=10' > $PK/foxd3_${STAGE}_top10.narrowPeak
head -1 $PK/foxd3_${STAGE}_top10.narrowPeak | awk -v OFS='\t' '
  { x=$9; e=int(x); if (e<x) e++; m=10^(e-x);
    print "coords","signalValue","-log10p","-log10q","q_value";
    printf "%s:%d-%d\t%s\t%s\t%s\t%.3fe-%d\n", $1,$2,$3,$7,$8,$9,m,e }' \
  > $QC/foxd3_${STAGE}_top_peak.tsv
cat $QC/foxd3_${STAGE}_top_peak.tsv

# 6. FRiP
log "6/7 FRiP"
CHIP_TOTAL=$(samtools view -@ $THREADS -c $CHIP)
INPUT_TOTAL=$(samtools view -@ $THREADS -c $INPUT)
IN_PEAKS=$(samtools view -@ $THREADS -c -L $NP $CHIP)
FRIP=$(awk -v a=$IN_PEAKS -v b=$CHIP_TOTAL 'BEGIN{printf "%.6f", a/b}')
printf "sample\treads_in_peaks\ttotal_reads\tFRiP\nfoxd3_%s_chip\t%s\t%s\t%s\n" \
  $STAGE $IN_PEAKS $CHIP_TOTAL $FRIP > $QC/foxd3_${STAGE}.frip.tsv
MAXQ=$(awk 'NR==1||$9>m{m=$9} END{print m}' $NP)
printf "stage\tchip_alignments\tinput_alignments\tpeaks\tFRiP\tmax_log10q\n%s\t%s\t%s\t%s\t%s\t%s\n" \
  $STAGE $CHIP_TOTAL $INPUT_TOTAL $NPEAKS $FRIP $MAXQ > $QC/foxd3_${STAGE}.summary.tsv
cat $QC/foxd3_${STAGE}.summary.tsv

# 7. Heatmap + profile вокруг summit
log "7/7 computeMatrix / plotHeatmap / plotProfile"
computeMatrix reference-point --referencePoint center \
  -R $PK/foxd3_${STAGE}_summits.bed \
  -S $SIG/foxd3_${STAGE}_chip.CPM.bw $SIG/foxd3_${STAGE}_input.CPM.bw \
  -b 2000 -a 2000 --binSize 20 -p $THREADS \
  -o $FIG/foxd3_${STAGE}.matrix.gz 2> logs/computeMatrix_${STAGE}.log
for T in pdf png; do
  plotHeatmap -m $FIG/foxd3_${STAGE}.matrix.gz -out $FIG/foxd3_${STAGE}.heatmap.$T \
    --samplesLabel "Foxd3 ChIP" "Input" --refPointLabel "summit" \
    --plotTitle "Foxd3 $STAGE" 2>> logs/plots_${STAGE}.log
  plotProfile -m $FIG/foxd3_${STAGE}.matrix.gz -out $FIG/foxd3_${STAGE}.profile.$T \
    --samplesLabel "Foxd3 ChIP" "Input" --refPointLabel "summit" \
    --plotTitle "Foxd3 $STAGE" --perGroup 2>> logs/plots_${STAGE}.log
done
log "DONE"
