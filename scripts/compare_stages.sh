#!/usr/bin/env bash
# Сравнение peaks четырёх стадий + сводная таблица.
# Запуск: bash scripts/compare_stages.sh
set -euo pipefail
cd ~/chipseq_seminar
show(){ if command -v column >/dev/null; then column -t "$1"; else cat "$1"; fi; }
OUT=results/stage_comparison; mkdir -p $OUT/beds
STAGES=(75ep 1-2ss 5-6ss 14ss)
declare -A NP=(
  [75ep]=results/peaks/foxd3_75ep_peaks.narrowPeak
  [1-2ss]=results/peaks/foxd3_1-2ss_peaks.narrowPeak
  [5-6ss]=teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.narrowPeak
  [14ss]=results/peaks/foxd3_14ss_peaks.narrowPeak )
declare -A XLS=(
  [75ep]=results/peaks/foxd3_75ep_peaks.xls
  [1-2ss]=results/peaks/foxd3_1-2ss_peaks.xls
  [5-6ss]=teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.xls
  [14ss]=results/peaks/foxd3_14ss_peaks.xls )

# отсортированные peaks (для intersect) и слитые интервалы (для jaccard)
for s in "${STAGES[@]}"; do
  sort -k1,1 -k2,2n ${NP[$s]} > $OUT/beds/$s.sorted.narrowPeak
  cut -f1-3 $OUT/beds/$s.sorted.narrowPeak | bedtools merge > $OUT/beds/$s.merged.bed
done

# ---- Сводная таблица 4 стадий ----
xls_val(){ grep -m1 "$2" "$1" 2>/dev/null | awk '{print $NF}' || echo NA; }
{
printf "stage\tchip_alignments\tinput_alignments\tchip_fragments_MACS\tinput_fragments_MACS\tpeaks\tFRiP\tmax_log10q\tmedian_width\tmedian_signalValue\n"
for s in "${STAGES[@]}"; do
  N=${NP[$s]}
  if [ "$s" = "5-6ss" ]; then
    CT=$(awk 'NR==2{print $3}' teacher/5-6ss/frip.tsv); IT=NA
    FR=$(awk 'NR==2{print $4}' teacher/5-6ss/frip.tsv)
  else
    read -r CT IT FR < <(awk 'NR==2{print $2,$3,$5}' results/qc/foxd3_${s}.summary.tsv)
  fi
  CF=$(xls_val ${XLS[$s]} "total fragments in treatment"); IF=$(xls_val ${XLS[$s]} "total fragments in control")
  NPK=$(wc -l < $N)
  MQ=$(awk 'NR==1||$9>m{m=$9} END{print m}' $N)
  MW=$(awk '{print $3-$2}' $N | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  MS=$(cut -f7 $N | sort -g | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $s $CT $IT ${CF:-NA} ${IF:-NA} $NPK $FR $MQ $MW $MS
done
} > $OUT/summary_4stages.tsv
show $OUT/summary_4stages.tsv

# ---- 6 попарных сравнений ----
printf "A\tB\tpeaks_A\tpeaks_B\tA_overlap_B\tB_overlap_A\tA_only\tB_only\tjaccard\n" > $OUT/pairwise.tsv
for ((i=0;i<4;i++)); do for ((j=i+1;j<4;j++)); do
  A=${STAGES[$i]}; B=${STAGES[$j]}
  FA=$OUT/beds/$A.sorted.narrowPeak; FB=$OUT/beds/$B.sorted.narrowPeak
  bedtools intersect -a $FA -b $FB -u > $OUT/shared_${A}_in_${B}.bed
  bedtools intersect -a $FA -b $FB -v > $OUT/only_${A}_vs_${B}.bed
  bedtools intersect -a $FB -b $FA -v > $OUT/only_${B}_vs_${A}.bed
  NA_=$(wc -l < $FA); NB_=$(wc -l < $FB)
  AB=$(wc -l < $OUT/shared_${A}_in_${B}.bed)
  BA=$(bedtools intersect -a $FB -b $FA -u | wc -l)
  AO=$(wc -l < $OUT/only_${A}_vs_${B}.bed); BO=$(wc -l < $OUT/only_${B}_vs_${A}.bed)
  J=$(bedtools jaccard -a $OUT/beds/$A.merged.bed -b $OUT/beds/$B.merged.bed | awk 'NR==2{print $3}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $A $B $NA_ $NB_ $AB $BA $AO $BO $J >> $OUT/pairwise.tsv
done; done
echo; show $OUT/pairwise.tsv

# ---- Jaccard матрица 4x4 ----
{ printf "stage"; for b in "${STAGES[@]}"; do printf "\t%s" $b; done; echo
  for a in "${STAGES[@]}"; do printf "%s" $a
    for b in "${STAGES[@]}"; do
      J=$(bedtools jaccard -a $OUT/beds/$a.merged.bed -b $OUT/beds/$b.merged.bed | awk 'NR==2{print $3}')
      printf "\t%s" $J; done; echo; done; } > $OUT/jaccard_matrix.tsv
echo; show $OUT/jaccard_matrix.tsv

# ---- Stage-specific peaks (нет пересечения ни с одной другой стадией) и общие для всех ----
printf "stage\tpeaks\tspecific\tshared_with_any_other\n" > $OUT/stage_specific.tsv
for s in "${STAGES[@]}"; do
  cat $(for o in "${STAGES[@]}"; do [ $o != $s ] && echo $OUT/beds/$o.merged.bed; done) \
    | sort -k1,1 -k2,2n | bedtools merge > $OUT/beds/others_$s.bed
  bedtools intersect -a $OUT/beds/$s.sorted.narrowPeak -b $OUT/beds/others_$s.bed -v > $OUT/specific_$s.narrowPeak
  N=$(wc -l < $OUT/beds/$s.sorted.narrowPeak); SP=$(wc -l < $OUT/specific_$s.narrowPeak)
  printf "%s\t%s\t%s\t%s\n" $s $N $SP $((N-SP)) >> $OUT/stage_specific.tsv
done
cp $OUT/beds/5-6ss.sorted.narrowPeak $OUT/core_all4.narrowPeak
for s in 75ep 1-2ss 14ss; do
  bedtools intersect -a $OUT/core_all4.narrowPeak -b $OUT/beds/$s.sorted.narrowPeak -u > $OUT/tmp && mv $OUT/tmp $OUT/core_all4.narrowPeak
done
echo; show $OUT/stage_specific.tsv
echo "5-6ss peaks present at all 4 stages: $(wc -l < $OUT/core_all4.narrowPeak)"

# ---- Графики ----
python3 - <<'PY'
import pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out="results/stage_comparison"
m=pd.read_csv(f"{out}/jaccard_matrix.tsv",sep="\t",index_col=0)
fig,ax=plt.subplots(figsize=(5,4.2))
im=ax.imshow(m.values,cmap="Blues",vmin=0,vmax=max(0.01,m.values[~(m.values==1)].max()))
ax.set_xticks(range(4)); ax.set_xticklabels(m.columns); ax.set_yticks(range(4)); ax.set_yticklabels(m.index)
for i in range(4):
  for j in range(4): ax.text(j,i,f"{m.values[i,j]:.3f}",ha="center",va="center",fontsize=9)
ax.set_title("Foxd3 peaks: Jaccard index"); fig.colorbar(im,ax=ax); fig.tight_layout()
fig.savefig(f"{out}/jaccard_heatmap.png",dpi=200); fig.savefig(f"{out}/jaccard_heatmap.pdf")
s=pd.read_csv(f"{out}/summary_4stages.tsv",sep="\t")
fig,axs=plt.subplots(1,3,figsize=(11,3.3))
for ax,col,t in zip(axs,["peaks","FRiP","median_signalValue"],["Number of peaks","FRiP","Median fold-enrichment"]):
  ax.bar(s["stage"],pd.to_numeric(s[col],errors="coerce"),color="#4C72B0"); ax.set_title(t)
fig.tight_layout(); fig.savefig(f"{out}/stage_summary_bars.png",dpi=200)
print("plots saved")
PY
