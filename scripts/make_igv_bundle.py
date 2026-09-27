#!/usr/bin/env python3
"""Собирает маленький пакет для IGV-Web:
- вырезает из BigWig только окна ±30 kb вокруг выбранных локусов;
- переименовывает хромосомы RefSeq (NC_007112.7 ...) в UCSC (chr1 ...),
  чтобы файлы совпадали со встроенным геномом Zebrafish GRCz11/danRer11 в IGV-Web;
- переименовывает хромосомы в narrowPeak.
Запуск: python3 scripts/make_igv_bundle.py   (из ~/chipseq_seminar)
"""
import os, pyBigWig

# RefSeq -> UCSC: NC_007112.7 = chr1 ... NC_007136.7 = chr25 (длины совпадают с danRer11)
REF2UCSC = {f"NC_{7112+i:06d}.7": f"chr{i+1}" for i in range(25)}
REF2UCSC["NC_002333.2"] = "chrM"

FLANK = 30000
LOCI = {  # имя: (RefSeq chrom, start, end)
    "L1_strong_5-6ss_chr3":      ("NC_007114.7", 14914145, 14914787),
    "L2_shared_all4_chr2":       ("NC_007113.7", 33275999, 33276609),
    "L2b_shared_all4_chr6":      ("NC_007117.7", 32206483, 32207024),
    "L3_specific_5-6ss_chr1":    ("NC_007112.7",   259734,   260218),
    "L3b_specific_5-6ss_chr17":  ("NC_007128.7", 27446917, 27447336),
    "L3c_biased_14ss_chr15":     ("NC_007126.7", 35583433, 35584048),
}

T = "teacher/5-6ss/signal"
TRACKS = {
    "75ep_chip.CPM":    "results/signal/foxd3_75ep_chip.CPM.bw",
    "75ep_input.CPM":   "results/signal/foxd3_75ep_input.CPM.bw",
    "75ep_log2":        "results/signal/foxd3_75ep.log2_IP_Input.bw",
    "1-2ss_chip.CPM":   "results/signal/foxd3_1-2ss_chip.CPM.bw",
    "1-2ss_input.CPM":  "results/signal/foxd3_1-2ss_input.CPM.bw",
    "1-2ss_log2":       "results/signal/foxd3_1-2ss.log2_IP_Input.bw",
    "5-6ss_chip.CPM":   f"{T}/foxd3_5-6ss_chip.nuclear.CPM.bw",
    "5-6ss_input.CPM":  f"{T}/foxd3_5-6ss_input.nuclear.CPM.bw",
    "5-6ss_log2":       f"{T}/foxd3_5-6ss_chip.log2_IP_Input.nuclear.bw",
    "14ss_chip.CPM":    "results/signal/foxd3_14ss_chip.CPM.bw",
    "14ss_input.CPM":   "results/signal/foxd3_14ss_input.CPM.bw",
    "14ss_log2":        "results/signal/foxd3_14ss.log2_IP_Input.bw",
}
PEAKS = {
    "75ep":  "results/peaks/foxd3_75ep_peaks.narrowPeak",
    "1-2ss": "results/peaks/foxd3_1-2ss_peaks.narrowPeak",
    "5-6ss": "teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.narrowPeak",
    "14ss":  "results/peaks/foxd3_14ss_peaks.narrowPeak",
}

OUT = "igv_bundle"
os.makedirs(OUT, exist_ok=True)

# окна по хромосомам (отсортированы, слиты)
wins = {}
for c, s, e in LOCI.values():
    wins.setdefault(c, []).append((max(0, s - FLANK), e + FLANK))
for c in wins:
    w = sorted(wins[c]); m = [list(w[0])]
    for a, b in w[1:]:
        if a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    wins[c] = m
chrom_order = sorted(wins, key=lambda c: int(REF2UCSC[c][3:]))

for name, path in TRACKS.items():
    src = pyBigWig.open(path)
    sizes = src.chroms()
    dst = pyBigWig.open(f"{OUT}/{name}.bw", "w")
    dst.addHeader([(REF2UCSC[c], sizes[c]) for c in chrom_order], maxZooms=10)
    for c in chrom_order:
        for a, b in wins[c]:
            iv = src.intervals(c, a, min(b, sizes[c])) or []
            if not iv: continue
            dst.addEntries([REF2UCSC[c]] * len(iv), [x[0] for x in iv],
                           ends=[x[1] for x in iv], values=[float(x[2]) for x in iv])
    dst.close(); src.close()
    print("bw  ", name)

for st, path in PEAKS.items():
    with open(path) as f, open(f"{OUT}/peaks_{st}.narrowPeak", "w") as o:
        for line in f:
            p = line.split("\t")
            if p[0] in REF2UCSC:
                p[0] = REF2UCSC[p[0]]; o.write("\t".join(p))
    print("peak", st)

with open(f"{OUT}/loci.txt", "w") as o:
    for n, (c, s, e) in LOCI.items():
        o.write(f"{n}\t{REF2UCSC[c]}:{s-5000}-{e+5000}\t(peak {c}:{s}-{e})\n")
print(open(f"{OUT}/loci.txt").read())
