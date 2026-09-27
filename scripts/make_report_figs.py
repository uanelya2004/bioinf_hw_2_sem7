import gzip, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,
 "axes.edgecolor":"#8a8983","axes.labelcolor":"#52514e","xtick.color":"#52514e","ytick.color":"#52514e",
 "axes.titlesize":10,"axes.titleweight":"bold","axes.titlecolor":"#0b0b0b","axes.titlelocation":"left"})
BLUE="#2a78d6"; ORANGE="#eb6834"; GRID="#e6e5e0"; INK="#0b0b0b"; INK2="#52514e"
O="report/fig/"; H="./"
import os; os.makedirs(O, exist_ok=True)
st=["75ep","1-2ss","5-6ss","14ss"]
chip=[23.86,12.82,18.85,41.55]; inp=[39.93,33.46,57.05,27.79]
peaks=[36,40,839,170]; frip=[0.0172,0.0140,0.236,0.162]
mats={"75ep":"results/figures/foxd3_75ep.matrix.gz","1-2ss":"results/figures/foxd3_1-2ss.matrix.gz",
      "5-6ss":"teacher/5-6ss/figures/foxd3_5-6ss.matrix.gz","14ss":"results/figures/foxd3_14ss.matrix.gz"}
prof={}
for s,p in mats.items():
    with gzip.open(H+p,'rt') as f:
        f.readline(); rows=[l.rstrip("\n").split('\t') for l in f]
    m=np.array([[float(x) if x!='nan' else np.nan for x in r[6:]] for r in rows])
    nw=np.array([r[0].startswith("NW") for r in rows])
    prof[s]=(np.nanmean(m[:,:200],0),np.nanmean(m[:,200:],0),m,nw)
ratio=[prof[s][0][98:102].mean()/prof[s][1][98:102].mean() for s in st]
print("ratio",ratio)

def lab(ax,bars,fmt):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2,b.get_height(),fmt(b.get_height()),ha="center",va="bottom",fontsize=8,color=INK)
def style(ax):
    ax.yaxis.grid(True,color=GRID,lw=0.8); ax.set_axisbelow(True); ax.tick_params(length=0)

# Fig 1: overview 1x4
fig,axs=plt.subplots(1,4,figsize=(11,2.9))
x=np.arange(4); w=0.38
b1=axs[0].bar(x-w/2-0.01,chip,w,color=BLUE,label="ChIP"); b2=axs[0].bar(x+w/2+0.01,inp,w,color=ORANGE,label="Input")
lab(axs[0],b1,lambda v:f"{v:.1f}"); lab(axs[0],b2,lambda v:f"{v:.1f}")
axs[0].set_title("Глубина, млн выравниваний"); axs[0].legend(frameon=False,fontsize=8,loc="upper left")
b=axs[1].bar(x,peaks,0.6,color=BLUE); lab(axs[1],b,lambda v:f"{v:.0f}"); axs[1].set_title("Число пиков (q < 0.01)")
b=axs[2].bar(x,frip,0.6,color=BLUE); lab(axs[2],b,lambda v:f"{v:.3f}"); axs[2].set_title("FRiP, %")
b=axs[3].bar(x,ratio,0.6,color=BLUE); lab(axs[3],b,lambda v:f"{v:.1f}"); axs[3].set_title("ChIP / Input в summit")
for ax in axs:
    style(ax); ax.set_xticks(x); ax.set_xticklabels(st)
    ax.set_ylim(0,ax.get_ylim()[1]*1.12)
fig.tight_layout(); fig.savefig(O+"fig_overview.png",dpi=220)

# Fig 2: average profiles, shared y
fig,axs=plt.subplots(1,4,figsize=(11,2.6),sharey=True)
xx=np.linspace(-2000,2000,200)
for ax,s in zip(axs,st):
    c,i,_,_=prof[s]
    ax.plot(xx,c,color=BLUE,lw=2,label="Foxd3 ChIP"); ax.plot(xx,i,color=ORANGE,lw=2,label="Input")
    ax.set_title(f"{s}  (n = {prof[s][2].shape[0]})"); style(ax); ax.set_xticks([-2000,0,2000]); ax.set_xticklabels(["−2 kb","summit","+2 kb"])
axs[0].set_ylabel("CPM (среднее по пикам)"); axs[0].legend(frameon=False,fontsize=8)
fig.tight_layout(); fig.savefig(O+"fig_profiles.png",dpi=220)

# Fig 3: overlap matrices
pair={("75ep","1-2ss"):(12,12),("75ep","5-6ss"):(22,22),("75ep","14ss"):(18,18),("1-2ss","5-6ss"):(30,30),("1-2ss","14ss"):(12,12),("5-6ss","14ss"):(40,40)}
jac={("75ep","1-2ss"):.1737,("75ep","5-6ss"):.0227,("75ep","14ss"):.0751,("1-2ss","5-6ss"):.0281,("1-2ss","14ss"):.0426,("5-6ss","14ss"):.0369}
n=dict(zip(st,peaks))
D=np.full((4,4),np.nan); J=np.full((4,4),np.nan)
for (a,b),(ab,ba) in pair.items():
    i,j=st.index(a),st.index(b); D[i,j]=100*ab/n[a]; D[j,i]=100*ba/n[b]; J[i,j]=J[j,i]=jac[(a,b)]
from matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list("b",["#eef4fc","#9ec3ef","#2a78d6","#0d3a73"])
fig,axs=plt.subplots(1,2,figsize=(10,3.9))
for ax,M,t,fmt in [(axs[0],D,"а) % пиков строки, задевающих пики столбца",lambda v:f"{v:.0f}%"),(axs[1],J,"б) Индекс Жаккара (bp)",lambda v:f"{v:.3f}")]:
    im=ax.imshow(np.ma.masked_invalid(M),cmap=cmap,vmin=0,vmax=np.nanmax(M))
    ax.set_xticks(range(4)); ax.set_xticklabels(st); ax.set_yticks(range(4)); ax.set_yticklabels(st)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    for i in range(4):
        for j in range(4):
            v=M[i,j]
            if np.isnan(v): ax.text(j,i,"—",ha="center",va="center",color="#8a8983")
            else: ax.text(j,i,fmt(v),ha="center",va="center",fontsize=9,color="white" if v>0.55*np.nanmax(M) else INK)
    ax.set_title(t,fontsize=9); fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04).outline.set_visible(False)
fig.tight_layout(); fig.savefig(O+"fig_overlap.png",dpi=220)

# Fig 4: peak location composition (chromosomes vs unplaced scaffolds NW_)
nw=[2,1,15,52]; chrom=[p-q for p,q in zip(peaks,nw)]
fig,ax=plt.subplots(figsize=(5.2,2.6))
y=np.arange(4)[::-1]
ax.barh(y,[100*c/p for c,p in zip(chrom,peaks)],color=BLUE,height=0.6,label="хромосомы 1–25")
ax.barh(y,[100*q/p for q,p in zip(nw,peaks)],left=[100*c/p+0.6 for c,p in zip(chrom,peaks)],color=ORANGE,height=0.6,label="неразмещённые скаффолды NW_")
for yy,q,p in zip(y,nw,peaks): ax.text(101,yy,f"{q}/{p} ({100*q/p:.0f}%)",va="center",fontsize=8,color=INK)
ax.set_yticks(y); ax.set_yticklabels(st); ax.set_xlim(0,125); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel("% пиков")
ax.xaxis.grid(True,color=GRID); ax.set_axisbelow(True); ax.tick_params(length=0)
ax.legend(frameon=False,fontsize=8,loc="upper center",bbox_to_anchor=(0.45,1.25),ncol=2)
fig.tight_layout(); fig.savefig(O+"fig_scaffolds.png",dpi=220)
print("ok")
