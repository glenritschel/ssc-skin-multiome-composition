import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np

SSC='#D97706'; HC='#3B6FD4'; INK='#1c1c1a'; MUT='#6b6b66'; GRID='#e3e3e0'; SURF='#fcfcfb'
D=[('SSC1',1296,5761,.647),('SSC2',274,3444,.496),('SSC3',139,2269,.282),('SSC4',5,2760,.079),
   ('SSC5',572,3820,.357),('SSC6',12,1880,.074),('SSC7',89,267,.685),('SSC8',170,2801,.197),
   ('SSC9',419,6835,.354),('SSC10',239,5000,.215),
   ('HC1',3435,4662,.982),('HC2',2061,3993,.886),('HC3',1423,5684,.685),('HC4',1027,5923,.293)]
D.sort(key=lambda r:r[1])
names=[r[0] for r in D]; nfib=np.array([r[1] for r in D],float)
share=np.array([r[1]/r[2] for r in D]); col=np.array([r[3] for r in D])
c=[HC if n.startswith('HC') else SSC for n in names]; y=np.arange(len(D))

fig,ax=plt.subplots(1,2,figsize=(11.4,5.4),facecolor=SURF,gridspec_kw={'width_ratios':[1.2,1]})
for a in ax: a.set_facecolor(SURF)

a=ax[0]
a.hlines(y,3,nfib,color=c,lw=2,zorder=2)
a.scatter(nfib,y,s=46,color=c,zorder=3,edgecolor=SURF,linewidth=1.4)
a.axvline(200,color=INK,lw=1.2,ls='--',zorder=4)
a.text(210,-1.15,'floor = 200, fixed 09-Sep before donor 2',fontsize=8,color=INK,va='center')
for i,(n,s) in enumerate(zip(nfib,share)):
    a.text(max(n*1.28, 300), i, f'{int(n):,}  ({s:.0%})', va='center', fontsize=8, color=MUT)
a.set_xscale('log'); a.set_xlim(3,11000)
a.set_xticks([10,100,1000,10000]); a.set_xticklabels(['10','100','1,000','10,000'],fontsize=9,color=MUT)
a.set_yticks(y); a.set_yticklabels(names,fontsize=9,color=INK); a.set_ylim(-2,len(D)-.3)
a.set_xlabel('fibroblast nuclei recovered  (log scale)',fontsize=9.5,color=MUT)
a.set_title('Five of ten dcSSc donors clear the floor',fontsize=11.5,color=INK,loc='left',pad=10)
a.xaxis.grid(True,color=GRID,lw=.8); a.set_axisbelow(True)
for s in ('top','right','left'): a.spines[s].set_visible(False)
a.spines['bottom'].set_color(GRID); a.tick_params(length=0)

b=ax[1]
b.scatter(col,share,s=60,color=c,edgecolor=SURF,linewidth=1.4,zorder=3)
off={'SSC4':(7,-11),'SSC6':(-4,7),'SSC10':(6,-11),'SSC8':(-30,2),'SSC3':(6,4),
     'SSC9':(6,-10),'SSC1':(6,-11),'HC3':(6,3)}
for x0,y0,n in zip(col,share,names):
    dx,dy=off.get(n,(7,3))
    b.annotate(n,(x0,y0),textcoords='offset points',xytext=(dx,dy),fontsize=7.5,color=MUT)
b.set_xlabel('COL1A1 detected, fraction of QC nuclei',fontsize=9.5,color=MUT)
b.set_ylabel('fibroblast share of QC nuclei',fontsize=9.5,color=MUT)
b.set_title('The transcript agrees with the label',fontsize=11.5,color=INK,loc='left',pad=10)
b.set_xlim(0,1.06); b.set_ylim(-.02,.85)
b.set_xticks([0,.25,.5,.75,1]); b.set_xticklabels(['0','25%','50%','75%','100%'],fontsize=9,color=MUT)
b.set_yticks([0,.2,.4,.6,.8]); b.set_yticklabels(['0','20%','40%','60%','80%'],fontsize=9,color=MUT)
b.grid(True,color=GRID,lw=.8); b.set_axisbelow(True)
for s in ('top','right'): b.spines[s].set_visible(False)
for s in ('bottom','left'): b.spines[s].set_color(GRID)
b.tick_params(length=0)

h=[plt.Line2D([],[],marker='o',ls='',color=SSC,ms=7,label='dcSSc (n=10)'),
   plt.Line2D([],[],marker='o',ls='',color=HC,ms=7,label='healthy control (n=4)')]
b.legend(handles=h,loc='lower right',frameon=False,fontsize=9,labelcolor=INK)
fig.text(.006,.014,'Median fibroblast share 6.1% dcSSc vs 38.3% HC, exact Mann-Whitney p = 0.014. '
 'SSC7 sits high on share only because it yielded 267 nuclei in total; its absolute count is 89.',
 fontsize=8,color=MUT)
fig.tight_layout(rect=[0,.045,1,1])
fig.savefig('composition_fig.png',dpi=220,facecolor=SURF)
fig.savefig('composition_fig.pdf',facecolor=SURF)
print('ok')
