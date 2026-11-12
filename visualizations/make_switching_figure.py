"""Reproduce the exact switching-volume figure with Matplotlib.

Proof sources, pinned before plotting:
https://github.com/Foundation-0/vitruvian3d/blob/a38ecc78910e4dc8dd8ed843ac7207341c0151a9/research/volume_layer_proof.md
https://github.com/Foundation-0/vitruvian3d/blob/65afe245178d6d0da1b748a5d5f52be5a7d865f4/research/volume_barrier.md

q is the moving-bridge coordinate, not a probe radius. Sampling below only
draws the proven formulas; the source proofs establish feasibility and bounds.
"""
from fractions import Fraction
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter

OUT=Path(__file__).resolve().parent
E=1.5*math.sqrt(22)-5
B=Fraction(270,6859)
QSTAR=1/math.sqrt(22)
M=9+1.5*math.sqrt(22)


def volume(q):
    return 9+q+4*math.sqrt(3-2*q*q)


def crossover_check(N):
    """Exact sign of E/N-B, reduced to a rational squared comparison."""
    rhs=Fraction(5)+N*B
    delta=Fraction(99,2)-rhs*rhs
    return (delta>0)-(delta<0)


def main():
    assert volume(1)==14 and volume(-1)==12
    assert abs(volume(QSTAR)-M)<1e-12
    assert crossover_check(51)==1 and crossover_check(52)==-1
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
        'axes.titlesize':12,'axes.titleweight':'normal','axes.spines.top':False,
        'axes.spines.right':False,'axes.edgecolor':'#4b5563','axes.labelcolor':'#20252b',
        'text.color':'#20252b','xtick.color':'#4b5563','ytick.color':'#4b5563',
        'svg.fonttype':'none','svg.hashsalt':'vitruvian3d-switching-volume-v1','pdf.fonttype':42})
    blue='#175b8e';orange='#b65a20';neutral='#68737d'
    fig,(left,right)=plt.subplots(1,2,figsize=(12.2,5.2))
    fig.subplots_adjust(left=.07,right=.98,bottom=.18,top=.88,wspace=.28)

    q=[1-2*k/1200 for k in range(1201)]
    left.plot(q,[volume(x) for x in q],color=blue,lw=2.4)
    left.axhline(14,color=neutral,lw=1,ls=(0,(4,3)),zorder=0)
    left.scatter([1,-1,QSTAR],[14,12,M],s=[30,30,38],color=blue,zorder=4)
    left.set(xlim=(1.08,-1.08),ylim=(11.45,17.2),xlabel='Bridge coordinate q (vacancy → FCC)',ylabel='Volume per six-site motif')
    left.set_title('A   Exact single-layer switching arc',loc='left',pad=14)
    left.set_xticks([1,.5,0,-.5,-1]);left.set_yticks([12,13,14,15,16,17])
    left.grid(axis='y',color='#d9dfe4',lw=.6,zorder=0)
    left.text(.96,13.7,'Vacancy endpoint\nV(1) = 14',ha='left',va='top',fontsize=9)
    left.text(-.98,11.82,'FCC endpoint\nV(−1) = 12',ha='right',va='bottom',fontsize=9)
    left.annotate('M = 9 + 1.5√22 ≈ 16.0356',xy=(QSTAR,M),xytext=(.85,16.76),
        arrowprops={'arrowstyle':'-','color':neutral,'lw':.8},ha='left',fontsize=9)
    left.text(-.15,15.75,r'$q_*=1/\sqrt{22}$',fontsize=10)
    left.annotate('',xy=(-.65,M),xytext=(-.65,14),arrowprops={'arrowstyle':'<->','color':orange,'lw':1.2})
    left.text(-.70,15.05,'E = M − 14\n≈ 2.0356',ha='left',va='center',color=orange,fontsize=9)
    left.text(.98,.96,r'$V(q)=9+q+4\sqrt{3-2q^2}$',transform=left.transAxes,ha='right',va='top',fontsize=11)

    ns=list(range(1,257));percent=[100*E/(14*n) for n in ns]
    lower=100*float(B)/14
    right.plot(ns,percent,color=blue,lw=2.4,label='Sequential layer construction')
    right.axhline(lower,color=orange,lw=1.3,ls=(0,(5,3)))
    right.scatter([52],[percent[51]],s=42,color=blue,zorder=5)
    right.set(xscale='log',yscale='log',xlim=(1,280),ylim=(.045,20),xlabel='Number of independently switching layers, N',ylabel='Peak expansion relative to initial volume (%)')
    right.set_title('B   Bulk expansion decreases as 1/N',loc='left',pad=14)
    right.xaxis.set_major_locator(FixedLocator([1,4,16,52,256]))
    right.xaxis.set_major_formatter(FuncFormatter(lambda x,_:str(int(x))))
    right.yaxis.set_major_locator(FixedLocator([.05,.1,.3,1,3,10]))
    right.yaxis.set_major_formatter(FuncFormatter(lambda y,_:f'{y:g}'))
    right.grid(which='major',color='#d9dfe4',lw=.6,zorder=0)
    right.text(.04,.50,r'$100E/(14N)$',transform=right.transAxes,color=blue,fontsize=11)
    right.text(1.2,lower*.77,'One-cell lower bound: B/14\nB = 270/6859;  0.281174%',color=orange,fontsize=8.7,va='top')
    right.annotate('N = 52: 0.279619%',xy=(52,percent[51]),xytext=(26,.72),
        arrowprops={'arrowstyle':'-','color':neutral,'lw':.8},fontsize=9)
    inset=right.inset_axes([.56,.58,.40,.35])
    inset.plot(range(48,57),[100*E/(14*n) for n in range(48,57)],color=blue,lw=1.7)
    inset.axhline(lower,color=orange,lw=1,ls=(0,(4,2)))
    inset.scatter([51,52],[100*E/(14*51),100*E/(14*52)],color=blue,s=20,zorder=4)
    inset.set(xlim=(49.8,54.2),ylim=(.267,.294),xticks=[50,51,52,53,54],yticks=[.27,.28,.29])
    inset.tick_params(labelsize=7.8,length=3)
    inset.set_title('First bound crossing: N = 52',fontsize=8.5,pad=5)
    inset.grid(color='#e3e7eb',lw=.5)

    fig.text(.07,.06,'Specified contact-preserving layer family. The N = 52 comparison uses a conservative one-cell lower bound, not its exact optimum.',fontsize=9)
    fig.text(.07,.025,'Kinematic volume bounds; no energy, force, fixed-box transition, or unrestricted supercell optimum is claimed.',fontsize=9,color=neutral)
    common={'Title':'Switching volume and layer-count scaling','Creator':'Matplotlib'}
    fig.savefig(OUT/'switching-volume.png',dpi=220,metadata={'Title':common['Title'],'Software':'Matplotlib'})
    fig.savefig(OUT/'switching-volume.svg',metadata={**common,'Date':None})
    svg=OUT/'switching-volume.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(OUT/'switching-volume.pdf',metadata={**common,'CreationDate':None,'ModDate':None,'Author':None})
    plt.close(fig)
    metrics={'formula':'V(q)=9+q+4sqrt(3-2q^2)','q_definition':'moving-bridge coordinate; not probe radius','q_star':'1/sqrt(22)',
        'maximum_volume_exact':'9+(3/2)sqrt(22)','maximum_volume_float':M,
        'excess_E_exact':'(3/2)sqrt(22)-5','excess_E_float':E,'one_cell_lower_excess_B_exact':'270/6859',
        'one_cell_lower_fractional_expansion_percent':lower,'N52_fractional_expansion_percent':percent[51],
        'first_integer_N_with_E_over_N_below_B':52,'exact_crossover_controls':{'E_over_51_greater_than_B':True,'E_over_52_less_than_B':True},
        'matplotlib_version':matplotlib.__version__,
        'proof_sources':['https://github.com/Foundation-0/vitruvian3d/blob/a38ecc78910e4dc8dd8ed843ac7207341c0151a9/research/volume_layer_proof.md',
            'https://github.com/Foundation-0/vitruvian3d/blob/65afe245178d6d0da1b748a5d5f52be5a7d865f4/research/volume_barrier.md']}
    (OUT/'switching-volume-values.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
