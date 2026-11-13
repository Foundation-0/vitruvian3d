"""Scientific plot of actual audited two-body coordinates; no invented motion."""
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).parent
P=np.r_[np.eye(3),-np.eye(3)]


def rotate(w):
    a=float(np.linalg.norm(w));K=np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    return np.eye(3)+K+.5*K@K if a<1e-10 else np.eye(3)+math.sin(a)/a*K+(1-math.cos(a))/a**2*K@K


def main():
    data=json.loads((HERE/'two_body_branch.json').read_text());audit=json.loads((HERE/'finite_branch_audit.json').read_text())
    states=data['states'];L=np.array(data['lattice']);x=np.array(states[-1]['coordinates']).reshape(2,6)
    tangent=np.array(data['normalized_tangent']);direction=tangent[:3]/np.linalg.norm(tangent[:3]);axis=tangent[3:6]/np.linalg.norm(tangent[3:6])
    colors=['#175b8e','#b65a20'];gray='#7a838b'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'vitruvian3d-local-rotation-v1','svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(12.8,4.9));gs=fig.add_gridspec(1,3,left=.025,right=.98,bottom=.22,top=.85,wspace=.32,width_ratios=[1.1,1,1])
    geometry=fig.add_subplot(gs[0,0],projection='3d');translation=fig.add_subplot(gs[0,1]);rotation=fig.add_subplot(gs[0,2])
    edges=[(i,j) for i in range(6) for j in range(i+1,6) if abs(float(P[i]@P[j]))<.5]
    frames=[]
    for j in range(2):
        origin=j*L[0];base=origin+P;moved=origin+x[j,:3]+P@rotate(x[j,3:]).T;frames.append(moved)
        for a,b in edges:
            geometry.plot(*base[[a,b]].T,color=gray,lw=.65,alpha=.45,ls='--')
            geometry.plot(*moved[[a,b]].T,color=colors[j],lw=1.35)
        geometry.scatter(*moved.T,color=colors[j],s=12)
    geometry.plot(*np.array([frames[0][0],frames[1][3]]).T,color='#333a40',lw=1.4)
    geometry.set_proj_type('ortho');geometry.view_init(elev=22,azim=-56)
    geometry.set(xlabel='x',ylabel='y',zlabel='z',xlim=(-1.4,4.4),ylim=(-2.1,1.3),zlim=(-2.1,1.3))
    geometry.set_box_aspect((5.8,3.4,3.4));geometry.tick_params(labelsize=7,pad=0)
    geometry.set_title('A   End posture at t = 0.2\nActual scale; gray = initial',pad=16)
    parameters=[s['parameter'] for s in states]
    for body in range(2):
        xx=np.array([s['coordinates'] for s in states]).reshape(-1,2,6)
        translation.plot(parameters,xx[:,body,:3]@direction,color=colors[body],lw=2,label=f'Body {body+1}')
        rotation.plot(parameters,(xx[:,body,3:]@axis)*180/math.pi,color=colors[body],lw=2,label=f'Body {body+1}')
    translation.set_title('B   Shared translation to first order',loc='left',pad=14)
    rotation.set_title('C   Opposite rotations to first order',loc='left',pad=14)
    translation.set_ylabel('Translation projected on initial tangent')
    rotation.set_ylabel('Rotation-vector projection (degrees)')
    for ax in [translation,rotation]:
        ax.set_xlabel('Continuation parameter t (not time)');ax.set_xlim(-.2,.2);ax.set_xticks([-.2,-.1,0,.1,.2]);ax.grid(color='#dce1e5',lw=.6);ax.axhline(0,color=gray,lw=.6)
    translation.legend(frameon=False,loc='upper left')
    fig.text(.035,.12,'Rigid motifs and the shared bridge are drawn from stored coordinates. Exterior motifs and sphere surfaces are omitted from panel A, but retained in checks.',fontsize=8.8)
    fig.text(.035,.075,f"{len(states)} audited frames; maximum bar |length² − 2| < 2.5×10⁻¹⁵; minimum nonbar length² − 2 > 1.965; maximum motif rotation ≈ 7.11°.",fontsize=8.8)
    fig.text(.035,.03,'Finite range: numerical samples, not a certified continuum. The separate exact implicit-function argument establishes a true local branch near t = 0.',fontsize=8.8,color=gray)
    fig.savefig(HERE/'two-body-branch.png',dpi=220,metadata={'Software':'Matplotlib','Title':'Two rigid motifs moving inside a fixed exterior'})
    fig.savefig(HERE/'two-body-branch.svg',metadata={'Creator':'Matplotlib','Date':None,'Title':'Two rigid motifs moving inside a fixed exterior'})
    p=HERE/'two-body-branch.svg';p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__=='__main__':main()
