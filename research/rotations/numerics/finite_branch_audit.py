"""Independent floating-point audit using Rodrigues rotations, not scipy Rotation.
This checks stored frames, not the continuum between them.
"""
import itertools
import json
import math
import hashlib
import importlib.metadata
import platform
from pathlib import Path
import numpy as np

HERE=Path(__file__).parent
P=np.r_[np.eye(3),-np.eye(3)]
PATCH=[(0,0,0),(1,0,0)]


def rotation(w):
    angle=float(np.linalg.norm(w));K=np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    if angle<1e-10:return np.eye(3)+K+.5*K@K
    return np.eye(3)+(math.sin(angle)/angle)*K+((1-math.cos(angle))/angle**2)*K@K


def main():
    data=json.loads((HERE/'two_body_branch.json').read_text());L=np.array(data['lattice'])
    assert np.max(np.abs(L-(3*np.eye(3)+np.array([[0,-1,-1],[1,0,-1],[1,1,0]])/math.sqrt(2))))<1e-14
    index={n:i for i,n in enumerate(PATCH)}
    starts=set()
    for n in PATCH:
        for i in range(3):
            starts.add((n,i));back=list(n);back[i]-=1;starts.add((tuple(back),i))
    bars=[]
    for n,i in starts:
        m=list(n);m[i]+=1;bars.append(((n,i),(tuple(m),i+3)))
    assert len(bars)==11
    barset={tuple(sorted(e)) for e in bars}
    fixed=[]
    for n in itertools.product(range(-2,4),repeat=3):
        if n not in index:
            for p in range(6):fixed.append(((n,p),np.array(n)@L+P[p]))
    fixedpoints=np.array([p for _,p in fixed]);fixedlabels=[s for s,_ in fixed]
    def positions(x):
        x=np.array(x).reshape(2,6)
        return {(n,p):np.array(n)@L+x[j,:3]+rotation(x[j,3:])@P[p] for j,n in enumerate(PATCH) for p in range(6)}
    def constraints(x):
        moving=positions(x)
        def point(site):return moving[site] if site in moving else np.array(site[0])@L+P[site[1]]
        return np.array([np.dot(d:=point(a)-point(b),d)-2 for a,b in bars])
    worstres=0;mingap=math.inf;farlower=math.inf;minsing=math.inf;maxstep=0
    previous=None
    for state in data['states']:
        x=np.array(state['coordinates']);moving=positions(x)
        res=float(np.abs(constraints(x)).max());worstres=max(worstres,res);assert res<1e-10
        for a,point in moving.items():
            delta=fixedpoints-point;ds=np.einsum('ij,ij->i',delta,delta)-2
            for k,b in enumerate(fixedlabels):
                if tuple(sorted((a,b))) in barset:ds[k]=np.inf
            mingap=min(mingap,float(ds.min()))
        keys=list(moving)
        for i,a in enumerate(keys):
            for b in keys[i+1:]:
                if a[0]==b[0] or tuple(sorted((a,b))) in barset:continue
                d=moving[a]-moving[b];mingap=min(mingap,float(d@d-2))
        farlower=min(farlower,8-max(float(np.linalg.norm(p)) for p in moving.values()))
        # Independently approximate the nonlinear augmented Jacobian to detect sampled singularities.
        h=1e-6;J=np.empty((11,12))
        for i in range(12):
            step=np.zeros(12);step[i]=h
            J[:,i]=(constraints(x+step)-constraints(x-step))/(2*h)
        aug=np.r_[J,np.array(data['normalized_tangent'])[None,:]]
        minsing=min(minsing,float(np.linalg.svd(aug,compute_uv=False)[-1]))
        if previous is not None:maxstep=max(maxstep,float(np.linalg.norm(x-previous)))
        previous=x
    assert mingap>1.9 and farlower>math.sqrt(2) and minsing>.1
    unique_parameters=len({s['parameter'] for s in data['states']})
    assert unique_parameters==len(data['states'])==81
    result={'stored_frames':len(data['states']),'unique_parameters':unique_parameters,'parameter_range':[data['states'][0]['parameter'],data['states'][-1]['parameter']],
        'independent_rodrigues_maximum_bar_squared_residual':worstres,'independent_minimum_nonbar_squared_gap':mingap,
        'omitted_image_distance_lower_bound':farlower,'minimum_sampled_augmented_jacobian_singular_value':minsing,
        'maximum_adjacent_coordinate_step':maxstep,
        'source_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['two_body_branch.py','two_body_branch.json','finite_branch_audit.py']},
        'library_versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':importlib.metadata.version('scipy')},
        'scope':'Independent floating-point frame audit. Central differences and sampled singular values are diagnostics, not an interval continuation proof. The exact IFT result separately establishes a local continuum near the initial frame.'}
    (HERE/'finite_branch_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
