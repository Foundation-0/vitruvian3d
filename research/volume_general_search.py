"""Exploratory five-variable sampled search; explicitly not a continuum certificate."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

P=np.vstack([np.eye(3),-np.eye(3)])
V=np.unique((P[:,None]-P[None,:]).reshape(-1,3),axis=0)
ns=[];vs=[]
for n in itertools.product(range(-2,3),repeat=3):
    if n==(0,0,0):continue
    if next(x for x in n if x)<0:continue
    for v in V:
        if sum(x*x for x in n)==1 and np.array_equal(v,-2*np.array(n)):continue
        ns.append(n);vs.append(v)
N=np.array(ns);D=np.array(vs)


def lattice(x,q):
    t1,p1,t2,p2,a=x
    r=np.sqrt(2-q*q)
    return np.array([[2+np.sqrt(2)*np.cos(t1),np.sqrt(2)*np.sin(t1)*np.cos(p1),np.sqrt(2)*np.sin(t1)*np.sin(p1)],
                     [np.sqrt(2)*np.sin(t2)*np.cos(p2),2+np.sqrt(2)*np.cos(t2),np.sqrt(2)*np.sin(t2)*np.sin(p2)],
                     [r*np.cos(a),q,2+r*np.sin(a)]])


def gaps(x,q):
    d=N@lattice(x,q)+D
    return np.sum(d*d,axis=1)-2


def main():
    rows=[];last=None
    bounds=[(0,np.pi/2),(-np.pi,-np.pi/2),(0,np.pi/2),(-np.pi/2,0),(0,np.pi/2)]
    for q in np.linspace(1,-1,41):
        z=(np.sqrt(3-2*q*q)-1)/2
        reference=np.array([np.pi/2,-3*np.pi/4,np.pi/2,-np.pi/4,np.arctan2(z,1+z)])
        best=None
        for initial in [reference] if last is None else [last,reference]:
            res=minimize(lambda x:np.linalg.det(lattice(x,q)),initial,method='SLSQP',bounds=bounds,
                         constraints=[{'type':'ineq','fun':lambda x:gaps(x,q)}],
                         options={'ftol':2e-10,'maxiter':250})
            if gaps(res.x,q).min()>-1e-6 and (best is None or res.fun<best.fun):best=res
        assert best is not None
        last=best.x;L=lattice(last,q)
        row={'q':float(q),'parameters':best.x.tolist(),'lattice':L.tolist(),'volume':float(np.linalg.det(L)),
             'minimum_squared_gap':float(gaps(best.x,q).min()),
             'sigma_min':float(np.linalg.svd(L,compute_uv=False)[-1]),'solver_success':bool(best.success)}
        rows.append(row);print(json.dumps(row),flush=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps({'scope':'Only sampled configurations, not verified interpolation or global optimization.',
        'peak_sample_volume':max(r['volume'] for r in rows),'rows':rows},indent=2)+'\n')


if __name__=='__main__':main()
