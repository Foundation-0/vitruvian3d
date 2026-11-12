"""Numerical search for lower-volume switching routes; samples are not certificates."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

P=np.vstack([np.eye(3),-np.eye(3)])
V=np.unique((P[:,None]-P[None,:]).reshape(-1,3),axis=0)
ns=[];vs=[]
for n in itertools.product([-1,0,1],repeat=3):
    if n==(0,0,0):continue
    if next(x for x in n if x)!=1:continue
    for v in V:
        if sum(x*x for x in n)==1 and np.array_equal(v,-2*np.array(n)):continue
        ns.append(n);vs.append(v)
N=np.array(ns);VEC=np.array(vs)


def lattice(theta,beta):
    a=np.sqrt(2)*np.cos(theta);t=np.sin(theta)
    return np.array([[2+a[0],-t[0],-t[0]],
                     [t[1],2+a[1],-t[1]],
                     [np.sqrt(2)*t[2]*np.cos(beta),np.sqrt(2)*t[2]*np.sin(beta),2+a[2]]])


def gaps(theta,beta):
    d=N@lattice(theta,beta)+VEC
    return np.sum(d*d,axis=1)-2


def main():
    rows=[];last=np.ones(3)*np.pi/2
    for beta_deg in np.linspace(45,-45,37):
        beta=np.deg2rad(beta_deg)
        best=None
        for initial in [last,np.ones(3)*np.deg2rad(65),np.ones(3)*np.deg2rad(40)]:
            result=minimize(lambda t:np.linalg.det(lattice(t,beta)),initial,method='SLSQP',
                            bounds=[(0,np.pi/2)]*3,
                            constraints=[{'type':'ineq','fun':lambda t:gaps(t,beta)}],
                            options={'maxiter':250,'ftol':1e-11})
            if gaps(result.x,beta).min()>-1e-7 and (best is None or result.fun<best.fun):best=result
        assert best is not None
        last=best.x
        L=lattice(last,beta)
        row={'beta_deg':float(beta_deg),'theta_deg':np.rad2deg(last).tolist(),
             'volume':float(np.linalg.det(L)),'minimum_squared_gap':float(gaps(last,beta).min()),
             'sigma_min':float(np.linalg.svd(L,compute_uv=False)[-1]),'solver_success':bool(best.success)}
        rows.append(row)
        print(json.dumps(row),flush=True)
    result={'scope':'Numerical optimized sampled configurations in a four-parameter restricted family; no path or full periodic non-overlap certificate.',
            'rows':rows,'peak_sample_volume':max(r['volume'] for r in rows)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
