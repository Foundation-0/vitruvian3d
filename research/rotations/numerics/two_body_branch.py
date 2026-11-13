"""Finite two-body branch with a fixed exterior: numerical companion evidence.

The exact IFT argument is separate. Every finite state and interval below is
numerical, not a certified continuation. No remote or existing file is changed.
"""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import root,linprog,least_squares
from scipy.spatial.transform import Rotation

HERE=Path(__file__).parent
P=np.r_[np.eye(3),-np.eye(3)]
C=np.array([[0,-1,-1],[1,0,-1],[1,1,0]],dtype=int)
PATCH=[(0,0,0),(1,0,0)]
AXES=np.eye(3,dtype=int)


def model(theta,patch=PATCH):
    angle=math.radians(theta)
    L=(2+math.sqrt(2)*math.cos(angle))*np.eye(3)+math.sin(angle)*C
    index={n:i for i,n in enumerate(patch)}
    starts=set()
    for n in patch:
        for i in range(3):
            starts.add((n,i));starts.add((tuple(np.array(n)-AXES[i]),i))
    bars=[]
    for n,i in sorted(starts):
        m=tuple(np.array(n)+AXES[i]);bars.append(((n,i),(m,i+3)))
    return L,index,bars


def transforms(x,patch=PATCH):
    x=np.array(x).reshape(len(patch),6)
    return x[:,:3],Rotation.from_rotvec(x[:,3:]).as_matrix()


def position(site,x,L,index,patch=PATCH):
    n,p=site
    if n not in index:return np.array(n)@L+P[p]
    j=index[n];u,R=transforms(x,patch)
    return np.array(n)@L+u[j]+R[j]@P[p]


def residual(x,L,index,bars):
    return np.array([np.dot(d:=position(a,x,L,index)-position(b,x,L,index),d)-2 for a,b in bars])


def initial_jacobian(L,index,bars,patch=PATCH):
    rows=[]
    for a,b in bars:
        d=np.array(a[0])@L+P[a[1]]-np.array(b[0])@L-P[b[1]]
        row=np.zeros(6*len(patch))
        for site,sgn in [(a,1),(b,-1)]:
            n,p=site
            if n in index:
                j=index[n];row[6*j:6*j+3]+=sgn*d
                row[6*j+3:6*j+6]+=sgn*np.cross(P[p],d)
        rows.append(row)
    return np.array(rows)


def pair_key(a,b):return tuple(sorted((a,b)))


def all_contact_rows(L,index,patch=PATCH):
    sites=set();rows=[];pairs=[]
    for n in patch:
        for p in range(6):
            a=(n,p)
            for shift in itertools.product(range(-1,2),repeat=3):
                m=tuple(np.array(n)+shift)
                for q in range(6):
                    b=(m,q)
                    if a==b:continue
                    key=pair_key(a,b)
                    if key in sites:continue
                    d=np.array(n)@L+P[p]-np.array(m)@L-P[q]
                    if abs(d@d-2)>1e-9:continue
                    sites.add(key)
                    row=np.zeros(6*len(patch))
                    for site,sgn in [(a,1),(b,-1)]:
                        if site[0] in index:
                            j=index[site[0]];row[6*j:6*j+3]+=sgn*d
                            row[6*j+3:6*j+6]+=sgn*np.cross(P[site[1]],d)
                    if np.linalg.norm(row)>1e-10:rows.append(row);pairs.append(key)
    return np.array(rows),pairs


def cone_probe(theta,patch):
    L,index,bars=model(theta,patch);J=initial_jacobian(L,index,bars,patch)
    A,pairs=all_contact_rows(L,index,patch)
    rank=int(np.linalg.matrix_rank(J,tol=1e-9));dof=6*len(patch)
    result={'angle':theta,'patch':[list(n) for n in patch],'motifs':len(patch),'bars':len(bars),'bar_rank':rank,'bar_nullity':dof-rank,'all_contact_linear_rank':int(np.linalg.matrix_rank(np.r_[J,A],tol=1e-9))}
    # Seek an exactly checkable strictly positive packing-contact stress.
    eq=np.c_[A.T,J.T]
    sol=linprog(np.r_[np.ones(len(A)),np.zeros(len(J))],A_eq=eq,b_eq=np.zeros(dof),bounds=[(1,None)]*len(A)+[(None,None)]*len(J),method='highs')
    result['positive_contact_stress_found']=bool(sol.success)
    if sol.success and theta==90:
        ints=np.rint(np.r_[A,J]).astype(int)
        assert np.max(np.abs(np.r_[A,J]-ints))<1e-12
        weights=[Fraction(float(x)).limit_denominator(100000) for x in sol.x]
        exact=[sum(Fraction(int(ints[i,j]))*weights[i] for i in range(len(weights))) for j in range(dof)]
        result['stress_exact_rational_equilibrium']=all(x==0 for x in exact)
        result['stress_contact_minimum']=str(min(weights[:len(A)]))
        result['contact_weights']=[str(x) for x in weights[:len(A)]]
        result['bar_multipliers']=[str(x) for x in weights[len(A):]]
        result['contact_rows_integer']=ints[:len(A)].tolist();result['bar_rows_integer']=ints[len(A):].tolist()
        result['contact_pairs']=pairs;result['bar_pairs']=bars
    return result


def clearance_data(x,L,index,bars):
    barset={pair_key(a,b) for a,b in bars}
    fixed=[]
    for n in itertools.product(range(-2,4),repeat=3):
        if n not in index:
            for p in range(6):fixed.append(((n,p),np.array(n)@L+P[p]))
    moving=[((n,p),position((n,p),x,L,index)) for n in PATCH for p in range(6)]
    best=(math.inf,None)
    for a,xa in moving:
        for b,xb in fixed:
            if pair_key(a,b) in barset:continue
            gap=float(np.dot(xa-xb,xa-xb)-2)
            if gap<best[0]:best=(gap,[a,b])
    for i,(a,xa) in enumerate(moving):
        for b,xb in moving[i+1:]:
            if a[0]==b[0] or pair_key(a,b) in barset:continue
            gap=float(np.dot(xa-xb,xa-xb)-2)
            if gap<best[0]:best=(gap,[a,b])
    u,_=transforms(x)
    maxsite=max(np.linalg.norm(y) for _,y in moving)
    # Outside [-2,3]^3, min |n|>=3. Since sigma_min(L)>=3 and |p|=1,
    # a global norm bound gives |nL+p-x|>=8-maxsite. This is conservative.
    omitted_bound=8-maxsite
    return {'minimum_nonbar_squared_gap':best[0],'minimum_pair':best[1],'maximum_motif_translation':float(np.linalg.norm(u,axis=1).max()),'maximum_rotation_radians':float(np.linalg.norm(np.array(x).reshape(2,6)[:,3:],axis=1).max()),'omitted_image_distance_lower_bound':omitted_bound,'maximum_site_norm':maxsite}


def main():
    L,index,bars=model(45);J=initial_jacobian(L,index,bars)
    U,S,V=np.linalg.svd(J);tangent=V[-1]
    if tangent[3]<0:tangent=-tangent
    rho=math.sqrt(2)
    u=np.array([2*(rho*rho+1),-2*(rho+1),-2*(rho-1)])/(rho*(rho*rho+3))
    analytic=np.r_[u,[1,1,-1],u,[-1,-1,1]];analytic/=np.linalg.norm(analytic)
    if analytic@tangent<0:analytic=-analytic
    assert np.max(np.abs(tangent-analytic))<1e-12
    def equations(x,t):return np.r_[residual(x,L,index,bars),tangent@x-t]
    states=[];terminations=[]
    zero=np.zeros(12)
    for sign in [-1,1]:
        previous=zero.copy();previous_t=0
        for k in range(41):
            t=sign*k*.005
            guess=previous+(t-previous_t)*tangent
            sol=root(lambda x:equations(x,t),guess,method='hybr',tol=1e-11)
            if np.max(np.abs(equations(sol.x,t)))>1e-9:
                sol=least_squares(lambda x:equations(x,t),guess,ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=300)
            if np.max(np.abs(equations(sol.x,t)))>1e-8:
                terminations.append({'parameter':t,'reason':'numerical solver did not close augmented equations','residual':float(np.max(np.abs(equations(sol.x,t))))});break
            x=sol.x;previous=x
            previous_t=t
            clear=clearance_data(x,L,index,bars)
            row={'parameter':t,'coordinates':x.tolist(),'maximum_bar_squared_length_residual':float(np.abs(residual(x,L,index,bars)).max()),**clear}
            states.append(row)
            if clear['minimum_nonbar_squared_gap']< -1e-8:break
    cone=[cone_probe(90,[(0,0,0)]),cone_probe(90,PATCH),cone_probe(45,[(0,0,0)]),cone_probe(45,PATCH)]
    unique={row['parameter']:row for row in states}
    result={'model':'Two adjacent rigid octahedra; all exterior motifs fixed;11incident contact bars retained. Initial skew open fold theta45. Finite branch states numerical, not a certified continuum.','lattice':L.tolist(),'bar_count':len(bars),'jacobian_singular_values':S.tolist(),'initial_jacobian_rank':int(np.linalg.matrix_rank(J)),'analytic_tangent_error':float(np.max(np.abs(tangent-analytic))),'normalized_tangent':tangent.tolist(),'states':[unique[t] for t in sorted(unique)],'state_count':len(unique),'terminations':terminations,'cone_controls':cone,'scope':'Original+negative branches sampled in .005 increments. Stored clearance is full local pair enumeration; omitted images have separate conservative norm bound. No claim of a finite certified motion until interval/implicit-function estimates are supplied.'}
    (HERE/'two_body_branch.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
    print(json.dumps({'rank':result['initial_jacobian_rank'],'tangent_error':result['analytic_tangent_error'],'parameter_range':[min(x['parameter'] for x in states),max(x['parameter'] for x in states)],'minimum_nonbar_squared_gap':min(x['minimum_nonbar_squared_gap'] for x in states),'maximum_bar_residual':max(x['maximum_bar_squared_length_residual'] for x in states),'maximum_translation':max(x['maximum_motif_translation'] for x in states),'maximum_rotation_radians':max(x['maximum_rotation_radians'] for x in states),'minimum_omitted_image_distance_bound':min(x['omitted_image_distance_lower_bound'] for x in states),'controls':[{k:v for k,v in x.items() if 'rows' not in k and 'weights' not in k and 'multipliers' not in k and 'pairs' not in k} for x in cone]},indent=2))


if __name__=='__main__':main()
