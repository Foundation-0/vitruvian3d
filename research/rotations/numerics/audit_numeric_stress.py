"""Independent exact vacancy-stress and numerical-frame audit.

No author implementation is imported. Pass a data-directory argument, or place
this alongside the published files. In the scratch tree it finds rotation_numeric.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import numpy as np
import sympy as s

HERE=Path(__file__).parent
P=np.vstack([np.eye(3,dtype=int),-np.eye(3,dtype=int)])
VAC=np.array([[2,-1,-1],[1,2,-1],[1,1,2]],dtype=int)


def canonical(site):
    n,p=site
    assert len(n)==3 and all(type(x) is int for x in n)
    assert type(p) is int and 0<=p<6
    return tuple(n),p


def edgekey(pair):return tuple(sorted(map(canonical,pair)))


def bridges(patch):
    result=set()
    for n in patch:
        for i in range(3):
            for source in [n,tuple(n[j]-(i==j) for j in range(3))]:
                target=tuple(source[j]+(i==j) for j in range(3))
                result.add(tuple(sorted(((source,i),(target,i+3)))))
    return result


def exact_rows(pairs,patch):
    rows=[]
    for pair in pairs:
        left,right=map(canonical,pair)
        x=VAC.T@s.Matrix(left[0])+s.Matrix(P[left[1]].tolist())
        y=VAC.T@s.Matrix(right[0])+s.Matrix(P[right[1]].tolist())
        d=x-y
        assert d.dot(d)==2
        row=[]
        for node in patch:
            # Rebuild each Jacobian column from an actual unit rigid-body velocity.
            for mode in range(6):
                velocity=s.zeros(3,1)
                for (n,p),sign in [(left,1),(right,-1)]:
                    if n!=node:continue
                    axis=s.eye(3).col(mode%3)
                    value=axis if mode<3 else axis.cross(s.Matrix(P[p].tolist()))
                    velocity+=sign*value
                row.append(d.dot(velocity))
        rows.append(row)
    return s.Matrix(rows)


def exact_controls(data):
    outputs=[]
    for control in data['cone_controls']:
        if control['angle']!=90:continue
        patch=[tuple(n) for n in control['patch']]
        assert len(set(patch))==len(patch)
        assert {edgekey(pair) for pair in control['bar_pairs']}==bridges(patch)
        assert len(control['bar_pairs'])==len(bridges(patch))
        A=exact_rows(control['contact_pairs'],patch)
        J=exact_rows(control['bar_pairs'],patch)
        assert A==s.Matrix(control['contact_rows_integer'])
        assert J==s.Matrix(control['bar_rows_integer'])
        w=s.Matrix([s.Rational(x) for x in control['contact_weights']])
        lam=s.Matrix([s.Rational(x) for x in control['bar_multipliers']])
        assert all(x>0 for x in w)
        assert A.T*w+J.T*lam==s.zeros(6*len(patch),1)
        assert A.col_join(J).rank()==6*len(patch)
        outputs.append({'motifs':len(patch),'genuine_prescribed_bar_count':J.rows,
                        'contact_rows':A.rows,'bar_rank':J.rank(),'stacked_rank':A.col_join(J).rank(),
                        'minimum_positive_weight':str(min(w)),
                        'exact_equilibrium_and_label_reconstruction_passed':True})
    assert sorted(x['motifs'] for x in outputs)==[1,2]
    return outputs


def quaternion_rotate(w,p):
    angle=float(np.linalg.norm(w))
    vector=np.zeros(3) if angle==0 else math.sin(angle/2)*w/angle
    scalar=math.cos(angle/2)
    return (scalar*scalar-vector@vector)*p+2*vector*(vector@p)+2*scalar*np.cross(vector,p)


def finite_frames(data):
    L=np.array(data['lattice']);patch=[(0,0,0),(1,0,0)];barset=bridges(patch)
    expected=3*np.eye(3)+(VAC-2*np.eye(3))/math.sqrt(2)
    assert np.max(np.abs(L-expected))<1e-14
    fixedlabels=[(n,p) for n in itertools.product(range(-2,4),repeat=3) if n not in patch for p in range(6)]
    fixed=np.array([np.array(n)@L+P[p] for n,p in fixedlabels])
    residual=0.;gap=math.inf;farlower=math.inf;phase_error=0.
    parameters=[row['parameter'] for row in data['states']]
    assert len(parameters)==len(set(parameters))==81 and parameters==sorted(parameters)
    for state in data['states']:
        x=np.array(state['coordinates']).reshape(2,6)
        moving={(n,p):np.array(n)@L+x[j,:3]+quaternion_rotate(x[j,3:],P[p]) for j,n in enumerate(patch) for p in range(6)}
        def point(site):return moving[site] if site in moving else np.array(site[0])@L+P[site[1]]
        for a,b in barset:
            d=point(a)-point(b);residual=max(residual,abs(d@d-2))
        for a,y in moving.items():
            delta=fixed-y;ds=np.einsum('ij,ij->i',delta,delta)-2
            for j,b in enumerate(fixedlabels):
                if tuple(sorted((a,b))) in barset:ds[j]=np.inf
            gap=min(gap,float(ds.min()))
        keys=list(moving)
        for a,b in itertools.combinations(keys,2):
            if a[0]==b[0] or tuple(sorted((a,b))) in barset:continue
            d=moving[a]-moving[b];gap=min(gap,float(d@d-2))
        farlower=min(farlower,8-max(float(np.linalg.norm(y)) for y in moving.values()))
        phase_error=max(phase_error,abs(np.array(data['normalized_tangent'])@x.ravel()-state['parameter']))
    assert residual<1e-10 and gap>1.9 and farlower>math.sqrt(2) and phase_error<1e-9
    return {'distinct_frames':81,'quaternion_maximum_bar_squared_residual':residual,
            'minimum_nonbar_squared_gap':gap,'omitted_image_lower_bound':farlower,
            'maximum_phase_constraint_residual':phase_error,
            'scope':'Floating-point stored-frame audit only; no interpolation/continuation certification.'}


def main():
    if len(sys.argv)>1:source=Path(sys.argv[1])
    else:source=HERE if (HERE/'two_body_branch.json').exists() else HERE.parent/'rotation_numeric'
    names=['two_body_branch.py','two_body_branch.json','vacancy_stress_verify.py','vacancy_stress_verify.json','finite_branch_audit.py','finite_branch_audit.json','README.md']
    hashes={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names}
    data=json.loads((source/'two_body_branch.json').read_text())
    result={'input_sha256':hashes,'author_code_imported':False,
            'exact_vacancy_controls':exact_controls(data),'finite_frame_audit':finite_frames(data),
            'inference':'Zero hard-contact tangent cone implies local isolation by normalized secants in a smooth rigid-pose chart. Only the specified vacancy one-/two-motif patches are covered.',
            'verdict':'PASS with stated local/numerical limits.'}
    (HERE/'audit_numeric_stress.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
