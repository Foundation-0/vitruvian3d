"""Stdlib-only replay of the vacancy hard-contact tangent-cone certificates."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).parent
L=((2,-1,-1),(1,2,-1),(1,1,2))
P=((1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1))


def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def point(site):
    n,p=site
    assert len(n)==3 and all(type(x) is int for x in n)
    assert type(p) is int and 0<=p<6
    return tuple(sum(n[i]*L[i][j] for i in range(3))+P[p][j] for j in range(3))


def pair_key(pair):return tuple(sorted((tuple(site[0]),site[1]) for site in pair))


def prescribed_bars(patch):
    expected=set()
    for n in patch:
        for i in range(3):
            forward=list(n);forward[i]+=1
            back=list(n);back[i]-=1
            expected.add(pair_key(((n,i),(forward,i+3))))
            expected.add(pair_key(((back,i),(n,i+3))))
    return expected


def row(pair,index):
    a,b=pair;x,y=point(a),point(b);d=tuple(s-t for s,t in zip(x,y))
    assert sum(t*t for t in d)==2
    answer=[0]*(6*len(index))
    for (n,p),sign in [(a,1),(b,-1)]:
        if tuple(n) in index:
            j=index[tuple(n)];coeff=(*d,*cross(P[p],d))
            for k in range(6):answer[6*j+k]+=sign*coeff[k]
    return answer


def rank(A):
    A=[[F(x) for x in r] for r in A];n=len(A[0]);pivot=0
    for j in range(n):
        i=next((i for i in range(pivot,len(A)) if A[i][j]),None)
        if i is None:continue
        A[pivot],A[i]=A[i],A[pivot];div=A[pivot][j]
        A[pivot]=[x/div for x in A[pivot]]
        for k in range(pivot+1,len(A)):
            if A[k][j]:
                q=A[k][j];A[k]=[x-q*y for x,y in zip(A[k],A[pivot])]
        pivot+=1
        if pivot==n:break
    return pivot


def main():
    data=json.loads((HERE/'two_body_branch.json').read_text());results=[];seen=set()
    allowed={((0,0,0),),((0,0,0),(1,0,0))}
    for c in data['cone_controls']:
        if c['angle']!=90:continue
        index={tuple(n):i for i,n in enumerate(c['patch'])}
        patch=tuple(tuple(n) for n in c['patch'])
        assert patch in allowed and patch not in seen and c['motifs']==len(patch)
        seen.add(patch)
        expected=prescribed_bars(patch)
        assert {pair_key(pair) for pair in c['bar_pairs']}==expected
        assert len(c['bar_pairs'])==len(expected)==(6 if len(patch)==1 else 11)
        A=[row(pair,index) for pair in c['contact_pairs']]
        J=[row(pair,index) for pair in c['bar_pairs']]
        assert A==c['contact_rows_integer'] and J==c['bar_rows_integer']
        w=list(map(F,c['contact_weights']));lam=list(map(F,c['bar_multipliers']))
        assert len(w)==len(A) and len(lam)==len(J) and min(w)>0
        dof=6*len(index)
        assert all(sum(w[i]*A[i][j] for i in range(len(A)))+sum(lam[i]*J[i][j] for i in range(len(J)))==0 for j in range(dof))
        assert rank(A+J)==dof
        results.append({'motifs':len(index),'degrees_of_freedom':dof,'bar_rank_exact':rank(J),'full_contact_rank_exact':dof,'positive_stress_equilibrium_exact':True,'sphere_contact_rows':len(A),'minimum_positive_weight':str(min(w)),'conclusion':'All retained-bar-compatible hard-contact first-order velocities are zero.'})
    assert seen==allowed and len(results)==2
    result={'verified':results,'proof':'For Jv=0 and Av>=0, positive weights w and equilibrium w^T A+lambda^T J=0 imply every Av row is zero. Full stacked rank then forces v=0.','scope':'Vacancy state, fixed exterior, rigid-body variables for the specified one- and two-motif patches. No statement about every finite patch follows.'}
    (HERE/'vacancy_stress_verify.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
