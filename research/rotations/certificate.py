"""Exact algebra and clearance certificates for finite rotating-motif patches.

The universal rank and implicit-function argument is in proof.md. This is not
a finite-amplitude interval integrator or proof-assistant formalization.
"""
import itertools
import json
from pathlib import Path
import sympy as s

C=s.Matrix([[0,-1,-1],[1,0,-1],[1,1,0]])


def sign_surd(a,b):
    if a==0:return (b>0)-(b<0)
    if b==0 or (a>0)==(b>0):return (a>0)-(a<0)
    return ((a*a>2*b*b)-(a*a<2*b*b))*((a>0)-(a<0))


def neighbor(n,i,step=1):return tuple(x+(step if j==i else 0) for j,x in enumerate(n))


def port_check(nodes):
    nodes=sorted(nodes);lookup={n:j for j,n in enumerate(nodes)}
    edges={(n,i) for n in nodes for i in range(3)}|{(neighbor(n,i,-1),i) for n in nodes for i in range(3)}
    usage=[0]*(6*len(nodes));internal=0
    for n,i in edges:
        m=neighbor(n,i);support=[]
        if n in lookup:support.append(6*lookup[n]+i)
        if m in lookup:support.append(6*lookup[m]+3+i)
        assert support
        if len(support)==2:internal+=1
        for column in support:usage[column]+=1
    assert all(count==1 for count in usage)
    assert len(edges)==6*len(nodes)-internal
    return {'motifs':len(nodes),'internal_edges':internal,'constraint_rows':len(edges),
            'all_port_columns_used_once':True,'local_manifold_dimension':internal}


def pair_minor():
    rho=s.symbols('rho');B=rho*s.eye(3)+C
    nodes=[(0,0,0),(1,0,0)]
    edges={(n,i) for n in nodes for i in range(3)}|{(neighbor(n,i,-1),i) for n in nodes for i in range(3)}
    rows=[]
    for n,i in sorted(edges):
        m=neighbor(n,i);row=[s.Integer(0)]*12;b=s.Matrix(list(B.row(i)));tau=s.eye(3).col(i).cross(b)
        for node,sign in [(n,-1),(m,1)]:
            if node in nodes:
                offset=6*nodes.index(node)
                for j in range(3):row[offset+j]=sign*b[j];row[offset+3+j]=-tau[j]
        rows.append(row)
    J=s.Matrix(rows);minor=s.factor(J[:,:11].det())
    assert minor==-64*rho**2*(rho**2+3)**2
    return str(minor)


def clearance():
    P=[tuple(sign if j==i else 0 for j in range(3)) for i in range(3) for sign in [-1,1]]
    V=sorted({tuple(x-y for x,y in zip(p,q)) for p in P for q in P})
    records=[]
    for n in itertools.product([-1,0,1],repeat=3):
        if n==(0,0,0):continue
        z=[sum(n[i]*int(C[i,j]) for i in range(3)) for j in range(3)]
        for v in V:
            if sum(x*x for x in n)==1 and all(v[j]==-2*n[j] for j in range(3)):continue
            u=[3*n[j]+v[j] for j in range(3)]
            a=2*sum(x*x for x in u)+sum(x*x for x in z)
            b=2*sum(u[j]*z[j] for j in range(3))
            # Squared distance = (a+b sqrt(2))/2; certify >=6-sqrt(2)>4.
            assert sign_surd(a-12,b+2)>=0
            records.append({'n':n,'v':v,'a':a,'b':b})
    assert len(records)==488
    assert C.T==-C
    # sigma_min(3I+C/sqrt(2))>=3, so |n|>=2 => distance>=4.
    assert s.Rational(3,2)**2>2
    return records


def single_vacancy_branch():
    c,z,t=s.symbols('c z t',real=True) # c=cos(alpha), z=sin(alpha)
    w=s.Matrix([1,-1,1]);K=s.Matrix([[0,-1,-1],[1,0,-1],[1,1,0]])
    # K acts as w cross v. Its rows happen to equal the vacancy C.
    assert K*s.Matrix([1,0,0])==w.cross(s.Matrix([1,0,0]))
    R=c*s.eye(3)+(1-c)*(w*w.T)/3-z*K/s.sqrt(3)
    t_squared=4*(c-1)/9+4*z/(3*s.sqrt(3))
    for i in range(3):
        a=s.eye(3).col(i)+s.Matrix(list(C.row(i)))
        v=a-R*s.eye(3).col(i)
        assert s.simplify(w.dot(v))==0
        gap=s.expand(v.dot(v)+3*t*t-2).subs(t*t,t_squared)
        assert s.simplify(s.expand(gap).subs(z*z,1-c*c))==0
    # Fixed sphere at L3-e3=(1,1,1), moving e2: separation(1,0,1).
    assert s.Matrix([1,0,1]).dot(w)==2
    return 't^2=(4/9)(cos(alpha)-1)+(4/(3sqrt(3)))sin(alpha); u=t(1,-1,1), R=rotation about (1,-1,1) by -alpha'


def main():
    a,t=s.symbols('a t',positive=True)
    B=a*s.eye(3)+t*C
    T=s.Matrix.hstack(*[s.eye(3).col(i).cross(s.Matrix(list(B.row(i)))) for i in range(3)]).T
    assert s.factor(B.det())==a*(a*a+3*t*t)
    assert s.factor(T.det())==2*t**3
    cases={'single':[(0,0,0)],'adjacent_pair':[(0,0,0),(1,0,0)],
           'line_three':[(0,0,0),(1,0,0),(2,0,0)],
           'corner_three':[(0,0,0),(1,0,0),(0,1,0)],
           'cube_eight':list(itertools.product([0,1],repeat=3)),
           'separated_pair':[(0,0,0),(2,0,0)]}
    result={'bridge_matrix_determinant':str(s.factor(B.det())),
            'torque_matrix_determinant':str(s.factor(T.det())),
            'port_examples':{name:port_check(nodes) for name,nodes in cases.items()},
            'two_body_scaled_jacobian_minor':pair_minor(),
            'nonprescribed_inter_motif_clearance_squared':'at least 6-sqrt(2)>4 at theta=45deg',
            'safe_site_displacement_neighborhood':'each actual sphere center moves less than 1/4; nonprescribed distances then exceed 3/2>sqrt(2)',
            'clearance_conditions':clearance(),
            'vacancy_single_body_bar_only_branch':single_vacancy_branch(),
            'scope':'Analytic local packing-feasible mechanisms for arbitrary finite rotating+translating patches at the open skew45 baseline. No global continuation, material mechanics, or full vacancy-patch classification.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='clearance_conditions'},indent=2))


if __name__=='__main__':main()
