"""Exact rational replay of fixed-boundary retained-bridge identities.

This checks identities on nontrivial finite/periodic fields; the universal
telescoping proof is in boundary_proof.md. No physical energy is modeled.
"""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path

B=((0,-1,-1),(1,0,-1),(1,1,0))
ZERO=(F(0),)*3


def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def neighbor(n,i,step=1):return tuple(x+(step if j==i else 0) for j,x in enumerate(n))


def finite_check(field):
    edges={(n,i) for n in field for i in range(3)}
    edges|={(neighbor(n,i,-1),i) for n in field for i in range(3)}
    energy=F(0);excess=F(0);linear=[F(0)]*3;sumdiff=[list(ZERO) for _ in range(3)]
    edge_excess=[]
    for n,i in edges:
        delta=sub(field.get(neighbor(n,i),ZERO),field.get(n,ZERO))
        strain=dot(add(B[i],delta),add(B[i],delta))-dot(B[i],B[i])
        energy+=dot(delta,delta);excess+=strain;linear[i]+=dot(B[i],delta)
        for j in range(3):sumdiff[i][j]+=delta[j]
        edge_excess.append(strain)
    assert excess==energy and all(x==0 for x in linear)
    assert all(all(x==0 for x in row) for row in sumdiff)
    return {'support_sites':len(field),'affected_or_incident_edges':len(edges),
            'total_squared_length_excess':str(excess),'discrete_dirichlet_sum':str(energy),
            'telescoping_linear_terms':[str(x) for x in linear],
            'max_squared_length_excess':str(max(edge_excess))}


def periodic_check(shape):
    nodes=list(itertools.product(*[range(n) for n in shape]))
    field={n:tuple(F(((j+2)*n[0]+(j+3)*n[1]+(j+5)*n[2])%7-3,10) for j in range(3)) for n in nodes}
    energy=F(0);excess=F(0);linear=F(0)
    for n in nodes:
        for i in range(3):
            m=tuple((n[j]+(i==j))%shape[j] for j in range(3))
            delta=sub(field[m],field[n]);new=add(B[i],delta)
            energy+=dot(delta,delta);excess+=dot(new,new)-dot(B[i],B[i]);linear+=dot(B[i],delta)
    assert linear==0 and energy==excess>0
    return {'periods':shape,'total_squared_length_excess':str(excess),'discrete_dirichlet_sum':str(energy)}


def main():
    cases=[]
    for R in [0,1,2]:
        field={n:tuple(F(((j+2)*n[0]+(j+3)*n[1]+(j+5)*n[2]+j+1)%9-4,20) for j in range(3))
               for n in itertools.product(range(-R,R+1),repeat=3)}
        result=finite_check(field)
        maximum_square=max(dot(u,u) for u in field.values())
        lower=6*maximum_square/F(R+1)
        assert F(result['discrete_dirichlet_sum'])>=lower
        result.update({'support_box_radius':R,'max_displacement_squared':str(maximum_square),
                       'line_path_dirichlet_lower_bound':str(lower)})
        cases.append(result)
    k=(1,-1,1)
    assert all(dot(row,k)==0 for row in B)
    epsilon=F(1,10)
    field={(0,0,0):tuple(epsilon*x for x in k)}
    example=finite_check(field)
    assert F(example['discrete_dirichlet_sum'])==18*epsilon**2
    assert F(example['max_squared_length_excess'])==3*epsilon**2
    # Every affected bridge has zero first-order change because B_i dot k=0.
    example['every_first_order_bridge_change_zero']=True
    example['every_affected_exact_squared_length_increase']=str(3*epsilon**2)
    example['displacement_direction']=list(k)
    result={'theorem':'No nontrivial finitely supported motif translation preserves every bridge length from a homogeneous baseline; fixed-period periodic translations are constant.',
            'identity':'sum_edges (|b_i+D_i u|^2-|b_i|^2)=sum_edges |D_i u|^2',
            'scope':'Fixed motif orientations and anchor offsets, all retained lattice bridges, homogeneous baseline, fixed exterior or fixed periodic cell. No sphere non-overlap assumption is needed.',
            'finite_rational_checks':cases,'periodic_rational_check':periodic_check((2,3,2)),
            'localized_infinitesimal_counterexample':example,
            'evidence_boundary':'Exact replay examples support the analytic telescoping proof; examples alone are not the universal proof. Dirichlet sum is a mathematical quadratic form, not a specified elastic energy.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
