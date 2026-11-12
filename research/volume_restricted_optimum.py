"""Exact algebra checks for minimax optimality with the first two rows fixed.

The topological retraction argument is in volume_restricted_optimum.md.
SymPy is already installed. No numerical optimizer is used for this result.
"""
import itertools
import json
from pathlib import Path
import sympy as s


def main():
    # c+a dot b >=0; final row z>=0 follows from single-row sphere contacts.
    constraints=[(1,-1,0,1),(1,0,-1,1),(1,1,-1,0),(1,1,1,0),(1,0,1,0),(0,0,0,1)]
    valid=set()
    for row1,row2 in itertools.combinations(constraints,2):
        A=s.Matrix([row1[1:],row2[1:]])
        if A.rank()!=2:continue
        rhs=-s.Matrix([row1[0],row2[0]])
        center=A.T*(A*A.T).inv()*rhs
        axis=s.Matrix(row1[1:]).cross(s.Matrix(row2[1:]))
        scale=s.simplify((2-center.dot(center))/axis.dot(axis))
        if scale<0:continue
        for sign in [-1,1]:
            point=center+sign*s.sqrt(scale)*axis
            gaps=[s.simplify(r[0]+s.Matrix(r[1:]).dot(point)) for r in constraints]
            if all(bool(g>=0) for g in gaps):valid.add(tuple(s.simplify(x) for x in point))
    expected={(s.Integer(1),s.Integer(1),s.Integer(0)),(s.Integer(1),s.Integer(-1),s.Integer(0)),
              (s.Rational(1,3),s.Rational(4,3),s.Rational(1,3)),
              (s.Integer(0),s.Integer(-1),s.Integer(1)),(s.Integer(-1),s.Integer(0),s.Integer(1))}
    assert valid==expected
    # Verify each necessary contact by direct squared-gap reduction.
    p,q,z=s.symbols('p q z',real=True)
    L=s.Matrix([[2,-1,-1],[1,2,-1],[p,q,2+z]])
    witnesses=[((-0,0,-1),(1,0,1),(1,-1,0,1),2),
               ((-1,0,-1),(2,0,0),(1,0,-1,1),2),
               ((-1,0,-1),(1,0,1),(1,1,-1,0),2),
               ((-1,-1,-1),(2,0,0),(1,1,1,0),2),
               ((0,-1,-1),(1,0,1),(1,0,1,0),4)]
    for n,v,row,factor in witnesses:
        d=s.Matrix(1,3,n)*L+s.Matrix(1,3,v)
        gap=s.expand(d.dot(d)-2)
        reduced=s.Poly(gap-(p*p+q*q+z*z-2),p,q,z).as_expr()
        assert s.expand(reduced-factor*(row[0]+row[1]*p+row[2]*q+row[3]*z))==0
    peak=9+s.Rational(3,2)*s.sqrt(22)
    edges={'AB':peak,'AC':8+s.Rational(9,2)*s.sqrt(2),'CE':peak,
           'ED':8+s.Rational(9,2)*s.sqrt(2),'DB':9+s.sqrt(34)}
    t=s.symbols('t',real=True)
    arc_checks={
        'AB':(9+t+4*s.sqrt(3-2*t*t),1/s.sqrt(22),-1,1),
        'AC':(8+3*t+3*s.sqrt(3-2*t*t),1/s.sqrt(2),s.Rational(1,3),1),
        'CE':(11+4*t+5*s.sqrt(1-2*t-2*t*t),-s.Rational(1,2)+2/s.sqrt(22),-1,s.Rational(1,3)),
        'ED':(9+2*t+5*s.sqrt(1-2*t-2*t*t),-s.Rational(1,2)+1/s.sqrt(18),-1,0),
        'DB':(9+3*t+5*s.sqrt(1-t*t),3/s.sqrt(34),0,1)}
    for label,(expression,critical,lower,upper) in arc_checks.items():
        assert bool(lower<critical) and bool(critical<upper)
        assert s.simplify(s.diff(expression,t).subs(t,critical))==0
        assert s.simplify(expression.subs(t,critical)-edges[label])==0
        assert bool(expression.subs(t,lower)>=12) and bool(expression.subs(t,upper)>=12)
    assert all(bool(value<=peak) for value in edges.values())
    h=s.Matrix([3,1,5])
    assert all(h.dot(s.Matrix(row[1:]))>0 for row in constraints)
    maximum_direction=s.sqrt(s.Rational(2,35))*h
    assert all(bool(s.simplify(row[0]+s.Matrix(row[1:]).dot(maximum_direction))>0) for row in constraints)
    result={'scope':'Exact minimax optimum with L1=(2,-1,-1), L2=(1,2,-1) fixed and only the third bridge moving; not unrestricted all-row optimum.',
            'necessary_spherical_region_constraints':[list(row) for row in constraints],
            'vertices':[[str(x) for x in point] for point in sorted(valid,key=lambda t:tuple(float(x) for x in t))],
            'boundary_edge_peak_volumes':{k:str(v) for k,v in edges.items()},
            'restricted_minimax_peak_exact':str(peak),'restricted_minimax_peak_float':float(peak),
            'arithmetic_checks_passed':True,
            'proof_boundary':'Algebra checks plus the explicit radial-retraction and boundary-path argument in volume_restricted_optimum.md; not formal proof-assistant output.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
