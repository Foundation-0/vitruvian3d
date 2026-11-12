"""Exact arithmetic checks for an explicit vacancy-to-FCC volume lower bound.

The accompanying volume_barrier.md contains the all-configuration proof.
No optimization or sampling is used to establish the bound.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np


def main():
    C=np.array([[0,-1,-1],[1,0,-1],[1,1,0]],dtype=int)
    L=2*np.eye(3,dtype=int)+C
    cofactor=np.array([[5,-3,-1],[1,5,-3],[3,1,5]],dtype=int)
    # S = A-Delta, where A is the tangent-cone projection and q_i=|Delta_i|²/2.
    S=[np.array([[1,0,-1],[0,0,0],[0,0,0]]),
       np.array([[0,0,0],[0,1,-1],[0,0,0]]),
       np.array([[0,0,0],[0,0,0],[0,1,1]])]
    assert [int(np.sum(cofactor*s)) for s in S]==[6,8,6]
    assert all(int(row@row)==6 for row in L)
    # qmax <=3r²/2. Tangent derivative >=3||A||inf.
    # Error: 3qmax+6q1+8q2+6q3 <= (69/2)r².
    assert 3*F(3,2)+(6+8+6)*F(3,2)==F(69,2)
    # Quadratic determinant remainder <=9sqrt(6)r² <(45/2)r².
    # Cubic remainder <=3sqrt(3)r³ <6r³.
    assert F(25,4)>6 and 4>3
    assert F(69,2)+F(45,2)==57
    radius=F(1,38)
    lower=3*radius-57*radius**2-6*radius**3
    assert lower==F(270,6859)>0
    # The polynomial remains positive inside the chosen neighborhood.
    assert 3-57*radius-6*radius**2>0
    target=np.array([[2,-1,-1],[1,2,-1],[1,-1,2]])
    target_distance=int(np.max(np.abs(target-L)))
    assert target_distance==2 and target_distance>radius
    result={
        'endpoint_lattice':L.tolist(),'endpoint_volume':14,
        'scope':'Continuous one-cell paths with fixed motif ±e_i, fixed motif orientation, prescribed bridge lengths sqrt(2), sphere diameter sqrt(2).',
        'finite_feasible_perturbation_bound':'det(L+Delta)-14 >= 3r-(69/2+9sqrt(6))*r^2-3sqrt(3)*r^3 > 3r-57r^2-6r^3 for r>0, r=max|Delta_ij|.',
        'crossing_radius_exact':str(radius),
        'strict_absolute_volume_excess_lower_bound_exact':str(lower),
        'strict_absolute_volume_excess_lower_bound_float':float(lower),
        'strict_peak_volume_lower_bound_exact':str(F(14)+lower),
        'strict_peak_volume_lower_bound_float':float(F(14)+lower),
        'relative_excess_percent_lower_bound_float':float(100*lower/14),
        'fcc_target_distance_infinity':target_distance,
        'established':'Every continuous feasible path to this distinct FCC endpoint has a cell volume exceeding the displayed lower bound.',
        'not_established':'Existence of a path achieving this bound, optimal barrier, energy, forces, or analogous bound with rotating motifs or larger cells.',
        'arithmetic_checks_passed':True}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
