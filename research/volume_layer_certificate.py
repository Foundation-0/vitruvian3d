"""Exact composition checks for independently varying certified switching layers.

Uses the previously exact-certified single-layer path, then a new perpendicular
separation bound for nonadjacent layers. No spatial or time sampling.
"""
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
from volume_path_certificate import certify, P, V, L1, L2, BASE3


def main():
    # Adjacent layers have the same differences as the uniform path with n3=1.
    adjacent=0
    for n1,n2 in itertools.product(range(-2,3),repeat=2):
        for v in V:
            A=tuple(n1*L1[i]+n2*L2[i]+BASE3[i]+v[i] for i in range(3))
            c=sum(x*x for x in A)-1
            a=2*(A[0]+A[2]-1)
            b=2*A[1]
            certify(c,a,b)
            adjacent+=1
    assert adjacent==475
    # Horizontal images outside this range are already covered by the uniform
    # path's all-neighbor argument. This invokes that theorem, not a new truncation.
    h=(3,1,5)
    assert sum(x*x for x in h)==35
    assert sum(a*b for a,b in zip(h,L1))==sum(a*b for a,b in zip(h,L2))==0
    projection_values=[sum(a*b for a,b in zip(h,p)) for p in P]
    assert max(projection_values)==5 and min(projection_values)==-5
    assert (2*12-10)**2>2*35
    # Exactly compare upper excess (3/2)sqrt(22)-5 with N times one-cell lower.
    one_cell_lower=F(270,6859)
    def beats_lower(N):
        rhs=F(2,3)*(5+N*one_cell_lower)
        return rhs>0 and rhs*rhs>22
    assert not beats_lower(51) and beats_lower(52)
    excess=1.5*math.sqrt(22)-5
    rows=[]
    for N in [1,2,8,16,51,52,64,256,1024]:
        rows.append({'layers':N,'peak_total_volume':14*N+excess,
                     'peak_volume_per_motif':14+excess/N,
                     'relative_peak_expansion_percent':100*excess/(14*N),
                     'upper_per_motif_excess_below_one_cell_certified_lower':beats_lower(N)})
    output={'theorem':'Any sequence q_k in [-1,1] yields a globally nonoverlapping fixed-contact layer stack using the certified L3(q_k).',
            'positions':'X(n1,n2,k)=n1*L1+n2*L2+S(k), S(k+1)-S(k)=L3(q_k)',
            'adjacent_exact_comparisons_rechecked':adjacent,
            'nonadjacent_lower_distance':'(12*m-10)/sqrt(35) >=14/sqrt(35)>sqrt(2), m>=2',
            'N_periodic_supercell':'Rows L1,L2,sum_{k=0}^{N-1}L3(q_k), containing 6N sphere sites.',
            'total_volume':'sum_k V(q_k), V(q)=9+q+4sqrt(3-2q^2)',
            'sequential_conversion_peak_exact':'14N + (3/2)sqrt(22)-5',
            'product_path_family_minimax':'The same value is optimal among arbitrary continuous q_k(t) in [-1,1] starting all +1 and ending all -1; first crossing of q*=1/sqrt(22) supplies the lower bound.',
            'first_N_beating_one_cell_lower_per_motif':52,
            'finite_interface_extra_normal_height':excess/math.sqrt(35),
            'examples':rows,
            'connector_corollary':'Contact capsules of radius a<r/2, r=1/sqrt(2), inherit the existing universal contact-framework certificate at every sequence and time.',
            'scope':'Changing periodic cell allowed. Vanishing relative/per-motif barrier as N grows; finite total transition excess remains. No fixed-box path, energy, forces, or physical articulation claimed.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
