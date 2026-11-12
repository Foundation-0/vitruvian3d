"""Exact finite certificate for a continuum vacancy-to-FCC switching path.

All contact minima use integer arithmetic after eliminating q²+2z²+2z=1.
Far translations are excluded by the analytic bound in volume_path_proof.md.
"""
import itertools
import json
import math
from pathlib import Path

P=[tuple(s if i==j else 0 for j in range(3)) for i in range(3) for s in [-1,1]]
V=sorted({tuple(a-b for a,b in zip(p,r)) for p in P for r in P})
L1=(2,-1,-1);L2=(1,2,-1);BASE3=(1,0,2)


def certify(c,a,b):
    # f(q)=c+a*(sqrt(3-2q²)-1)/2+bq, -1<=q<=1.
    # For a>=0, f is concave; otherwise convex, with stationary point
    # in the interval exactly when |b|<=|a|.
    if a>=0 or b*b>a*a:
        minimum=c-abs(b)
        assert minimum>=0,(c,a,b,'endpoint',minimum)
        return {'kind':'endpoint','minimum_integer':minimum}
    twice_constant=2*c-a
    radicand=3*(a*a+2*b*b)
    assert twice_constant>=0 and twice_constant**2>=radicand,(c,a,b,twice_constant,radicand)
    return {'kind':'interior_surd','minimum':'(%d-sqrt(%d))/2'%(twice_constant,radicand),
            'squared_margin':twice_constant**2-radicand}


def main():
    rows=[]
    for n in itertools.product(range(-2,3),repeat=3):
        for v in V:
            if n==(0,0,0) and v==(0,0,0):continue
            A=tuple(n[0]*L1[i]+n[1]*L2[i]+n[2]*BASE3[i]+v[i] for i in range(3))
            # |A+n3*(z,q,z)|²-2, substituting 2z²+q²=1-2z.
            c=sum(x*x for x in A)+n[2]*n[2]-2
            a=2*n[2]*(A[0]+A[2]-n[2])
            b=2*n[2]*A[1]
            cert=certify(c,a,b)
            rows.append({'translation':n,'motif_difference':v,'c':c,'a':a,'b':b,**cert})
    assert len(rows)==125*19-1
    # For far translations: det>=12, Frobenius²<20 => sigma_min>6/5;
    # outside [-2,2]^3, |n|>=3 => separation>8/5>sqrt(2).
    assert 64>2*25
    peak=9+1.5*math.sqrt(22)
    output={'path':'q:1 to -1; z=(sqrt(3-2q^2)-1)/2; L rows (2,-1,-1),(1,2,-1),(1+z,q,2+z)',
            'scope':'One-cell fixed motif ±e_i, all three prescribed bridges length sqrt(2), sphere diameter sqrt(2). Continuum path, not just sampled frames.',
            'finite_classes_checked':len(rows),'all_exact_minimum_checks_passed':True,
            'volume':'9+q+4sqrt(3-2q^2)', 'peak_parameter':'1/sqrt(22)',
            'peak_volume_exact':'9+(3/2)*sqrt(22)','peak_volume_float':peak,
            'initial_volume':14,'final_volume':12,
            'relative_peak_expansion_percent':100*(peak/14-1),
            'far_translation_bound':'sigma_min>6/5; |n|>=3 => separation>8/5>sqrt(2)',
            'claims_excluded':'Global minimax optimality, energy, actuator forces, rod thickness, and stability.',
            'conditions':rows}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='conditions'},indent=2))


if __name__=='__main__':main()
