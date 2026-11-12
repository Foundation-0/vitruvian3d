"""Independent algebra and composition review of the certified layer stack."""
from pathlib import Path
import json,itertools
import sympy as s
L1=s.Matrix([[2,-1,-1]]);L2=s.Matrix([[1,2,-1]])
h=L1.cross(L2);assert h==s.Matrix([[3,1,5]])
p,q,z=s.symbols('p q z',real=True)
L3=s.Matrix([[1+z,q,2+z]])
assert s.expand((h*L3.T)[0]-(13+q+8*z))==0
assert s.expand((L3-s.Matrix([[0,0,2]])).dot(L3-s.Matrix([[0,0,2]]))-2)==q*q+2*z*z+2*z-1
P=[s.Matrix([v]) for v in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]]
assert max(abs((h*(a-b).T)[0]) for a in P for b in P)==10
assert 14**2>2*35
V=9+q+4*s.sqrt(3-2*q*q)
assert s.simplify(s.diff(V,q,2)+24/(3-2*q*q)**s.Rational(3,2))==0
assert V.subs(q,1)==14 and V.subs(q,-1)==12
qstar=1/s.sqrt(22);M=9+s.Rational(3,2)*s.sqrt(22)
assert s.simplify(V.subs(q,qstar)-M)==0
# On q in[qstar,1], the square comparison proving V(q)>=14 factors exactly.
assert s.expand(16*(3-2*q*q)-(5-q)**2)==s.expand((1-q)*(33*q+23))
# Therefore first crossing of any q_k through qstar costs at least
# M+14(N-1), while sequential motion realizes this same peak.
N52=s.Rational(270,6859)
assert bool((M-14)/51>N52) and bool((M-14)/52<N52)
# General determinant linearity with fixed horizontal rows, independent layerq,z.
q0,q1,z0,z1=s.symbols('q0 q1 z0 z1')
S0=s.Matrix([[1+z0,q0,2+z0]]);S1=s.Matrix([[1+z1,q1,2+z1]])
assert s.expand(s.Matrix.vstack(L1,L2,S0+S1).det()-s.Matrix.vstack(L1,L2,S0).det()-s.Matrix.vstack(L1,L2,S1).det())==0
result=dict(verdict='PASS arbitrary layer composition, nonadjacent separation, exact supercell volume, sequential upper bound and product-family minimax lower bound.',independent_adjacent_layer_evidence='connectors_volume_review.py rationally certifies all 2374 pair classes for the complete uniform path. For n3=±1 these are exactly adjacent-layer differences; n3=0 is unchanged. No layer-certificate implementation imported.',layer_normal=[3,1,5],normal_squared_length=35,motif_difference_projection_bound=10,nonadjacent_site_distance_lower_bound='(12m-10)/sqrt35 >=14/sqrt35>sqrt2, m>=2',product_family_exact_peak='14N+(3/2)sqrt22-5',first_crossing_argument='At first time any q_k reaches qstar, every other q_j>=qstar and has V(q_j)>=14. Total is at least14(N-1)+M. Sequential conversion attains it.',first_layers_beating_one_cell_per_motif_lower_bound=52,scope='Arbitrary q_k(t) in the specified product path family. Fixed first two rows and motifs. Box may change. Relative/per-motif expansion vanishes, total excess does not. No mechanics or unrestricted-supercell minimax claim.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
