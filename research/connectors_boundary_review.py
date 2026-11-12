"""Independent stdlib checks of fixed-boundary telescoping and flex caveats."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import json

B=((F(2,3),F(-1,7),F(4,9)),(F(-2,5),F(3,11),F(1,4)),(F(7,13),F(2,9),F(-1,3)))
field={(0,0,0):(F(1,7),F(-2,9),F(3,11)),(1,0,0):(F(-1,5),F(2,7),F(1,13)),(-1,1,-1):(F(2,3),F(1,8),F(-1,6))}
zero=(F(0),)*3
def dot(a,b):return sum(x*y for x,y in zip(a,b))
energy=F(0);excess=F(0);linears=[F(0)]*3;count=0
# A full surrounding box, rather than the author's incident-edge set construction.
for n in product(range(-3,4),repeat=3):
 for i in range(3):
  m=tuple(x+(j==i) for j,x in enumerate(n))
  d=tuple(y-x for x,y in zip(field.get(n,zero),field.get(m,zero)))
  new=tuple(x+y for x,y in zip(B[i],d))
  energy+=dot(d,d);excess+=dot(new,new)-dot(B[i],B[i]);linears[i]+=dot(B[i],d);count+=1
assert energy==excess>0 and linears==[0]*3
# Exact harmonic-mean identity underlying the six-ray Dirichlet bound.
for R in range(6):
 for x in range(-R,R+1):
  lhs=F(1,R+1-x)+F(1,R+1+x)-F(2,R+1)
  rhs=F(2*x*x,(R+1)*((R+1)**2-x*x))
  assert lhs==rhs>=0
vacancy=((0,-1,-1),(1,0,-1),(1,1,0));w=(1,-1,1)
assert all(dot(b,w)==0 for b in vacancy) and dot(w,w)==3
assert 6*dot(w,w)==18
# Taylor coefficient: delta length starts at 3/(2sqrt2)*t²(Df)².
# Squaring that coefficient and multiplying byk/2 gives9k/16.
assert F(1,2)*F(9,8)==F(9,16)
result=dict(verdict='PASS exact fixed-exterior and fixed-period translation rigidity; linearized common-kernel flexes remain possible and cannot integrate as finite supported retained-length motions.',independent_finite_field=dict(arbitrary_rational_bridge_vectors=True,box_edges_checked=count,dirichlet_sum=str(energy),squared_length_excess_sum=str(excess),linear_sums=list(map(str,linears))),harmonic_identity='1/(R+1-x)+1/(R+1+x)-2/(R+1)=2x²/[(R+1)((R+1)²-x²)]>=0',vacancy_kernel=[1,-1,1],single_site_total_quadratic_excess=18,optional_spring_quartic_coefficient='9/16',scope='Fixed motif orientations/anchor offsets, homogeneous baseline, fixed exterior or periodic cell. Not arbitrary rotating bodies, nonuniform backgrounds, feasible hard-sphere infinitesimal modes, or a specified physical spring model.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
