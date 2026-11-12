"""Independent symbolic audit of the finite volume projection proof."""
import json
from pathlib import Path
import sympy as s
L=s.Matrix([[2,-1,-1],[1,2,-1],[1,1,2]])
x=s.symbols('d11 d12 d13 d21 d22 d23 d31 d32 d33');D=s.Matrix(3,3,x)
q1,q2,q3=s.symbols('q1 q2 q3',nonnegative=True)
qs=[q1,q2,q3]
bridge={D[0,2]:-D[0,1]+q1,D[1,2]:D[1,0]+q2,D[2,1]:-D[2,0]-q3}
S=s.Matrix([[q1,0,-q1],[0,q2,-q2],[0,q3,q3]])
A=(D+S).subs(bridge)
expected=s.Matrix([[D[0,0]+q1,D[0,1],-D[0,1]],[D[1,0],D[1,1]+q2,D[1,0]],[D[2,0],-D[2,0],D[2,2]+q3]])
assert A==expected
# Six genuine distance constraints selected directly from contact vectors.
contacts=[(0,s.Matrix([[1,-1,0]])),(0,s.Matrix([[1,0,-1]])),(1,s.Matrix([[0,1,-1]])),(1,s.Matrix([[1,1,0]])),(2,s.Matrix([[0,1,1]])),(2,s.Matrix([[1,0,1]]))]
P=[s.Matrix([p]) for p in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]]
diffs={tuple(p-t) for p in P for t in P}
g=[]
for i,d in contacts:
 assert (d*d.T)[0]==2 and tuple(d-L.row(i)) in diffs
 g.append(s.expand((d*D.row(i).T)[0]+qs[i]).subs(bridge))
# Nonnegative linear combinations of contact inequalities plus nonnegative q's.
cone=[A[0,0]-D[0,1],A[0,0]+D[0,1],A[1,1]-D[1,0],A[1,1]+D[1,0],A[2,2]-D[2,0],A[2,2]+D[2,0]]
assert [s.expand(a-b) for a,b in zip(cone,[g[0],g[1]+q1,g[2]+q2,g[3],g[4]+q3,g[5]])]==[0]*6
cof=L.cofactor_matrix();assert cof==s.Matrix([[5,-3,-1],[1,5,-3],[3,1,5]])
cofS=s.expand(sum(cof[i,j]*S[i,j] for i in range(3) for j in range(3)))
assert cofS==6*q1+8*q2+6*q3
cofA=s.expand(sum(cof[i,j]*A[i,j] for i in range(3) for j in range(3)))
assert cofA==5*(A[0,0]+A[1,1]+A[2,2])-2*D[0,1]-2*D[1,0]+2*D[2,0]
# Independent determinant expansion checks the three homogeneous orders.
t=s.symbols('t');detpoly=s.Poly((L+t*D).det(),t)
assert detpoly.nth(0)==14
assert s.expand(detpoly.nth(1)-sum(cof[i,j]*D[i,j] for i in range(3) for j in range(3)))==0
assert detpoly.nth(3)==D.det()
Q=detpoly.nth(2)
assert s.expand(Q-(s.Matrix.vstack(D.row(0),D.row(1),L.row(2)).det()+s.Matrix.vstack(D.row(0),L.row(1),D.row(2)).det()+s.Matrix.vstack(L.row(0),D.row(1),D.row(2)).det()))==0
assert all((L.row(i)*L.row(i).T)[0]==6 for i in range(3))
assert s.Rational(69,2)+9*s.sqrt(6)<57
assert 3*s.sqrt(3)<6
r=s.Rational(1,38);lower=3*r-57*r*r-6*r**3
assert lower==s.Rational(270,6859)
result=dict(review='PASS finite projection, actual contact witnesses, tangent inequality, determinant expansion and exact lower bound. No author implementation imported.',cofactor=[list(map(int,cof.row(i))) for i in range(3)],projection_contact_slack=['0','q1','q2','0','q3','0'],determinant_quadratic_polynomial=str(s.expand(Q)),strict_absolute_barrier_lower_bound=str(lower),strict_peak_lower_bound=str(14+lower),limitations='Bound restricted to stated one-cell fixed-motif bridge family. No claim about energy, sharpness, larger cells or motif rotations.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

# Distinct continuum upper-path verification: rationalize the ellipse and check
# elementary integer quadratics, avoiding the author's surd minimization.
from itertools import product
import numpy as np
w=s.symbols('w',real=True)
den=4*w*w-8*w+6;Q=4*w*w-12*w+6;Z=4*w-4*w*w
assert s.expand(Q*Q+2*Z*Z+2*Z*den-den*den)==0
assert s.simplify((Q/den).subs(w,0))==1 and s.simplify((Q/den).subs(w,1))==-1
# Monotonicity checked by exact derivative numerator below.
derivative_numerator=s.expand(s.diff(Q,w)*den-Q*s.diff(den,w))
assert derivative_numerator==16*w*w-24
# The upper ellipse arc is traced once: denominator>0, Z>=0, q decreases.
Dpoly=np.array([6,-8,4],dtype=np.int64);Qpoly=np.array([6,-12,4],dtype=np.int64);Zpoly=np.array([0,4,-4],dtype=np.int64)
row1=(2,-1,-1);row2=(1,2,-1)
classes=0;interior_count=0;certificates=[]
for n in product(range(-2,3),repeat=3):
 for v in diffs:
  if n==(0,0,0) and v==(0,0,0):continue
  coords=[]
  for j in range(3):
   initial=n[0]*row1[j]+n[1]*row2[j]+v[j]
   extra=[Dpoly+Zpoly,Qpoly,2*Dpoly+Zpoly][j]
   coords.append(initial*Dpoly+n[2]*extra)
  polynomial=sum((np.convolve(c,c) for c in coords),start=np.zeros(5,dtype=np.int64))-2*np.convolve(Dpoly,Dpoly)
  remainder=polynomial.copy();g=np.zeros(3,dtype=np.int64)
  for power in [4,3,2]:
   assert remainder[power]%4==0
   factor=remainder[power]//4;g[power-2]=factor
   remainder[power-2:power+1]-=factor*Dpoly
  assert not np.any(remainder)
  cc,bb,aa=map(int,g)
  assert cc>=0 and aa+bb+cc>=0
  interior=(aa>0 and 0 < -bb < 2*aa)
  if interior:
   assert 4*aa*cc-bb*bb>=0
   interior_count+=1
  classes+=1
  certificates.append(dict(n=n,v=list(map(int,v)),quadratic_coefficients_low_to_high=[cc,bb,aa],interior_vertex=interior))
assert classes==2374
result['upper_path_independent_rational_certificate']=dict(parameter='w in[0,1], D=4w²-8w+6, q=(4w²-12w+6)/D, z=(4w-4w²)/D',finite_comparisons=classes,interior_quadratic_vertex_checks=interior_count,method='Expand actual rational coordinate distance numerators, divide by positive denominator, check endpoint values and integer quadratic discriminant at any interior minimum.',all_passed=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
Path(__file__).with_name('connectors_volume_path_conditions.json').write_text(json.dumps(certificates,indent=2)+'\n')
print(json.dumps(result['upper_path_independent_rational_certificate'],indent=2))
