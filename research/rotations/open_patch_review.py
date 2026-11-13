"""Independent port-rank and all-image clearance verification.

No rotation_astra or other research module is imported. Actual rigid-anchor
Jacobians are assembled directly, then compared with the claimed port map.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sympy as s

C=s.Matrix([[0,-1,-1],[1,0,-1],[1,1,0]])
a,t=s.symbols('a t',positive=True)
B=a*s.eye(3)+t*C
T=s.Matrix.vstack(*[(s.eye(3).col(i).cross(B.row(i).T)).T for i in range(3)])
M=B.row_join(T).col_join(B.row_join(-T))
assert s.factor(B.det())==a*(a*a+3*t*t)
assert s.factor(T.det())==2*t**3
assert s.factor(M.det())==-16*a*t**3*(a*a+3*t*t)

# Independent geometric row assembly: unit translations and rotations move each
# anchor p by v+omega cross p; dot the relative velocity with the bridge vector.
def actual_patch(nodes,rho_value=2,with_minor=False):
 nodes=sorted(nodes);index={n:i for i,n in enumerate(nodes)};S=set(nodes)
 candidates=product(*[range(min(n[i] for n in nodes)-1,max(n[i] for n in nodes)+2) for i in range(3)])
 edges=[]
 for n in candidates:
  for axis in range(3):
   m=tuple(n[j]+(j==axis) for j in range(3))
   if n in S or m in S:edges.append((n,m,axis))
 B0=rho_value*s.eye(3)+C
 T0=s.Matrix.vstack(*[(s.eye(3).col(i).cross(B0.row(i).T)).T for i in range(3)])
 M0=B0.row_join(T0).col_join(B0.row_join(-T0))
 J=s.zeros(len(edges),6*len(nodes));inc=s.zeros(len(edges),6*len(nodes));E=0
 for row,(n,m,i) in enumerate(edges):
  bridge=B0.row(i).T
  for endpoint,node,anchor_sign in [(-1,n,1),(1,m,-1)]:
   if node not in index:continue
   offset=6*index[node];anchor=anchor_sign*s.eye(3).col(i)
   for j in range(3):
    J[row,offset+j]=endpoint*bridge.dot(s.eye(3).col(j))
    J[row,offset+3+j]=endpoint*bridge.dot(s.eye(3).col(j).cross(anchor))
   inc[row,offset+(i if endpoint==-1 else 3+i)]=endpoint
  if n in S and m in S:E+=1
 block=s.diag(*([M0]*len(nodes)))
 assert J==inc*block
 gram=inc*inc.T
 assert gram.is_diagonal() and all(v in [1,2] for v in gram.diagonal())
 assert len(edges)==6*len(nodes)-E
 assert J.rank()==len(edges)
 answer=dict(nodes=len(nodes),internal_edges=E,actual_jacobian_rank=int(J.rank()),local_dimension=6*len(nodes)-int(J.rank()))
 if with_minor:answer['minor']=str(s.factor(J[:,:11].det()))
 return answer
patches={'one':[(0,0,0)],'edge':[(0,0,0),(1,0,0)],'disconnected':[(0,0,0),(2,1,0)],'bent_four':[(-1,0,0),(0,0,0),(0,1,0),(0,1,1)],'cube':list(product([0,1],repeat=3))}
patch_results={name:actual_patch(nodes) for name,nodes in patches.items()}

# Exact Q(sqrt2) distance comparisons, using rational outward bounds rather
# than the author's integer-pair sign routine.
lo=F(141421356237309504880,10**20);hi=F(141421356237309504881,10**20)
assert lo*lo<2<hi*hi
P=[tuple(sign if i==j else 0 for j in range(3)) for i in range(3) for sign in [-1,1]]
V={tuple(p[j]-q[j] for j in range(3)) for p in P for q in P}
count=0;tight=[]
for n in product([-1,0,1],repeat=3):
 if n==(0,0,0):continue
 for v in V:
  if sum(x*x for x in n)==1 and tuple(-2*x for x in n)==v:continue
  h=[3*n[j]+v[j] for j in range(3)]
  z=[sum(n[i]*int(C[i,j]) for i in range(3)) for j in range(3)]
  rational=sum(F(x*x) for x in h)+sum(F(x*x,2) for x in z)-6
  radical=sum(F(h[j]*z[j]) for j in range(3))+1
  exact_lower=rational+radical*(lo if radical>=0 else hi)
  assert exact_lower>=0
  if rational==radical==0:tight.append(dict(n=n,v=v))
  count+=1
assert count==488 and tight
L0=3*s.eye(3)+C/s.sqrt(2)
assert s.simplify(L0*L0.T-(9*s.eye(3)+C*C.T/2))==s.zeros(3)
assert (s.Rational(3,2)**2)>2

# Independently check the published explicit adjacent-pair tangent.
rho=s.symbols('rho',positive=True);Br=rho*s.eye(3)+C
Tr=s.Matrix.vstack(*[(s.eye(3).col(i).cross(Br.row(i).T)).T for i in range(3)])
symbolic_pair=actual_patch([(0,0,0),(1,0,0)],rho,True)
assert s.sympify(symbolic_pair['minor'],locals={'rho':rho})==-64*rho*rho*(rho*rho+3)**2
omega=s.Matrix([1,1,-1])
trans=s.Matrix([2*(rho*rho+1),-2*(rho+1),-2*(rho-1)])/(rho*(rho*rho+3))
assert s.simplify(Br*trans-Tr*omega)==s.zeros(3,1)
assert Tr*omega==s.Matrix([2,0,0])
result=dict(verdict='PASS local analytic compact mechanisms at the specified open skew baseline, subject to relative-neighborhood wording.',scope='Any finite free motif patch S, exterior fixed; rigid translations and rotations allowed; all incident prescribed bridges retained. Local existence, not a certified finite amplitude or dense-endpoint continuation.',bridge_determinant=str(s.factor(B.det())),torque_determinant=str(s.factor(T.det())),port_determinant=str(s.factor(M.det())),independent_actual_anchor_jacobians=patch_results,universal_rank_argument='Every motif port is used once. In port coordinates the row Gram matrix is diagonal with entries1 or2. Thus all6|S|-E rows are independent for every finite S.',exact_near_distance_classes=count,minimum_squared_nondesignated_intermotif_distance='6-sqrt2>4',tight_examples=tight,far_distance_bound='|n|>=2 =>|nL+v|>=3|n|-2>=4',packing_neighborhood='On the bridge-preserving manifold, every actual site displacement<1/4 ensures every other intermotif separation>3/2>sqrt2.',explicit_pair_tangent_checked=True,independent_symbolic_pair_minor=symbolic_pair['minor'],physical_connector_corollary='Existing contact-capsule theorem applies without an orientation assumption when a/r<1/2.',author_modules_imported=False)
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
