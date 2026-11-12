"""Independent symbolic audit. Does not import the discovery implementation."""
from itertools import product, permutations
from math import comb
from pathlib import Path
import json
import sympy as S
u=S.symbols('u'); rt=S.sqrt(2)
C=S.Matrix([[0,1,1],[-1,0,1],[-1,-1,0]])
P=[S.Matrix(1,3,list(x)) for x in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]]
D={tuple(p-q) for p in P for q in P}
Ns=list(product(range(-1,2),repeat=3))
cache={}; minimum=None; negative=[]; count=0
for nn in Ns:
 n=S.Matrix(1,3,nn)
 for vv in D:
  if nn==(0,0,0) and vv==(0,0,0):continue
  v=S.Matrix(1,3,vv)
  # Derive directly from rationalized cos(theta), sin(theta), not supplied quartic.
  coords=n*((2*(1+u*u)+rt*(1-u*u))*S.eye(3)+2*u*C)+(1+u*u)*v
  expr=S.Poly(S.expand((coords*coords.T)[0]-2*(1+u*u)**2),u)
  co=tuple(expr.nth(i) for i in range(5))
  if co not in cache:
   bs=[S.expand(sum(co[i]*S.Rational(comb(k,i),comb(4,i)) for i in range(k+1))) for k in range(5)]
   assert all(b.is_nonnegative is not None for b in bs)
   cache[co]=bs
  bs=cache[co]
  for k,b in enumerate(bs):
   if not b.is_nonnegative:negative.append(dict(n=nn,v=vv,k=k,b=str(b)))
  count+=1
assert not negative
L=2*S.eye(3)+C; inv=L.inv(); residue=S.Matrix([5,1,3]); assert L*residue==S.Matrix([14,0,0])
# Membership independently checked by exact inverse, against parity+residue predicates.
checks=0
for z in product(range(-4,5),repeat=3):
 z=S.Matrix(1,3,z)
 member=any(all(x.q==1 for x in (z-p)*inv) for p in P)
 claimed=(sum(z)%2==1 and int((z*residue)[0])%7!=0)
 assert member==claimed
 checks+=1
# Independently form conjugacy orbits by graph relabeling plus coordinate gauge signs.
matrices=[]
for bits in product([-1,1],repeat=6):
 M=S.zeros(3);it=iter(bits)
 for i in range(3):
  for j in range(3):
   if i!=j:M[i,j]=next(it)
 matrices.append(M)
def key(M):return tuple(M)
unseen=set(map(key,matrices));orbits=[]
while unseen:
 a=min(unseen);M=S.Matrix(3,3,a);orbit=set()
 for perm in permutations(range(3)):
  for signs in product([-1,1],repeat=3):
   B=S.Matrix(3,3,lambda i,j:signs[i]*signs[j]*M[perm[i],perm[j]])
   orbit.add(key(B))
 unseen-=orbit
 # Find exact endpoint coincidences with distinct motif representatives.
 E=2*S.eye(3)+M; collisions=[]
 if E.det()!=0:
  Ei=E.inv()
  for i,p in enumerate(P):
   for j,q in enumerate(P):
    if i==j:continue
    shift=(p-q)*Ei
    if all(x.q==1 for x in shift):collisions.append(dict(p=list(map(int,p)),q=list(map(int,q)),n=list(map(int,shift))))
 else:collisions=[dict(singular=True)]
 orbits.append(dict(size=len(orbit),matrix=[list(map(int,M.row(i))) for i in range(3)],det=int(E.det()),endpoint_coincidences=collisions))
assert sorted((o['size'],o['det']) for o in orbits)==sorted([(4,4),(24,6),(12,8),(4,0),(12,12),(8,14)])
assert sum(o['size'] for o in orbits if not o['endpoint_coincidences'])==20
roots=[z for z in product([-1,0,1],repeat=3) if sum(x*x for x in z)==2]
root_residues=[int((S.Matrix(1,3,z)*residue)[0])%7 for z in roots]
assert sorted(root_residues)==[1,1,2,2,3,3,4,4,5,5,6,6]
report=dict(independent_method='SymPy direct coordinate polynomial expansion; degree-four Bernstein conversion through exact binomial rationals; no discovery code imported.',skew_continuum=dict(images=27,distance_conditions=count,distinct_polynomials=len(cache),coefficients_checked=5*count,negative_coefficients=negative,far_image_proof='C skew => |nL|²=d²|n|²+t²|nC|² >=4|n|². For |n|>=2, |nL+v|>=2>sqrt2; remaining integer n lie in[-1,1]^3.'),endpoint=dict(determinant=int(L.det()),fcc_index=int(L.det()/2),coset_inverse_checks=checks,residue_vector=[5,1,3],occupied_condition='sum(z) odd and dot(z,[5,1,3]) nonzero mod7',density='sqrt(2)*pi/7',fcc_neighbor_step_residues=sorted(root_residues),contacts_per_occupied_site=10,vacancy_neighbors_occupied=12,saturated=False),independent_orbits=orbits,limits='Skew full continuum verified independently. Full-fcc continuum accepted from repo proof plus orbit equivalence, not recertified here. Finite inverse membership checks supplement exact index/coset argument, rather than prove infinite equality by enumeration.')
# Exact identification with Koch-Fischer International Tables C, Table9.1.1.2,type20.
hex_basis=S.Matrix([[-2,1,3],[-1,-3,-2],[2,-2,2]])
origin=S.Matrix([[2,0,-1]])
def mod1(t):return tuple(x%1 for x in t)
actual={mod1((p+k*L.row(0)-origin)*hex_basis.inv()) for k in range(3) for p in P}
reference=set();x,y=S.Rational(3,7),S.Rational(1,7)
for k in range(3):
 for sign in [-1,1]:
  for a,b,c in [(0,0,0),(S.Rational(2,3),S.Rational(1,3),S.Rational(1,3)),(S.Rational(1,3),S.Rational(2,3),S.Rational(2,3))]:
   reference.add(mod1([sign*x+a,sign*y+b,c]))
 x,y=-y,x-y
assert len(actual)==len(reference)==18 and actual==reference
assert hex_basis*hex_basis.T==S.Matrix([[14,-7,0],[-7,14,0],[0,0,12]])
report['prior_art_identification']=dict(source='https://onlinelibrary.wiley.com/iucr/itc/Cb/ch9o1v0001/',source_detail='Koch and Fischer, International Tables for Crystallography C(2006),Table9.1.1.2,type20; text explicitly describes deletion of1/7 of fcc sites.',hex_basis=[list(map(int,hex_basis.row(i))) for i in range(3)],origin=list(map(int,origin)),gram=[[14,-7,0],[-7,14,0],[0,0,12]],reference_wyckoff='R-3 18f, x=3/7,y=1/7,z=0,c/a=sqrt42/7',actual_sites=sorted([list(map(str,z)) for z in actual]),exact_reference_orbit_equal=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
