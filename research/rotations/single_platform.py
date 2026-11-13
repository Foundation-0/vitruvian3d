"""Independent exact audit of one rotating six-bar octahedral platform.

No existing research implementation is imported. All universal conclusions
use the accompanying analytic proof, with symbolic identities checked here.
"""
from pathlib import Path
import json
import sympy as s
L=s.Matrix([[2,-1,-1],[1,2,-1],[1,1,2]])
b=s.Matrix([-1,1,-1]);K=s.Matrix([[0,1,1],[-1,0,1],[-1,-1,0]])
A=L.T-s.eye(3)
assert A==s.eye(3)+K and K.T==-K and K*b==s.zeros(3,1)
assert A*A.T==4*s.eye(3)-b*b.T
h,t=s.symbols('h t',real=True)
c=(1-h*h)/(1+h*h);ss=2*h/(1+h*h)
R=c*s.eye(3)+(1-c)*(b*b.T)/3+ss*K/s.sqrt(3)
assert s.simplify(R.T*R-s.eye(3))==s.zeros(3)
assert s.simplify(R*b-b)==s.zeros(3,1)
T=8*h*(s.sqrt(3)-h)/(9*(1+h*h))
u=t*b
for i in range(3):
 ei=s.eye(3).col(i);ai=A.col(i)
 assert s.simplify(u.dot(R*ei-ai))==0
 for sign in [-1,1]:
  d=u+sign*(R*ei-ai)
  gap=s.Poly(s.expand(d.dot(d)-2),t)
  assert gap.degree()==2
  assert s.simplify(gap.nth(2)*T+gap.nth(1)*t+gap.nth(0))==0
# At zero translation the three orientation equations have invertible Jacobian.
M=s.Matrix.vstack(*[(s.eye(3).col(i).cross(A.col(i))).T for i in range(3)])
assert M.det()==2
p=s.Matrix([1,0,0]);q=s.Matrix([1,1,-1])
assert q==L.row(1).T-s.Matrix([0,1,0]) and (p-q).dot(p-q)==2
for sign in [-1,1]:
 d=u+sign*(R*p-q)
 expected=3*t*t-sign*4*t+s.Rational(8,3)*(1-c)
 assert s.simplify(d.dot(d)-2-expected)==0
# Both expressions are nonnegative for0<=h<=1/sqrt3, and strictly so
# before its upper endpoint where stated in the proof.
slack=s.factor(3*T-s.Rational(8,3)*(1-c))
assert s.simplify(slack-8*h*(s.sqrt(3)-3*h)/(3*(1+h*h)))==0
max_translation_slack=s.factor(s.Rational(4,9)-T)
assert s.simplify(max_translation_slack-4*(s.sqrt(3)*h-1)**2/(9*(1+h*h)))==0
# Smooth local parameterization by signed axial translation t.
h_local=9*t*t/(4*s.sqrt(3)+s.sqrt(48-72*t*t-81*t**4))
assert s.limit(h_local/(t*t),t,0)==3*s.sqrt(3)/8
assert s.simplify((8+9*t*t)*h_local*h_local-8*s.sqrt(3)*h_local+9*t*t)==0
result=dict(platform='L rows(2,-1,-1),(1,2,-1),(1,1,2), fixed neighbors; bodyvertices u+R(+-e_i), fixed bridge anchors+-a_i with a_i=L_i-e_i.',bar_only_self_motion=True,axis=[-1,1,-1],rational_rotation_parameter='h=tan(phi/2), c=(1-h²)/(1+h²), sin(phi)=2h/(1+h²)',translation='u=t*(-1,1,-1), t²=8h(sqrt3-h)/(9(1+h²))',bar_family_domain='0<=h<=sqrt3, corresponding to0<=phi<=120degrees; t has either sign.',local_smooth_branch='h=9t²/[4sqrt3+sqrt(48-72t²-81t⁴)], |t|<=2/3, 0<=phi<=60degrees',small_translation_rotation='phi=(3sqrt3/4)t²+O(t⁴)',classification='Any solution with u!=0 has u parallel to axis and R fixes axis. At u=0, identity orientation is isolated by a three-equation Jacobian of determinant2.',nonbar_contact_pair=dict(body_vertex=[1,0,0],fixed_site=[1,1,-1],fixed_site_label='n=(0,1,0), motif point(0,-1,0)',opposite_pair='body vertex(-1,0,0), fixed site(-1,-1,1)',minimum_squared_gap='3t²-4|t|+(8/3)(1-cos(phi))',strict_collision_interval='0<phi<60degrees on the identity branch. Then (8/3)(1-cos(phi))<=3t² and0<|t|<2/3, giving gap<=6t²-4|t|<0.'),hard_sphere_local_self_motion=False,scope='One moving rigid octahedral motif, all neighbors fixed, all six bridges preserved. Local obstruction at identity only; not a multi-motif rotational rigidity theorem. No mechanical or historical novelty claim.',symbolic_checks_passed=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
