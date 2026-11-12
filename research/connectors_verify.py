"""Replay saved improved routes without running their numerical optimizer.

Exact connector inequalities are checked separately by connectors_exact.py.
This route replay uses floating-point distances and claims no interval proof.
"""
from pathlib import Path
import itertools,json,math
import numpy as np


def main():
 data=json.loads(Path(__file__).with_name('connectors_routes.json').read_text())
 theta=math.radians(46.1113417712)
 C=np.array([[0,1,1],[-1,0,1],[-1,-1,0]])
 L=(2+math.sqrt(2)*math.cos(theta))*np.eye(3)+math.sin(theta)*C
 P=np.r_[np.eye(3),-np.eye(3)]
 N=np.array(list(itertools.product(range(-4,5),repeat=3)))
 sites=((N@L)[:,None,:]+P).reshape(-1,3)
 root=None;winds=[];minimum=float('inf');segments=0
 for loop in data['results']:
  points=np.array(loop['euclidean_route']);wind=np.array(loop['winding'])
  if root is None:root=points[0]
  assert np.linalg.norm(points[0]-root)<1e-12
  assert np.linalg.norm(points[-1]-points[0]-wind@L)<1e-12
  frac=points@np.linalg.inv(L)
  omitted_gap=min(5-float(frac.max()),5+float(frac.min()))
  omitted_bound=np.linalg.svd(L,compute_uv=False)[-1]*omitted_gap-1
  assert omitted_bound>2
  for a,b in zip(points,points[1:]):
   delta=b-a;vv=float(delta@delta);assert vv>0
   rel=sites-a;t=np.clip(rel@delta/vv,0,1)
   nearest=float(np.linalg.norm(rel-t[:,None]*delta,axis=1).min())
   assert nearest>math.sqrt(2)+.0035
   minimum=min(minimum,nearest);segments+=1
  winds.append(wind)
 assert abs(round(np.linalg.det(winds)))==1
 result=dict(saved_loops=len(winds),segments=segments,minimum_site_distance=minimum,minimum_margin_above_equal_probe=minimum-math.sqrt(2),minimum_accepted_margin=.0035,route_family='Matched-volume skew configuration at 46.1113417712 degrees',arithmetic='Floating-point full-image projection with analytic omitted-image bound; not an interval certificate.')
 print(json.dumps(result,indent=2))


if __name__=='__main__':main()
