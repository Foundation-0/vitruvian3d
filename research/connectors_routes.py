"""Optimize independently established skew winding routes for positive clearance."""
from pathlib import Path
import itertools,json
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
R=1/np.sqrt(2);theta=np.deg2rad(46.1113417712)
C=np.array([[0,1,1],[-1,0,1],[-1,-1,0]])
L=(2+np.sqrt(2)*np.cos(theta))*np.eye(3)+np.sin(theta)*C
P=np.r_[np.eye(3),-np.eye(3)]
N=np.array(list(itertools.product(range(-4,5),repeat=3)))
sites=((N@L)[:,None,:]+P).reshape(-1,3)
tree=cKDTree(sites)
def ds_segment(x,y,ss):
 d=y-x;rel=ss-x;u=np.clip(rel@d/(d@d),0,1)
 return np.einsum('ij,ij->i',rel-u[:,None]*d,rel-u[:,None]*d)
def clearance(x,y):
 mid=(x+y)/2;ss=sites[tree.query_ball_point(mid,np.linalg.norm(y-x)/2+2)]
 return np.sqrt(ds_segment(x,y,ss).min())
def main():
 data=json.loads(Path(__file__).with_name('transport_independent.json').read_text())
 case=next(r for r in data['results'] if r['case']=='skew' and r['resolution']==96)
 results=[]
 for witness in case['witness_routes']:
  route=np.array(witness['fractional_lifted_route'])@L
  # Greedy visible shortcuts, preserving original loop endpoints and homology.
  selected=[route[0]];i=0
  while i<len(route)-1:
   nexti=i+1
   for j in range(len(route)-1,i,-1):
    if clearance(route[i],route[j])>np.sqrt(2)+1e-7:nexti=j;break
   selected.append(route[nexti]);i=nexti
  initial=np.array(selected);k=len(initial)-2
  # Include all sites within a generous movement bound of each initial segment.
  candidates=[]
  for x,y in zip(initial,initial[1:]):
   near=sites[tree.query_ball_point((x+y)/2,np.linalg.norm(y-x)/2+2)]
   candidates.append(near[ds_segment(x,y,near)<(np.sqrt(2)+.55)**2])
  def unpack(z):return np.r_[initial[0:1],z[:-1].reshape(-1,3),initial[-1:]]
  def con(z):
   rr=unpack(z);D=z[-1]
   return np.concatenate([ds_segment(x,y,ss)-D*D for x,y,ss in zip(rr,rr[1:],candidates)])
  z=np.r_[initial[1:-1].reshape(-1),min(clearance(x,y) for x,y in zip(initial,initial[1:]))]
  bounds=[(x-.2,x+.2) for x in z[:-1]]+[(np.sqrt(2),np.sqrt(2)+.1)]
  opt=minimize(lambda z:-z[-1],z,method='SLSQP',bounds=bounds,constraints=[dict(type='ineq',fun=con)],options=dict(maxiter=400,ftol=1e-12))
  rr=unpack(opt.x);minimum=min(clearance(x,y) for x,y in zip(rr,rr[1:]))
  assert minimum>np.sqrt(2)
  fractional=rr@np.linalg.inv(L)
  assert np.max(np.abs(fractional[-1]-fractional[0]-np.array(witness['winding'])))<1e-12
  omitted_fractional_gap=min(5-fractional.max(),5+fractional.min())
  omitted_lower=np.linalg.svd(L,compute_uv=False)[-1]*omitted_fractional_gap-1
  assert omitted_lower>minimum
  all_site_minimum=min(float(np.sqrt(ds_segment(x,y,sites).min())) for x,y in zip(rr,rr[1:]))
  assert abs(all_site_minimum-minimum)<1e-12
  result=dict(winding=witness['winding'],segments=len(rr)-1,solver_success=bool(opt.success),solver_message=str(opt.message),minimum_site_distance=float(minimum),probe_clearance_above_equal_radius=float(minimum-np.sqrt(2)),euclidean_route=rr.tolist(),fractional_route=(rr@np.linalg.inv(L)).tolist(),candidate_constraint_minimum=float(con(opt.x).min()),omitted_site_distance_lower_bound=float(omitted_lower),independent_full_image_projection_minimum=all_site_minimum)
  results.append(result);print(json.dumps(result),flush=True)
  Path(__file__).with_suffix('.json').write_text(json.dumps(dict(method='Greedy visibility simplification followed by constrained waypoint optimization; analytic segment clearance. Local numerical optimization, no optimality claim.',results=results),indent=2)+'\n')

if __name__=="__main__":
 main()
