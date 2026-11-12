"""Numerical periodic Voronoi void probe. Spherical obstacle sites only; no bars.
Source commit 2d374fff2a59cc9331566733043a4bb25c1040bf. Not a formal proof.
"""
import itertools,json
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi,cKDTree
P=np.r_[np.eye(3),-np.eye(3)]
C=np.array([[0,-1,-1],[1,0,-1],[1,-1,0]])
r=1/np.sqrt(2)
rows=[]
replication_checks=[]
for deg in [0,15,30,45,60,75,85,89,90]:
 t=np.deg2rad(deg); L=(2+np.sqrt(2)*np.cos(t))*np.eye(3)+np.sin(t)*C
 ns=np.array(list(itertools.product(range(-2,3),repeat=3)))
 sites=(ns@L)[:,None,:]+P
 sites=sites.reshape(-1,3)
 vor=Voronoi(sites); verts=vor.vertices; f=verts@np.linalg.inv(L)
 vv=verts[np.all((f>=-1e-8)&(f<1-1e-8),axis=1)]
 tree=cKDTree(sites); ds=tree.query(vv)[0]
 # Independent larger-image Voronoi comparison at all sampled angles.
 ns3=np.array(list(itertools.product(range(-3,4),repeat=3)))
 sites3=((ns3@L)[:,None,:]+P).reshape(-1,3)
 verts3=Voronoi(sites3).vertices; frac3=verts3@np.linalg.inv(L)
 vv3=verts3[np.all((frac3>=-1e-8)&(frac3<1-1e-8),axis=1)]
 max3=float(cKDTree(sites3).query(vv3)[0].max())
 assert abs(ds.max()-max3)<1e-8
 replication_checks.append(dict(angle=deg,maximum_difference=float(abs(ds.max()-max3))))
 order=np.argsort(-ds); best=vv[order[0]]
 maxima=[]
 for i in order:
  if ds[i] < ds[order[0]]-1e-7: break
  if not any(np.linalg.norm((vv[i]-v)@np.linalg.inv(L)-np.round((vv[i]-v)@np.linalg.inv(L)))<1e-7 for v in maxima): maxima.append(vv[i])
 rows.append(dict(angle=deg,largest_void_radius=float(ds.max()-r),ratio_to_obstacle_radius=float((ds.max()-r)/r),motif_origin_clearance=float(tree.query([0,0,0])[0]-r),cell_midpoint_clearance=float(tree.query(np.ones(3)@L/2)[0]-r),largest_void_fractional=(best@np.linalg.inv(L)).tolist(),largest_void_multiplicity=len(maxima),number_of_voronoi_vertices=len(vv)))
assert abs(rows[-1]['ratio_to_obstacle_radius']-(np.sqrt(2)-1))<1e-7
result=dict(method='Voronoi vertices in fundamental parallelepiped, sphere centers replicated over [-2,2]^3. Floating point, sampled angles; no path accessibility claim; bars have zero thickness and are omitted.',results=rows,replication_checks=replication_checks,central_cage=dict(room_radius=float(1-r),triangular_gate_radius=float(np.sqrt(2/3)-r),gate_ratio_to_obstacle_radius=float(2/np.sqrt(3)-1),scope='Exact upper bound for spherical probes exiting the rigid motif interior. Does not establish a periodic traversal route.'))
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
