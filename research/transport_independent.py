"""Independent periodic uniform-grid transport check; no Voronoi or union-find."""
import itertools,json,time
from pathlib import Path
from collections import deque
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
P=np.r_[np.eye(3),-np.eye(3)]
CASES=[('original',40.,np.array([[0,-1,-1],[1,0,-1],[1,-1,0]])),('skew',46.1113417712,np.array([[0,1,1],[-1,0,1],[-1,-1,0]]))]
def winding_basis(vectors):
 b=[]
 for v in vectors:
  v=np.asarray(v,dtype=np.int64)
  if not np.any(v):continue
  if not b:b.append(v)
  elif len(b)==1 and np.any(np.cross(b[0],v)):b.append(v)
  elif len(b)==2 and np.dot(np.cross(b[0],b[1]),v)!=0:b.append(v);break
 return [v.tolist() for v in b]
def run(name,degree,C,m):
 start=time.time();t=np.deg2rad(degree);L=(2+np.sqrt(2)*np.cos(t))*np.eye(3)+np.sin(t)*C
 ns=np.array(list(itertools.product(range(-2,3),repeat=3)))
 sites=((ns@L)[:,None,:]+P).reshape(-1,3);tree=cKDTree(sites)
 sigma=np.linalg.svd(L,compute_uv=False)[-1]
 omitted_bound=sigma*(2-1/m)-1
 assert omitted_bound>np.sqrt(2)
 # Half-grid offset avoids relying on vertices exactly on periodic boundaries.
 ijk=np.indices((m,m,m),dtype=np.int32).reshape(3,-1).T
 points=((ijk+.5)/m)@L
 dist=tree.query(points,workers=-1)[0]
 free=dist>np.sqrt(2)+1e-10
 ids=np.flatnonzero(free);mapping=np.full(len(points),-1,dtype=np.int32);mapping[ids]=np.arange(len(ids))
 sources=[];targets=[];wraps=[];mins=[];uncertain_count=0
 for axis,stride in enumerate([m*m,m,1]):
  endijk=ijk[ids].copy();endijk[:,axis]=(endijk[:,axis]+1)%m
  endid=(endijk[:,0]*m+endijk[:,1])*m+endijk[:,2]
  valid=free[endid];a=ids[valid];b=endid[valid]
  step=L[axis]/m;length=np.linalg.norm(step)
  safe=(np.minimum(dist[a],dist[b])-length/2>np.sqrt(2)+1e-10)
  uncertain=np.flatnonzero(~safe);uncertain_count+=len(uncertain)
  # Exact segment/site projection for every edge not certified by Lipschitz bound.
  for k in uncertain:
   x=points[a[k]];mid=x+step/2
   possible=tree.query_ball_point(mid,np.sqrt(2)+length/2+1e-9)
   if not possible:safe[k]=True;continue
   rel=sites[possible]-x
   tau=np.clip(rel@step/(step@step),0,1)
   minimum=np.min(np.linalg.norm(rel-tau[:,None]*step,axis=1))
   if minimum>np.sqrt(2)+1e-10:safe[k]=True
  a=a[safe];b=b[safe];w=np.zeros((len(a),3),dtype=np.int32);w[:,axis]=(ijk[a,axis]==m-1)
  sources.append(mapping[a]);targets.append(mapping[b]);wraps.append(w)
 U=np.concatenate(sources);V=np.concatenate(targets);W=np.concatenate(wraps)
 # CSR adjacency stores signed edge IDs. BFS assigns integer covering-space cells.
 edgeidx=np.arange(len(U),dtype=np.int64)+1
 graph=coo_matrix((np.r_[edgeidx,-edgeidx],(np.r_[U,V],np.r_[V,U])),shape=(len(ids),len(ids))).tocsr()
 parent=np.full(len(ids),-1,dtype=np.int32)
 seen=np.zeros(len(ids),bool);potential=np.zeros((len(ids),3),dtype=np.int32);labels=np.full(len(ids),-1,dtype=np.int32);sizes=[]
 for root in range(len(ids)):
  if seen[root]:continue
  label=len(sizes);seen[root]=True;labels[root]=label;q=deque([root]);size=0
  while q:
   x=q.popleft();size+=1
   for loc in range(graph.indptr[x],graph.indptr[x+1]):
    y=graph.indices[loc]
    if seen[y]:continue
    eid=int(graph.data[loc]);direction=1 if eid>0 else -1
    parent[y]=x
    potential[y]=potential[x]+direction*W[abs(eid)-1]
    seen[y]=True;labels[y]=label;q.append(y)
  sizes.append(size)
 cycles=potential[U]+W-potential[V]
 nonzero=np.flatnonzero(np.any(cycles!=0,axis=1));components=[]
 for label in sorted(set(labels[U[nonzero]].tolist())):
  mask=labels[U[nonzero]]==label;bs=winding_basis(cycles[nonzero[mask]])
  components.append(dict(vertices=sizes[label],rank=len(bs),basis=bs))
 witnesses=[]
 if m==96:
  # Reconstruct independent cycles and recertify every segment against a larger site set.
  ns3=np.array(list(itertools.product(range(-3,4),repeat=3)))
  sites3=((ns3@L)[:,None,:]+P).reshape(-1,3);tree3=cKDTree(sites3)
  chosen=[]
  for e in nonzero:
   if len(winding_basis(chosen+[cycles[e]]))==len(chosen):continue
   chosen.append(cycles[e])
   def ancestors(node):
    chain=[int(node)]
    while parent[chain[-1]]>=0:chain.append(int(parent[chain[-1]]))
    return chain
   left=ancestors(U[e])[::-1];right=ancestors(V[e]);wind=cycles[e]
   fracroute=[(ijk[ids[n]]+.5)/m+potential[n] for n in left]
   fracroute += [(ijk[ids[n]]+.5)/m+potential[n]+wind for n in right]
   minsep=float('inf')
   for f,g in zip(fracroute,fracroute[1:]):
    segment_start=(f-np.floor(f))@L;step=(g-f)@L;ln=np.linalg.norm(step)
    near=tree3.query_ball_point(segment_start+step/2,np.sqrt(2)+ln/2+.1)
    if not near:
     minsep=min(minsep,np.sqrt(2)+.1);continue
    rel=sites3[near]-segment_start;tau=np.clip(rel@step/(step@step),0,1)
    sep=float(np.min(np.linalg.norm(rel-tau[:,None]*step,axis=1)))
    minsep=min(minsep,sep)
   assert minsep>np.sqrt(2)+1e-10
   assert np.max(np.abs(fracroute[-1]-fracroute[0]-wind))<1e-12
   witnesses.append(dict(winding=wind.tolist(),segments=len(fracroute)-1,minimum_site_distance=minsep,clearance_above_required_distance=float(minsep-np.sqrt(2)),fractional_lifted_route=np.array(fracroute).tolist()))
   if len(chosen)==3:break
 return dict(case=name,angle=degree,resolution=m,determinant=float(np.linalg.det(L)),free_grid_nodes=len(ids),accepted_edges=len(U),exact_edge_checks=uncertain_count,max_winding_rank=max([c['rank'] for c in components],default=0),winding_components=components,witness_routes=witnesses,omitted_image_distance_lower_bound=float(omitted_bound),seconds=round(time.time()-start,2))
if __name__=='__main__':
 rows=[]
 for m in [32,64,96]:
  for name,degree,C in CASES:
   row=run(name,degree,C,m);rows.append(row);print(json.dumps(row),flush=True)
   Path(__file__).with_suffix('.json').write_text(json.dumps(dict(method='Uniform fractional midpoint grid, periodic six-neighbor edges, analytic segment clearance and BFS covering-space cycle potentials. A positive winding is a valid polygonal path. Missing winding on finite grid does not exclude another path.',results=rows),indent=2)+'\n')
