"""Periodic sphere-void traversal from Voronoi edges and integer winding.
Floating-point geometry; analytic segment clearances; no connector obstacles.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi

P = np.r_[np.eye(3),-np.eye(3)]
R = 1/math.sqrt(2)
CS = {'fcc':np.array([[0,-1,-1],[1,0,-1],[1,-1,0]]),
      'vacancy':np.array([[0,1,1],[-1,0,1],[-1,-1,0]])}


def rank(vectors):
    basis=[]
    for v in vectors:
        v=np.array(v,dtype=np.int64)
        if not np.any(v): continue
        if len(basis)==0: basis.append(v)
        elif len(basis)==1:
            if np.any(np.cross(basis[0],v)): basis.append(v)
        elif len(basis)==2:
            if int(np.dot(np.cross(basis[0],basis[1]),v)): return [*basis,v]
    return basis


class WindingUF:
    def __init__(self,n):
        self.parent=list(range(n)); self.offset=[np.zeros(3,dtype=np.int64) for _ in range(n)]; self.basis=[[] for _ in range(n)]
    def find(self,x):
        if self.parent[x]!=x:
            old=self.parent[x]; root,delta=self.find(old)
            self.offset[x]=self.offset[x]+delta; self.parent[x]=root
        return self.parent[x],self.offset[x].copy()
    def add(self,u,v,shift):
        ru,pu=self.find(u); rv,pv=self.find(v)
        if ru==rv:
            self.basis[ru]=rank(self.basis[ru]+[pv-pu-shift])
        else:
            self.parent[rv]=ru; self.offset[rv]=np.array(shift)+pu-pv
            self.basis[ru]=rank(self.basis[ru]+self.basis[rv])
    def maxrank(self):
        return max((len(self.basis[i]) for i in range(len(self.parent)) if self.parent[i]==i),default=0)


def probe(deg,C,rep=2):
    theta=math.radians(deg)
    L=(2+math.sqrt(2)*math.cos(theta))*np.eye(3)+math.sin(theta)*C
    ns=np.array(list(itertools.product(range(-rep,rep+1),repeat=3)))
    sites=((ns@L)[:,None,:]+P).reshape(-1,3)
    vor=Voronoi(sites)
    verts=vor.vertices; frac=verts@np.linalg.inv(L)
    # Snap only periodic boundary roundoff; never displace interior vertices.
    frac[np.abs(frac-np.rint(frac))<1e-7]=np.rint(frac[np.abs(frac-np.rint(frac))<1e-7])
    shifts=np.floor(frac).astype(int); wrapped=frac-shifts
    keys=[tuple(np.round(f,7)) for f in wrapped]
    vertex_map={}; edges={}; lengths=[]; details={}
    for face in vor.ridge_vertices:
        if len(face)<2: continue
        for a,b in zip(face,face[1:]+face[:1]):
            if a<0 or b<0: continue
            mid=(frac[a]+frac[b])/2
            if not np.all((mid>=-1e-7)&(mid<1-1e-7)): continue
            x,y=verts[a],verts[b]; delta=y-x; norm=float(delta@delta)
            if norm<1e-16: continue
            u=vertex_map.setdefault(keys[a],len(vertex_map)); v=vertex_map.setdefault(keys[b],len(vertex_map))
            shift=shifts[b]-shifts[a]
            forward=(u,v,*shift.tolist()); backward=(v,u,*(-shift).tolist()); ek=min(forward,backward)
            if ek in edges: continue
            alpha=np.clip((sites-x)@delta/norm,0,1)
            residual=sites-x-alpha[:,None]*delta
            ds=np.einsum('ij,ij->i',residual,residual)
            nearest=int(np.argmin(ds)); closest=x+alpha[nearest]*delta
            clearance=float(np.sqrt(ds[nearest])-R)
            edges[ek]=clearance; lengths.append(math.sqrt(norm))
            support=np.linalg.norm(sites-closest,axis=1)
            details[ek]={'segment_parameter':float(alpha[nearest]),'bottleneck_cartesian':closest.tolist(),'nearest_sites':sites[np.abs(support-(clearance+R))<1e-6].tolist()}
    ordered=sorted(edges.items(),key=lambda x:-x[1]); uf=WindingUF(len(vertex_map)); thresholds={}; witness={}
    for _,group in itertools.groupby(ordered,key=lambda x:round(x[1],9)):
        items=list(group)
        for (u,v,*s),w in items: uf.add(u,v,np.array(s,dtype=int))
        r=uf.maxrank()
        for dim in range(1,r+1):
            if dim not in thresholds:
                thresholds[dim]=min(w for _,w in items)
                roots=[i for i in range(len(uf.parent)) if uf.parent[i]==i and len(uf.basis[i])>=dim]
                witness[dim]=[v.tolist() for v in uf.basis[roots[0]]]
    assert 3 in thresholds, ('no winding rank3',deg,rep,len(edges))
    return {'angle':deg,'replication':rep,'vertices':len(vertex_map),'edges':len(edges),'rank_thresholds':{str(k):{'radius':v,'ratio':v/R,'winding_basis':witness[k],'one_bottleneck_edge':details[min(edges,key=lambda e:abs(edges[e]-v))]} for k,v in thresholds.items()},'edge_clearance_range':[min(edges.values()),max(edges.values())],'max_edge_length':max(lengths)}


def matched_density_controls():
    from scipy.optimize import brentq
    def determinant(deg,C):
        theta=math.radians(deg)
        return float(np.linalg.det((2+math.sqrt(2)*math.cos(theta))*np.eye(3)+math.sin(theta)*C))
    controls=[]
    for a in [30,40,45,50,60,75]:
        volume=determinant(a,CS['fcc'])
        b=brentq(lambda deg:determinant(deg,CS['vacancy'])-volume,0,90,xtol=1e-11)
        original=probe(a,CS['fcc'],2); vacancy=probe(b,CS['vacancy'],2)
        difference=abs(determinant(b,CS['vacancy'])-volume)
        assert difference<1e-10
        row={'original_angle':a,'vacancy_angle':b,'cell_volume':volume,'packing_fraction':2*math.sqrt(2)*math.pi/volume,'volume_match_absolute_error':difference,'original':original['rank_thresholds'],'vacancy':vacancy['rank_thresholds']}
        if a==40:
            checks=[]
            for rep in [3,4]:
                for route,deg,base in [('fcc',a,original),('vacancy',b,vacancy)]:
                    test=probe(deg,CS[route],rep)
                    d=max(abs(test['rank_thresholds'][k]['radius']-base['rank_thresholds'][k]['radius']) for k in ['1','2','3'])
                    assert d<1e-7
                    checks.append({'route':route,'replication':rep,'maximum_threshold_radius_difference':d})
            row['replication_checks']=checks
            assert original['rank_thresholds']['3']['ratio']<1<original['rank_thresholds']['2']['ratio']
            assert vacancy['rank_thresholds']['3']['ratio']>1
        controls.append(row)
    return controls


def main():
    results=[]
    for route,C in CS.items():
        for deg in [0,15,30,45,60,75,85,89,90]:
            row=probe(deg,C,2); row['route']=route
            check=probe(deg,C,3)
            diffs={k:abs(row['rank_thresholds'][k]['radius']-check['rank_thresholds'][k]['radius']) for k in row['rank_thresholds']}
            assert max(diffs.values())<1e-7,(route,deg,diffs)
            row['replication3_max_threshold_difference']=max(diffs.values())
            row['replication3_graph_counts']=[check['vertices'],check['edges']]
            results.append(row)
            print(route,deg,{k:round(v['ratio'],8) for k,v in row['rank_thresholds'].items()},flush=True)
    fcc=next(x for x in results if x['route']=='fcc' and x['angle']==90)
    assert abs(fcc['rank_thresholds']['3']['ratio']-(2/math.sqrt(3)-1))<1e-7
    assert abs(results[0]['rank_thresholds']['3']['ratio']-(math.sqrt(4+2*math.sqrt(2))-1))<1e-7
    extra=[]
    for route,C in CS.items():
        for deg in [0,45,89,90]:
            check=probe(deg,C,4)
            baseline=next(x for x in results if x['route']==route and x['angle']==deg)
            difference=max(abs(check['rank_thresholds'][k]['radius']-baseline['rank_thresholds'][k]['radius']) for k in ['1','2','3'])
            assert difference<1e-7
            extra.append({'route':route,'angle':deg,'replication4_max_difference':difference})
    grid=[]
    for deg in range(91):
        a=probe(deg,CS['fcc'],2);b=probe(deg,CS['vacancy'],2)
        x={k:v['ratio'] for k,v in a['rank_thresholds'].items()}; y={k:v['ratio'] for k,v in b['rank_thresholds'].items()}
        assert abs(x['1']-x['2'])<1e-7
        assert max(y.values())-min(y.values())<1e-7
        assert abs(x['2']-y['3'])<1e-7
        grid.append({'angle':deg,'fcc_rank2':x['2'],'fcc_rank3':x['3'],'vacancy_rank3':y['3']})
    from scipy.optimize import brentq
    unit_transitions={}
    for dim in ['2','3']:
        root=brentq(lambda deg:probe(deg,CS['fcc'],2)['rank_thresholds'][dim]['ratio']-1,20,70,xtol=1e-8)
        unit_transitions['rank'+dim+'_loss_angle']=root
    # UF controls: loops without a cell shift cannot imply traversal.
    u=WindingUF(2);u.add(0,1,np.zeros(3,dtype=int));u.add(1,0,np.zeros(3,dtype=int));assert u.maxrank()==0
    u.add(0,1,np.array([1,0,0]));assert u.maxrank()==1
    u.add(0,1,np.array([0,1,0]));assert u.maxrank()==2
    u.add(0,1,np.array([0,0,1]));assert u.maxrank()==3
    Path(__file__).with_suffix('.json').write_text(json.dumps({'method':'Periodic Voronoi edge network; clearance minimized analytically over each full edge against all replicated sphere centers; union-find potentials detect integer winding rank. Numerical geometry, no bars.','rank_definition':'Maximum rank of integer winding subgroup within ONE connected periodic component; directions from disconnected components are never pooled.','results':results,'replication4_checks':extra,'one_degree_grid':grid,'unit_probe_transitions':unit_transitions,'matched_density_controls':matched_density_controls(),'controls':{'fcc_expected_ratio':2/math.sqrt(3)-1,'fcc_pass':True,'open_cubic_expected_ratio':math.sqrt(4+2*math.sqrt(2))-1,'open_cubic_pass':True,'union_find_rank_controls_pass':True}},indent=2)+'\n')


if __name__=='__main__': main()
