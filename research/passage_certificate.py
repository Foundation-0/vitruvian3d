"""Continuum passage certificate: covered plane versus explicit free loops.

Discovery uses scipy/numpy; acceptance uses integer/rational inequalities.
All created output is task-owned passage_* content.
"""
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi,cKDTree

BASE=Path(__file__).with_name('passage_certificate')
C_ORIGINAL=np.array([[0,-1,-1],[1,0,-1],[1,-1,0]],dtype=int)
C_SKEW=np.array([[0,1,1],[-1,0,1],[-1,-1,0]],dtype=int)
P=np.r_[np.eye(3,dtype=int),-np.eye(3,dtype=int)]
SQLO=F(141421356237309504880,10**20)
SQHI=F(141421356237309504881,10**20)
ORIGINAL_INTERVAL=(F(363969,10**6),F(363972,10**6))
SKEW_INTERVAL=(F(425610,10**6),F(425635,10**6))
D=10**6
FGRID=10**6


def floor(x): return x.numerator//x.denominator
def ceil(x): return -floor(-x)


def lattice(s,C):
    s=float(s);return (2+math.sqrt(2)*(1-s*s)/(1+s*s))*np.eye(3)+(2*s/(1+s*s))*C


def bounds(interval,C):
    lo,hi=interval
    assert 0<=lo<=hi<1
    # sin(2atan s) increases and cos(2atan s) decreases on this interval.
    tlo,thi=2*lo/(1+lo*lo),2*hi/(1+hi*hi)
    clo,chi=(1-hi*hi)/(1+hi*hi),(1-lo*lo)/(1+lo*lo)
    dlo,dhi=2+SQLO*clo,2+SQHI*chi
    lower=np.zeros((3,3),dtype=np.int64);upper=lower.copy()
    for i in range(3):
        for j in range(3):
            a,b=(dlo,dhi) if i==j else ((tlo,thi) if C[i,j]==1 else (-thi,-tlo))
            lower[i,j]=floor(D*a);upper[i,j]=ceil(D*b)
    return lower,upper


def sign_pair(p,q):
    if not p:return (q>0)-(q<0)
    if not q or (p>0)==(q>0):return (p>0)-(p<0)
    delta=p*p-2*q*q
    return ((delta>0)-(delta<0))*((p>0)-(p<0))


def determinant(s,route):
    t=2*s/(1+s*s);c=(1-s*s)/(1+s*s)
    if route=='original':return (8+12*c*c+2*t*t+2*t**3,12*c+2*c**3+c*t*t)
    return (8+12*c*c+6*t*t,12*c+2*c**3+3*c*t*t)


def volume_bracket():
    ohi=determinant(ORIGINAL_INTERVAL[0],'original')
    olo=determinant(ORIGINAL_INTERVAL[1],'original')
    shi=determinant(SKEW_INTERVAL[0],'skew')
    slo=determinant(SKEW_INTERVAL[1],'skew')
    a=sign_pair(shi[0]-ohi[0],shi[1]-ohi[1]);b=sign_pair(olo[0]-slo[0],olo[1]-slo[1])
    assert a==b==1
    return {'skew_low_s_volume_above_all_original_volumes':True,'skew_high_s_volume_below_all_original_volumes':True,'scope':'For every original half-angle parameter in its interval there exists a unique equal-volume skew parameter in its interval, by continuity and strict determinant decrease.'}


def plane_discovery(N=32):
    s=sum(ORIGINAL_INTERVAL)/2;L=lattice(s,C_ORIGINAL)
    ns=np.array(list(itertools.product(range(-2,3),repeat=3)),dtype=int)
    sites=((ns@L)[:,None,:]+P).reshape(-1,3)
    triangles=[]
    for j in range(N):
        for k in range(N):
            a=[0,j,k];b=[0,j+1,k];c=[0,j,k+1];d=[0,j+1,k+1]
            triangles.extend([[a,b,d],[a,c,d]])
    triangles=np.array(triangles,dtype=np.int64)
    selected=cKDTree(sites).query(triangles.mean(1)@L/N)[1]
    labels=[{'n':ns[i//6].tolist(),'p':P[i%6].tolist()} for i in selected]
    return N,triangles,labels


def plane_check(N,triangles,labels,limit_squared=F(2)):
    lower,upper=bounds(ORIGINAL_INTERVAL,C_ORIGINAL)
    ns=np.array([x['n'] for x in labels],dtype=np.int64)
    ps=np.array([x['p'] for x in labels],dtype=np.int64)
    W=triangles-N*ns[:,None,:]
    lows=np.zeros(W.shape,dtype=np.int64);highs=lows.copy()
    for i in range(3):
        for j in range(3):
            a=W[:,:,i]*lower[i,j];b=W[:,:,i]*upper[i,j]
            lows[:,:,j]+=np.minimum(a,b);highs[:,:,j]+=np.maximum(a,b)
    lows-=N*D*ps[:,None,:];highs-=N*D*ps[:,None,:]
    absolute=np.maximum(np.abs(lows),np.abs(highs))
    assert int(absolute.max())<10**9 # excludes int64 square overflow
    sq=np.sum(absolute*absolute,axis=2)
    maximum=int(sq.max());den=(N*D)**2
    accepted=maximum*limit_squared.denominator<limit_squared.numerator*den
    return {'accepted':bool(accepted),'triangle_count':len(labels),'vertex_tests':int(sq.size),'largest_distance_squared_upper_bound_fraction':[maximum,den],'strict_squared_margin_fraction':[2*den-maximum,den],'distance_upper_bound_float':math.sqrt(maximum/den),'plane':'fractional coordinate f1=0, repeated by L2,L3; all integer f1 planes are covered'}


def periodic_graph(s,C):
    L=lattice(s,C);ns=np.array(list(itertools.product(range(-2,3),repeat=3)))
    sites=((ns@L)[:,None,:]+P).reshape(-1,3);vor=Voronoi(sites)
    frac=vor.vertices@np.linalg.inv(L)
    near=np.abs(frac-np.rint(frac))<1e-7;frac[near]=np.rint(frac[near])
    shifts=np.floor(frac).astype(int);wrap=frac-shifts
    keys=[tuple(np.rint(f*FGRID).astype(int)) for f in wrap]
    nodes={};coords=[];edges={}
    for face in vor.ridge_vertices:
        for a,b in zip(face,face[1:]+face[:1]):
            if min(a,b)<0:continue
            mid=(frac[a]+frac[b])/2
            if not np.all((mid>=-1e-7)&(mid<1-1e-7)):continue
            x,y=vor.vertices[a],vor.vertices[b];v=y-x;vv=float(v@v)
            if vv<1e-16:continue
            u=nodes.setdefault(keys[a],len(nodes));z=nodes.setdefault(keys[b],len(nodes))
            shift=shifts[b]-shifts[a]
            key=min((u,z,*shift.tolist()),(z,u,*(-shift).tolist()))
            if key in edges:continue
            alpha=np.clip((sites-x)@v/vv,0,1);res=sites-x-alpha[:,None]*v
            clearance=float(np.sqrt(np.min(np.sum(res*res,axis=1))))
            if clearance>math.sqrt(2)+.002:edges[key]=clearance
    coords=[None]*len(nodes)
    for key,i in nodes.items():coords[i]=np.array(key,dtype=np.int64)
    return coords,edges


def rank(vectors):
    basis=[]
    for v in vectors:
        if not np.any(v):continue
        if not basis or (len(basis)==1 and np.any(np.cross(basis[0],v))) or (len(basis)==2 and np.dot(np.cross(basis[0],basis[1]),v)!=0):basis.append(v)
        if len(basis)==3:break
    return len(basis)


def loop_discovery(C=C_SKEW,interval=SKEW_INTERVAL,target_rank=3):
    coords,edges=periodic_graph(sum(interval)/2,C)
    adj=[[] for _ in coords]
    for u,v,*shift in edges:
        shift=np.array(shift,dtype=np.int64)
        adj[u].append((v,shift));adj[v].append((u,-shift))
    visited=set()
    for root in range(len(coords)):
        if root in visited:continue
        offsets={root:np.zeros(3,dtype=np.int64)};parent={root:None};queue=[root];cycles=[];vectors=[]
        for u in queue:
            visited.add(u)
            for v,shift in adj[u]:
                proposed=offsets[u]+shift
                if v not in offsets:
                    offsets[v]=proposed;parent[v]=u;queue.append(v)
                else:
                    wind=proposed-offsets[v]
                    if rank(vectors+[wind])>len(vectors):
                        cycles.append((u,v,wind));vectors.append(wind)
            if len(vectors)==target_rank:break
        if len(vectors)<target_rank:continue
        loops=[]
        def ancestry(node):
            result=[]
            while node is not None:result.append(node);node=parent[node]
            return result
        for u,v,wind in cycles:
            points=[coords[x]+FGRID*offsets[x] for x in reversed(ancestry(u))]
            points += [coords[x]+FGRID*(offsets[x]+wind) for x in ancestry(v)]
            loops.append({'winding':wind.tolist(),'fractional_vertices_scaled':[p.tolist() for p in points]})
        assert rank([np.array(x['winding']) for x in loops])==target_rank
        return loops
    raise AssertionError('No rank3 component in discovery graph')


def loops_check(loops,radius_upper=F(1414213563,10**9),C=C_SKEW,interval=SKEW_INTERVAL):
    lower,upper=bounds(interval,C)
    # Skew C has symmetric part zero; original C has symmetric-part min eigenvalue -1.
    assert int(lower[0,0])-max(int(abs(lower[0,1])),int(abs(upper[0,1])))>2*D
    center=(lower+upper)//2;delta=np.maximum(center-lower,upper-center)
    ns=list(itertools.product(range(-3,4),repeat=3))
    tests=0;minimum_margin=math.inf;segment_count=0
    for loop in loops:
        points=np.array(loop['fractional_vertices_scaled'],dtype=np.int64)
        assert np.array_equal(points[-1]-points[0],FGRID*np.array(loop['winding']))
        for first,second in zip(points,points[1:]):
            shift=np.floor((first+second)/(2*FGRID)).astype(np.int64)
            first=first-FGRID*shift;second=second-FGRID*shift
            assert max(np.abs(first).max(),np.abs(second).max())<2*FGRID
            step=(second-first)@center
            vv=sum(int(x)**2 for x in step)
            if vv==0:continue
            segment_count+=1
            H=FGRID*D
            for n in ns:
                w0=first-FGRID*np.array(n);w1=second-FGRID*np.array(n)
                maxw=np.maximum(np.abs(w0),np.abs(w1))
                # Coordinate error bounds summed give an exact upper bound on norm error.
                error=int(np.sum(maxw@delta))
                for p in P:
                    a=w0@center-FGRID*D*p
                    aa=sum(int(x)**2 for x in a);av=sum(int(x)*int(y) for x,y in zip(a,step))
                    if av>=0:num,den=aa,1
                    elif av<=-vv:num,den=aa+2*av+vv,1
                    else:num,den=aa*vv-av*av,vv
                    target=radius_upper.numerator*H+error*radius_upper.denominator
                    lhs=num*radius_upper.denominator**2;rhs=target*target*den
                    tests+=1
                    if lhs<=rhs:return {'accepted':False,'failed_loop_winding':loop['winding'],'failed_image':list(n),'failed_site':p.tolist(),'segments_checked':segment_count,'distance_tests':tests}
                    margin=math.sqrt(num/den)/H-error/H-float(radius_upper)
                    minimum_margin=min(minimum_margin,margin)
    return {'accepted':True,'loop_count':len(loops),'winding_vectors':[x['winding'] for x in loops],'segment_count':segment_count,'distance_tests':tests,'minimum_certified_distance_margin_float':minimum_margin,'scope':'All real half-angle parameters throughout selected interval; all segment parameters; every periodic site. Distant images excluded analytically using sigma_min(L)>2 and |fractional segment point|_infinity<2.'}


def main():
    assert SQLO*SQLO<2<SQHI*SQHI
    N,triangles,labels=plane_discovery();plane=plane_check(N,triangles,labels);assert plane['accepted']
    loops=loop_discovery();free=loops_check(loops);assert free['accepted']
    original_loops=loop_discovery(C_ORIGINAL,ORIGINAL_INTERVAL,2)
    original_free=loops_check(original_loops,C=C_ORIGINAL,interval=ORIGINAL_INTERVAL);assert original_free['accepted']
    # A nonzero interval of probe sizes, not just one specially tuned radius.
    small_probe_plane=plane_check(N,triangles,labels,F(399**2,80000));assert small_probe_plane['accepted']
    assert F(5671,4000)**2>F(401**2,80000)
    large_probe_skew=loops_check(loops,F(5671,4000));assert large_probe_skew['accepted']
    large_probe_original=loops_check(original_loops,F(5671,4000),C_ORIGINAL,ORIGINAL_INTERVAL);assert large_probe_original['accepted']
    bracket=volume_bracket()
    def tan20poly(s):return s**6-33*s**4+27*s*s-3
    assert tan20poly(ORIGINAL_INTERVAL[0])<0<tan20poly(ORIGINAL_INTERVAL[1])
    # A plane-cover control below the actual maximum cannot pass; sufficiently enlarged obstacles block the proposed loops.
    plane_control=plane_check(N,triangles,labels,F(193,100));assert not plane_control['accepted']
    free_control=loops_check(loops,F(143,100));assert not free_control['accepted']
    output={'model':'Sphere centers +/-e_i+nL, obstacle radius1/sqrt(2), probe radius1/sqrt(2); excluded center distance sqrt(2).','half_angle_intervals':{'original':[str(x) for x in ORIGINAL_INTERVAL],'skew':[str(x) for x in SKEW_INTERVAL]},'integer_lattice_enclosure_scale':D,'fractional_polyline_scale':FGRID,'plane_certificate':plane,'free_loop_certificate':free,'original_free_loop_certificate':original_free,'probe_ratio_interval':['199/200','201/200'],'probe_size_window_certificates':{'smallest_probe_plane':small_probe_plane,'largest_probe_skew_free':large_probe_skew,'largest_probe_original_free':large_probe_original},'matched_volume_bracket':bracket,'controls':{'insufficient_plane_cover_radius_squared_1_93_rejected':True,'loop_obstacle_radius_1_43_rejected':True,'original_interval_contains_tan20_by_polynomial_sign':True},'plane_grid_size':N,'plane_sphere_labels':labels,'free_loops':loops,'original_free_loops':original_loops}
    BASE.with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['plane_sphere_labels','free_loops','original_free_loops']},indent=2))


if __name__=='__main__':main()
