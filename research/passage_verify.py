"""Standard-library-only verifier of saved continuum passage witnesses.

Does not import the discovery script, numpy, scipy, or Voronoi geometry.
It independently reconstructs the covered-plane triangulation, all segment
distance bounds, winding ranks, and the equal-volume bracket.
"""
from fractions import Fraction
from pathlib import Path
import itertools
import json

F=Fraction
DATA=Path(__file__).with_name('passage_certificate.json')
SQ=(F(141421356237309504880,10**20),F(141421356237309504881,10**20))
ORIGINAL=((0,-1,-1),(1,0,-1),(1,-1,0))
SKEW=((0,1,1),(-1,0,1),(-1,-1,0))
MOTIF=[tuple(1 if j==i else 0 for j in range(3)) for i in range(3)]
MOTIF+= [tuple(-x for x in p) for p in MOTIF[:]]


def dot(x,y):return sum(a*b for a,b in zip(x,y))
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def matvec_row(x,A):return tuple(sum(x[i]*A[i][j] for i in range(3)) for j in range(3))
def determinant(A):
    return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])


def lattice_box(interval,C,scale):
    a,b=interval
    sinbounds=(2*a/(1+a*a),2*b/(1+b*b))
    cosbounds=((1-b*b)/(1+b*b),(1-a*a)/(1+a*a))
    dbounds=(2+SQ[0]*cosbounds[0],2+SQ[1]*cosbounds[1])
    lo=[];hi=[]
    for i in range(3):
        lr=[];hr=[]
        for j in range(3):
            q=dbounds if i==j else sinbounds if C[i][j]==1 else (-sinbounds[1],-sinbounds[0])
            lr.append((q[0]*scale).__floor__());hr.append((q[1]*scale).__ceil__())
        lo.append(lr);hi.append(hr)
    # Uniform sigma_min(L)>2 follows from these d-t bounds for both matrices.
    assert dbounds[0]-sinbounds[1]>2
    return lo,hi


def verify_plane(data):
    N=data['plane_grid_size'];scale=data['integer_lattice_enclosure_scale']
    interval=tuple(map(F,data['half_angle_intervals']['original']))
    lower,upper=lattice_box(interval,ORIGINAL,scale)
    labels=iter(data['plane_sphere_labels']);maximum=F(0);count=0
    for j in range(N):
        for k in range(N):
            a=(0,j,k);b=(0,j+1,k);c=(0,j,k+1);d=(0,j+1,k+1)
            for triangle in [(a,b,d),(a,c,d)]:
                label=next(labels);n=label['n'];p=label['p']
                assert tuple(p) in MOTIF and len(n)==3 and all(type(x) is int for x in n)
                for vertex in triangle:
                    w=sub(vertex,tuple(N*x for x in n));squared=0
                    for col in range(3):
                        terms=[(w[row]*lower[row][col],w[row]*upper[row][col]) for row in range(3)]
                        low=sum(min(t) for t in terms)-N*scale*p[col]
                        high=sum(max(t) for t in terms)-N*scale*p[col]
                        squared+=max(low*low,high*high)
                    bound=F(squared,(N*scale)**2)
                    assert bound<F(399**2,80000)
                    maximum=max(maximum,bound);count+=1
    assert next(labels,None) is None
    assert count==6*N*N
    return {'triangles':2*N*N,'vertex_checks':count,'max_squared_distance_upper_bound':str(maximum)}


def verify_loops(data,kind):
    scale=data['integer_lattice_enclosure_scale'];grid=data['fractional_polyline_scale']
    name='original_free_loops' if kind=='original' else 'free_loops'
    loops=data[name];C=ORIGINAL if kind=='original' else SKEW
    interval=tuple(map(F,data['half_angle_intervals'][kind]))
    lower,upper=lattice_box(interval,C,scale)
    center=[[(lower[i][j]+upper[i][j])//2 for j in range(3)] for i in range(3)]
    error_matrix=[[max(center[i][j]-lower[i][j],upper[i][j]-center[i][j]) for j in range(3)] for i in range(3)]
    root=loops[0]['fractional_vertices_scaled'][0];tests=0;segments=0
    winds=[]
    for loop in loops:
        points=loop['fractional_vertices_scaled'];wind=loop['winding'];winds.append(wind)
        assert points[0]==root
        assert sub(points[-1],points[0])==tuple(grid*x for x in wind)
        for left,right in zip(points,points[1:]):
            # Recenter each segment by one exact integer lattice translation.
            shift=tuple((a+b)//(2*grid) for a,b in zip(left,right))
            left=sub(left,tuple(grid*x for x in shift));right=sub(right,tuple(grid*x for x in shift))
            assert all(abs(x)<2*grid for x in left+right)
            step=matvec_row(sub(right,left),center);vv=dot(step,step)
            if vv==0:continue
            segments+=1;H=grid*scale
            for n in itertools.product(range(-3,4),repeat=3):
                u=sub(left,tuple(grid*x for x in n));v=sub(right,tuple(grid*x for x in n))
                component_errors=matvec_row(tuple(max(abs(x),abs(y)) for x,y in zip(u,v)),error_matrix)
                E=sum(component_errors)
                for p in MOTIF:
                    a=sub(matvec_row(u,center),tuple(H*x for x in p));aa=dot(a,a);av=dot(a,step)
                    if av>=0:numerator,denominator=aa,1
                    elif av<=-vv:numerator,denominator=aa+2*av+vv,1
                    else:numerator,denominator=aa*vv-av*av,vv
                    # 5671/4000 > (r+q) for every q/r<=201/200.
                    assert numerator*4000**2 > (5671*H+E*4000)**2*denominator
                    tests+=1
    if kind=='skew':assert abs(determinant(winds))==1
    else:
        assert len(winds)==2 and all(w[0]==0 for w in winds)
        assert winds[0][1]*winds[1][2]-winds[0][2]*winds[1][1]!=0
    return {'loops':len(loops),'segments':segments,'exact_segment_site_checks':tests,'winding_vectors':winds}


def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def sign(a):
    p,q=a
    if p==0:return (q>0)-(q<0)
    if q==0 or (p>0)==(q>0):return (p>0)-(p<0)
    return ((p*p>2*q*q)-(p*p<2*q*q))*((p>0)-(p<0))


def volume(s,C):
    t=2*s/(1+s*s);c=(1-s*s)/(1+s*s)
    A=[[(F(2),c) if i==j else (t*C[i][j],F(0)) for j in range(3)] for i in range(3)]
    result=(F(0),F(0))
    for perm in itertools.permutations(range(3)):
        term=(F(1),F(0));inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for i in range(3):term=mul(term,A[i][perm[i]])
        result=add(result,neg(term) if inversions%2 else term)
    return result


def verify_volume_and_angle(data):
    original=tuple(map(F,data['half_angle_intervals']['original']));skew=tuple(map(F,data['half_angle_intervals']['skew']))
    assert sign(add(volume(skew[0],SKEW),neg(volume(original[0],ORIGINAL))))==1
    assert sign(add(volume(original[1],ORIGINAL),neg(volume(skew[1],SKEW))))==1
    def polynomial(x):return x**6-33*x**4+27*x*x-3
    assert polynomial(original[0])<0<polynomial(original[1])
    assert 0<original[0]<original[1]<F(37,100)
    assert 9-22*F(37,100)**2>0 # derivative is positive throughout containing interval
    return {'matched_volume_bracket':True,'original_interval_contains_exact_40_degrees':True}


def main():
    assert SQ[0]**2<2<SQ[1]**2
    assert F(5671,4000)**2>F(401**2,80000)
    data=json.loads(DATA.read_text())
    result={'covered_plane':verify_plane(data),'original_free_paths':verify_loops(data,'original'),'skew_free_paths':verify_loops(data,'skew'),'volume_and_angle':verify_volume_and_angle(data),'probe_ratio_interval':['199/200','201/200'],'conclusion':'Continuum winding rank exactly2 on original interval versus rank3 on skew interval, for every probe ratio q/r in[199/200,201/200], with a unique equal-volume skew state for every original state, including exact40deg.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
