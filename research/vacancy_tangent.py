"""First-order tangent cone at the fixed-motif vacancy endpoint. No finite path claim."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction


def main():
    C=np.array([[0,-1,-1],[1,0,-1],[1,1,0]],dtype=int)
    L=2*np.eye(3,dtype=int)+C
    P=np.vstack([np.eye(3,dtype=int),-np.eye(3,dtype=int)])
    V=np.unique((P[:,None]-P[None,:]).reshape(-1,3),axis=0)
    classes=[]
    # C is skew, hence |nL|^2=4|n|^2+|nC|^2. If |n|>=2,
    # |nL+v|>=2|n|-2>=2>sqrt(2). The 26 neighbors suffice exactly.
    assert np.array_equal(C.T,-C)
    for n in itertools.product([-1,0,1],repeat=3):
        n=np.array(n)
        if not np.any(n): continue
        for v in V:
            d=n@L+v
            if d@d==2:
                classes.append({'n':n.tolist(),'v':v.tolist(),'d':d.tolist(),
                                'row':np.outer(n,d).ravel().tolist()})
    G=np.unique(np.array([c['row'] for c in classes]),axis=0)
    E=np.zeros((3,9),dtype=int)
    for i in range(3):E[i,3*i:3*i+3]=C[i]
    cofactor=np.array([[5,-3,-1],[1,5,-3],[3,1,5]],dtype=int).ravel()
    # Bar constraints parameterize deltaL as [[x,y,-y],[z,w,z],[u,-u,v]].
    B=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,-1,0,0,0,0],
                [0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],
                [0,0,0,0,1,0],[0,0,0,0,-1,0],[0,0,0,0,0,1]],dtype=int)
    assert not np.any(E@B)
    H=np.array([[1,1,0,0,0,0],[1,-1,0,0,0,0],
                [0,0,1,1,0,0],[0,0,-1,1,0,0],
                [0,0,0,0,1,1],[0,0,0,0,-1,1]],dtype=int)
    reduced=G@B
    assert all(any(np.array_equal(h,row) for row in reduced) for h in H)
    sums=[np.array(c)@H for c in itertools.product([0,1],repeat=6)]
    assert all(any(np.array_equal(row,total) for total in sums) for row in reduced)
    assert np.array_equal(cofactor@B,[5,-2,-2,5,2,5])
    nonbar=np.array([r for r in G if not any(np.array_equal(r,e) for e in E)])
    twice_lambda=np.array([3,3,2,3,3,2,3,3,2])
    assert np.array_equal(nonbar.T@twice_lambda+E.sum(axis=0),2*cofactor)
    bounds=[(-1,1)]*9
    result=linprog(cofactor,A_ub=-G,b_ub=np.zeros(len(G)),A_eq=E,b_eq=np.zeros(3),bounds=bounds,method='highs')
    assert result.success
    rational=[Fraction(float(x)).limit_denominator(10000) for x in result.x]
    def exact_dot(row):return sum(int(a)*b for a,b in zip(row,rational))
    assert all(exact_dot(row)>=0 for row in G)
    assert all(exact_dot(row)==0 for row in E)
    exact_derivative=exact_dot(cofactor)
    # A uniform omitted-neighbor bound: outside the box |n|>=4, and
    # |nL+v| >= sigma_min(L)*4-2 > sqrt(2).
    singular_min=float(np.linalg.svd(L,compute_uv=False)[-1])
    assert 4*singular_min-2>np.sqrt(2)
    neutral=[]
    for j in range(9):
        for sign in [-1,1]:
            objective=np.zeros(9);objective[j]=-sign
            res=linprog(objective,A_ub=np.vstack([-G,cofactor]),b_ub=np.zeros(len(G)+1),
                        A_eq=E,b_eq=np.zeros(3),bounds=bounds,method='highs')
            neutral.append({'coordinate':j,'sign':sign,'maximum':float(sign*res.x[j]) if res.success else None})
    output={'endpoint_lattice':L.tolist(),'bridge_vectors':C.tolist(),
            'determinant':14,'determinant_gradient':cofactor.reshape(3,3).tolist(),
            'directed_active_classes':len(classes),'unique_inequality_rows':len(G),
            'box_bounded_min_det_derivative_float':float(result.fun),
            'exact_candidate_det_derivative':str(exact_derivative),
            'exact_candidate_delta_L':[[str(x) for x in rational[i:i+3]] for i in range(0,9,3)],
            'exact_candidate_min_active_pair_dot':str(min(exact_dot(row) for row in G)),
            'exact_candidate_bar_dot_products':[str(exact_dot(row)) for row in E],
            'max_each_coordinate_with_nonpositive_derivative':neutral,
            'exact_reduced_cone':'deltaL=[[x,y,-y],[z,w,z],[u,-u,v]], x>=abs(y), w>=abs(z), v>=abs(u)',
            'exact_positive_bound':'d(detL)=5(x+w+v)-2y-2z+2u >=3(x+w+v)>=3*max(abs(deltaL entries))',
            'exact_nonbar_dual_weights_twice':twice_lambda.tolist(),
            'exact_bar_dual_weights_twice':[1,1,1],
            'exact_dual_identity_verified':True,
            'exact_active_set_completeness':'C skew => |nL|>=2|n|; |v|<=2; |n|>=2 => separation>=2>sqrt(2)',
            'outside_contact_search_lower_bound_float':4*singular_min-2,
            'active_classes':classes,'inequality_rows':G.tolist(),
            'scope':'Necessary first-order cone only, fixed motif, fixed bridge lengths to first order. Feasible tangent does not establish a finite feasible path.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['active_classes','inequality_rows']},indent=2))


if __name__=='__main__':main()
