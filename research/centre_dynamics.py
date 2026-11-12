"""Bounded experiments on face-pyramid Chebyshev-centre iteration.

Actual coplanar hull triangles are merged into polygonal faces. Origin is the
solid volume centroid; radius normalization is RMS over hull vertices. Largest
ball ties choose the unique minimum-Euclidean-norm centre. No novelty claim.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull, distance
from scipy.optimize import linprog
from centre_extraction import tetra_incenter


def faces(points, tolerance):
    hull = ConvexHull(points)
    groups = []
    for tri, equation in zip(hull.simplices, hull.equations):
        for group in groups:
            if np.max(np.abs(equation - group['equation'])) < tolerance:
                group['vertices'].update(tri.tolist())
                break
        else:
            groups.append({'equation': equation, 'vertices': set(tri.tolist())})
    return hull, groups


def normalize(points):
    hull = ConvexHull(points)
    points = points[hull.vertices]
    hull = ConvexHull(points)
    reference = points.mean(axis=0)
    centroid_sum, volume_sum = np.zeros(3), 0.
    for tri in hull.simplices:
        vertices = points[tri]
        volume = abs(np.linalg.det(vertices-reference)) / 6
        centroid_sum += volume * (vertices.sum(axis=0)+reference) / 4
        volume_sum += volume
    centroid = centroid_sum / volume_sum
    points = points-centroid
    return points / np.sqrt(np.mean(np.sum(points**2, axis=1)))


def minimum_norm_center(A, b, r, tolerance=2e-8):
    rhs = b-r
    candidates = []
    if np.max(-rhs) < tolerance:
        candidates.append(np.zeros(3))
    for size in [1, 2, 3]:
        for subset in itertools.combinations(range(len(A)), size):
            sub = A[list(subset)]
            if np.linalg.matrix_rank(sub, tol=1e-10) < size:
                continue
            candidate = np.linalg.lstsq(sub, rhs[list(subset)], rcond=None)[0]
            if np.max(A @ candidate-rhs) < tolerance:
                candidates.append(candidate)
    if not candidates:
        raise RuntimeError('No minimum-norm optimal centre candidate')
    return min(candidates, key=lambda x: np.dot(x,x))


def center(pyramid, tolerance):
    if len(pyramid) == 4:
        return tetra_incenter(pyramid)
    _, groups = faces(pyramid, tolerance)
    A = np.array([g['equation'][:3] for g in groups])
    b = -np.array([g['equation'][3] for g in groups])
    result = linprog([0,0,0,-1], A_ub=np.column_stack([A,np.ones(len(A))]),
                     b_ub=b, bounds=[(None,None)]*3+[(0,None)], method='highs',
                     options={'primal_feasibility_tolerance':1e-9,
                              'dual_feasibility_tolerance':1e-9})
    if not result.success:
        raise RuntimeError(result.message)
    return minimum_norm_center(A,b,result.x[3])


def step(points, tolerance):
    _, groups = faces(points,tolerance)
    return normalize(np.array([center(np.vstack([np.zeros(3),points[sorted(g['vertices'])]]),tolerance)
                               for g in groups]))


def hausdorff(a,b):
    d = distance.cdist(a,b)
    return float(max(d.min(axis=0).max(), d.min(axis=1).max()))


def support_gap(a,b):
    directions=np.vstack([a,b,np.array(list(itertools.product([-1.,1.],repeat=3)))])
    directions=directions/np.linalg.norm(directions,axis=1)[:,None]
    return float(np.max(np.abs(np.max(a @ directions.T,axis=0)-np.max(b @ directions.T,axis=0))))


def diagnostics(points, tolerance):
    hull, groups = faces(points,tolerance)
    eig = np.linalg.eigvalsh(points.T @ points/len(points))
    radii=np.linalg.norm(points,axis=1)
    return {'vertices':len(points),'faces':len(groups),
            'face_sizes':sorted(len(g['vertices']) for g in groups),
            'covariance_eigenvalues':eig.tolist(), 'vertex_radius_cv':float(radii.std()/radii.mean()),
            'volume':float(hull.volume)}


def trajectory(label, initial, tolerance, count=10, cap=260):
    states=[normalize(initial)]
    rows=[]
    for index in range(count+1):
        points=states[-1]
        row={'iteration':index,**diagnostics(points,tolerance)}
        if len(states)>=3:
            row['two_step_hausdorff_fixed_frame']=hausdorff(points,states[-3])
        rows.append(row)
        if index==count or row['faces']>cap:
            break
        states.append(step(points,tolerance))
    return {'label':label,'coplanar_tolerance':tolerance,'trajectory':rows,
            'initial_points':states[0].tolist(),'final_points':states[-1].tolist(),
            'stopped_at_complexity_cap':rows[-1]['faces']>cap}


def main():
    signs=np.array(list(itertools.product([-1.,1.],repeat=3)))
    octa=np.vstack([np.eye(3),-np.eye(3)])
    rng=np.random.default_rng(20261112)
    noise=rng.normal(size=signs.shape)
    cases=[('cube',signs),('unequal_axis_octa',octa @ np.diag([2,1,3])),
           ('shear_octa_small',octa @ np.array([[1,0,0],[.02,1,0],[0,0,1]])),
           ('shear_octa_large',octa @ np.array([[2,0,0],[.4,1,0],[0,0,3]])),
           ('sheared_cube',signs @ np.array([[1,0,0],[.1,1,0],[0,0,1]])),
           ('random_octa',octa+.01*rng.normal(size=octa.shape)),
           ('random_cube',signs+.01*noise)]
    results=[trajectory(label,p,1e-8) for label,p in cases]
    robustness=[trajectory('shear_octa_small',cases[2][1],tol,6)
                for tol in [1e-6,1e-10]]
    reference=step(normalize(signs),1e-10)
    discontinuity=[]
    for eps in [1e-2,1e-3,1e-4,1e-5,1e-6]:
        p=normalize(signs+eps*noise)
        q=step(p,1e-10)
        discontinuity.append({'perturbation':eps,'input_hausdorff':hausdorff(p,normalize(signs)),
                              'output_hausdorff':hausdorff(q,reference),
                              'convex_body_support_gap_lower_bound':support_gap(q,reference),
                              'output_vertices':len(q)})
    parity=[]
    for eps in [.49,.1,1e-2,1e-4,1e-6]:
        p=normalize(signs*(1+eps*np.prod(signs,axis=1))[:,None])
        q=step(p,1e-10)
        parity.append({'perturbation':eps,'input_vertex_hausdorff':hausdorff(p,normalize(signs)),
                       'output_vertex_hausdorff':hausdorff(q,reference),
                       'convex_body_support_gap_lower_bound':support_gap(q,reference),
                       'output_vertices':len(q),'points':q.tolist()})
    tetra_cycle=trajectory('tetrahedral_parity_perturbed_cube',signs*(1+.01*np.prod(signs,axis=1))[:,None],1e-9,8)
    parity_q=np.array(parity[-1]['points'])
    parity_q2=step(parity_q,1e-9)
    radii=np.linalg.norm(parity_q2,axis=1)
    expected_ratio=(np.sqrt(6)-1)/(np.sqrt(6)+1)
    parity_q5=parity_q.copy()
    for _ in range(4):
        parity_q5=step(parity_q5,1e-9)
    result={'definition':'Volume-centroid origin; true coplanar faces merged; largest-ball centre with minimum-norm tie-break; hull-vertex RMS normalization.',
            'numerics':'SciPy ConvexHull + analytic tetra incentres; otherwise HiGHS LP and active-set minimum-norm tie-break. Floating-point evidence, no general convergence claim.',
            'runs':results,'tolerance_checks':robustness,'cube_discontinuity_probe':discontinuity,
            'parity_symmetric_discontinuity_probe':parity,
            'tetrahedral_cycle':tetra_cycle,
            'tetrahedral_cycle_checks':{'observed_two_radius_ratio':float(radii.min()/radii.max()),
                                        'exact_ratio':float(expected_ratio),
                                        'four_step_vertex_hausdorff':hausdorff(parity_q,parity_q5)}}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'runs':[{ 'label':r['label'],'counts':[(x['vertices'],x['faces']) for x in r['trajectory']],
                              'last_two_step':r['trajectory'][-1].get('two_step_hausdorff_fixed_frame')}
                             for r in results], 'cube_discontinuity_probe':discontinuity},indent=2))


if __name__=='__main__':
    main()
