"""Task-owned audit of all 64 signed off-diagonal coordination matrices.

Uses numpy for discovery and Python integer arithmetic for certificates.
No change is made to the original motion probe.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np

OUT = Path(__file__).with_suffix('.json')
POS = [(i, j) for i in range(3) for j in range(3) if i != j]
MOTIF = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=int)
DIFF = np.array(sorted(set(tuple(p-q) for p in MOTIF for q in MOTIF)), dtype=int)
IMAGES = np.array([n for n in itertools.product(range(-5,6), repeat=3) if sum(x*x for x in n) <= 33], dtype=int)
ROTATIONS = []
for perm in itertools.permutations(range(3)):
    for signs in itertools.product([-1,1], repeat=3):
        P = np.zeros((3,3), dtype=int)
        for i,j in enumerate(perm):
            P[i,j] = signs[i]
        ROTATIONS.append(P)


def key(C):
    return ''.join('+' if C[i,j] > 0 else '-' for i,j in POS)


def invariants(C):
    pairs = [int(C[i,j]*C[j,i]) for i,j in [(0,1),(0,2),(1,2)]]
    cycles = [int(C[0,1]*C[1,2]*C[2,0]),int(C[0,2]*C[2,1]*C[1,0])]
    return {'symmetric_pair_count':pairs.count(1),'reciprocal_products':pairs,'directed_triangle_products':cycles,'det_C':sum(cycles),'sum_reciprocal_products':sum(pairs)}


def det3(A):
    a,b,c = A.tolist()
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


MATRICES = {}
for signs in itertools.product([-1,1], repeat=6):
    C = np.zeros((3,3), dtype=int)
    for ij,x in zip(POS, signs):
        C[ij] = x
    MATRICES[key(C)] = C


def pair_sign(p, q):
    """Exact sign of p+q sqrt(2), for integer p,q."""
    if p == 0:
        return (q > 0)-(q < 0)
    if q == 0 or (p > 0) == (q > 0):
        return (p > 0)-(p < 0)
    d = p*p-2*q*q
    return ((d > 0)-(d < 0)) * ((p > 0)-(p < 0))


def exact_bernstein(C):
    """Every distance quartic: integer pairs for twelve times Bernstein coeffs."""
    negative = []
    count = 0
    for n0 in IMAGES.tolist():
        n = np.array(n0, dtype=int)
        z = n @ C
        for v0 in DIFF.tolist():
            if n0 == [0,0,0] and v0 == [0,0,0]:
                continue
            h = 2*n+np.array(v0, dtype=int)
            hh, nn, zz = int(h@h), int(n@n), int(z@z)
            hn, hz, nz = int(h@n), int(h@z), int(n@z)
            coeffs = [(hh+2*nn-2,2*hn), (4*hz,4*nz), (4*zz+2*hh-4*nn-4,0), (4*hz,-4*nz), (hh+2*nn-2,-2*hn)]
            for k in range(5):
                weights = [12,3*k,2*math.comb(k,2),3*math.comb(k,3),12*math.comb(k,4)]
                p = sum(w*c[0] for w,c in zip(weights,coeffs))
                q = sum(w*c[1] for w,c in zip(weights,coeffs))
                count += 1
                if pair_sign(p,q) < 0:
                    negative.append({'n':n0,'v':v0,'coefficient':k,'pair_times_12':[p,q]})
                    return {'all_nonnegative':False,'coefficients_checked':count,'first_negative':negative[0]}
    return {'all_nonnegative':True,'coefficients_checked':count,'distance_conditions':count//5}


def numerical_probe(C):
    best = None
    # Each n,v represents a pair of motif sites in periodic cells.
    for deg in range(0,91):
        theta = math.radians(deg)
        L = (2+math.sqrt(2)*math.cos(theta))*np.eye(3)+math.sin(theta)*C
        distances = IMAGES @ L
        delta = distances[:,None,:]+DIFF[None,:,:]
        sq = np.einsum('ijk,ijk->ij',delta,delta)
        sq[np.all(IMAGES==0,axis=1)[:,None] & np.all(DIFF==0,axis=1)[None,:]] = np.inf
        idx = np.unravel_index(np.argmin(sq),sq.shape)
        minimum = float(sq[idx])
        if best is None or minimum < best['distance_squared']:
            best = {'angle_degrees':deg,'distance_squared':minimum,'n':IMAGES[idx[0]].tolist(),'v':DIFF[idx[1]].tolist()}
    return best


def endpoint_contacts(C):
    L = 2*np.eye(3,dtype=int)+C
    contacts = []
    collision = None
    for p in MOTIF:
        degree = 0
        for n in IMAGES:
            for q in MOTIF:
                if np.all(n==0) and np.all(p==q):
                    continue
                v = n@L+q-p
                d = int(v@v)
                if d == 2:
                    degree += 1
                if d < 2 and collision is None:
                    collision = {'n':n.tolist(),'p':p.tolist(),'q':q.tolist(),'distance_squared':d}
        contacts.append(degree)
    return {'determinant':det3(L),'contact_degrees':contacts,'collision':collision}


def main():
    unseen = set(MATRICES)
    orbits = []
    while unseen:
        rep = min(unseen)
        C = MATRICES[rep]
        orbit = sorted(set(key(P@C@P.T) for P in ROTATIONS))
        proper = sorted(set(key(P@C@P.T) for P in ROTATIONS if det3(P) == 1))
        unseen.difference_update(orbit)
        eig = np.linalg.eigvalsh((C+C.T)/2)
        inv = invariants(C)
        result = {'representative':rep,'matrix':C.tolist(),'members':orbit,'orbit_size':len(orbit),'proper_rotation_orbit_size':len(proper),'proper_and_full_orbits_equal':proper==orbit,'invariants':inv,'symmetric_part_eigenvalues':eig.tolist(),'endpoint':endpoint_contacts(C),'sampled_minimum':numerical_probe(C)}
        # For these matrices min eig(S) >= -sqrt(2), so sigma_min L >= 2-sqrt(2).
        # Excluded class has eigenvalue -2 and a singular endpoint.
        if not (inv['symmetric_pair_count'] == 3 and inv['det_C'] == -2):
            result['uniform_singular_value_lower_bound'] = '2-sqrt(2)'
            result['bernstein'] = exact_bernstein(C)
        else:
            result['bernstein'] = {'all_nonnegative':False,'reason':'singular endpoint with exact overlap witness'}
        orbits.append(result)
    original = np.array([[0,-1,-1],[1,0,-1],[1,-1,0]])
    flipped = original.copy(); flipped[0,1] *= -1
    lookup = {m:o for o in orbits for m in o['members']}
    assert lookup[key(original)]['bernstein']['all_nonnegative']
    assert lookup[key(flipped)]['endpoint']['collision']['distance_squared'] == 0
    assert sum(o['orbit_size'] for o in orbits) == 64
    assert sum(o['orbit_size'] for o in orbits if o['bernstein']['all_nonnegative']) == 20
    assert all(o['proper_and_full_orbits_equal'] for o in orbits)
    # Independently compare the algebraic quartic to direct geometric distances.
    rng = np.random.default_rng(5729)
    residual = 0.0
    for C in MATRICES.values():
        for _ in range(20):
            n = IMAGES[int(rng.integers(len(IMAGES)))]; v = DIFF[int(rng.integers(len(DIFF)))]; s = float(rng.uniform())
            a0 = (2+math.sqrt(2))*n+v; a1 = 2*n@C; a2 = (2-math.sqrt(2))*n+v
            f = np.dot(a0+a1*s+a2*s*s,a0+a1*s+a2*s*s)-2*(1+s*s)**2
            theta = 2*math.atan(s)
            L = (2+math.sqrt(2)*math.cos(theta))*np.eye(3)+math.sin(theta)*C
            direct = (1+s*s)**2*(np.dot(n@L+v,n@L+v)-2)
            residual = max(residual,abs(float(f-direct)))
    assert residual < 1e-10
    skew = np.array([[0,1,1],[-1,0,1],[-1,-1,0]])
    skewL = 2*np.eye(3,dtype=int)+skew
    residue = np.array([5,1,3]); missing = np.array([1,-1,1])
    assert det3(skewL) == 14
    assert np.all((skewL@residue)%7 == 0)
    assert set((MOTIF@residue%7).tolist()) == set(range(1,7))
    assert int(missing.sum())%2 == 1 and int(missing@residue)%7 == 0
    data = {'model':'L=(2+sqrt(2) cos(theta)) I+sin(theta) C; motif=+-unit axes; radius=1/sqrt(2)', 'finite_image_bound':'For every nonsingular-endpoint orbit symmetric-part min eigenvalue >= -sqrt(2). Thus |n L+v| >= (2-sqrt(2))|n|-2. If |n|^2 >=34, this exceeds sqrt(2), since (2+sqrt(2))/(2-sqrt(2))=3+2sqrt(2), whose square=17+12sqrt(2)<34. Enumerate |n|^2<=33.','image_count':len(IMAGES),'motif_difference_count':len(DIFF),'orbits':orbits,'controls':{'original_key':key(original),'original_accepted':True,'one_sign_flip_key':key(flipped),'one_sign_flip_rejected_by_exact_endpoint_overlap':True,'quartic_geometry_comparisons':1280,'max_float_residual':residual,'orbit_partition_count':64},'skew_endpoint':{'representative_C':skew.tolist(),'residue_vector_mod7':residue.tolist(),'missing_odd_fcc_site':missing.tolist(),'occupied_residues':[1,2,3,4,5,6],'vacancy_residue':0,'fcc_occupancy_fraction':'6/7','packing_density':'sqrt(2)*pi/7'}}
    OUT.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'orbits':[{'key':x['representative'],'size':x['orbit_size'],'det':x['endpoint']['determinant'],'minimum':x['sampled_minimum']['distance_squared'],'bernstein':x['bernstein']['all_nonnegative'],'proper_equals_full':x['proper_and_full_orbits_equal']} for x in orbits]},indent=2))


if __name__ == '__main__':
    main()
