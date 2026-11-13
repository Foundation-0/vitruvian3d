"""Independent rational checks for the transport note; no author-code imports."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
FOLLOWUP = HERE.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root', type=Path, required=True)
parser.add_argument('--output', type=Path, default=HERE / 'independent-bounds.json')
args = parser.parse_args()
SOURCE = args.source_root.resolve()
PIN = 'e7445d56b9ecfd97b408653a2726f3bab70a6e4c'
FROZEN = {
    'THEOREM.md': '9cf373dd9c9c5451fe76b38a04258abfa6f1ac978aade05e6c3250b8a4a28d7d',
    'check.py': '6f7e8bec76180aec8f651412ec212a2d1c22678881aa9967c028bf36474f2ba8',
    'certificate.json': '29c20b3be2b6b2764764090e080e585c17f1e30ea1c6aa2a142d0d1c80744ac4',
}
for name, digest in FROZEN.items():
    assert hashlib.sha256((FOLLOWUP/name).read_bytes()).hexdigest() == digest
claimed = json.loads((FOLLOWUP/'certificate.json').read_text())
for name, digest in claimed['source_sha256'].items():
    assert hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() == digest
    pinned_bytes = subprocess.check_output(['git', 'show', f'{PIN}:{name}'], cwd=SOURCE)
    assert hashlib.sha256(pinned_bytes).hexdigest() == digest

data = json.loads((SOURCE/'research/passage_certificate.json').read_text())
sqrt_lo = F(141421356237309504880, 10**20)
sqrt_hi = F(141421356237309504881, 10**20)
assert sqrt_lo**2 < 2 < sqrt_hi**2
geometry = {}
for kind, interval in data['half_angle_intervals'].items():
    a, b = map(F, interval)
    assert 0 < a < b < 1
    dmin = 2 + sqrt_lo*(1-b*b)/(1+b*b)
    dmax = 2 + sqrt_hi*(1-a*a)/(1+a*a)
    tmax = 2*b/(1+b*b)
    assert dmin-tmax > 2
    assert dmax < 4 and 0 < tmax < 1
    assert 3*dmax*dmax + 6*tmax*tmax < 54 < 64
    row_norm_sq_upper = dmax*dmax + 2*tmax*tmax
    assert row_norm_sq_upper < 18
    assert 18 < F(9,2)**2 and 18*F(9,2) == 81
    key = 'original_free_loops' if kind == 'original' else 'free_loops'
    lengths = []
    for loop in data[key]:
        vertices = loop['fractional_vertices_scaled']
        delta = [vertices[-1][i] - vertices[0][i] for i in range(3)]
        assert delta == [data['fractional_polyline_scale']*w for w in loop['winding']]
        bound = F(8, data['fractional_polyline_scale'])*sum(
            sum(abs(x-y) for x,y in zip(p,q)) for p,q in zip(vertices,vertices[1:]))
        assert 0 < bound < 20
        lengths.append(str(bound))
    geometry[kind] = {'d_min':str(dmin), 'd_max':str(dmax), 't_max':str(tmax),
                      'sigma_min_lower':str(dmin-tmax),
                      'row_norm_squared_upper':str(row_norm_sq_upper), 'length_upper':lengths}

Co = [[0,-1,-1],[1,0,-1],[1,-1,0]]
Cv = [[0,1,1],[-1,0,1],[-1,-1,0]]
sym_o = [[F(Co[i][j]+Co[j][i],2) for j in range(3)] for i in range(3)]
sym_v = [[F(Cv[i][j]+Cv[j][i],2) for j in range(3)] for i in range(3)]
assert sym_o == [[0,0,0],[0,0,-1],[0,-1,0]]
assert sym_v == [[0,0,0],[0,0,0],[0,0,0]]
# The original symmetric form is -2*y*z >= -(y*y+z*z).
# This gives the claimed accretivity bound d-t on all vectors.

def det3(A):
    return sum((1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else -1)*A[0][i]*A[1][j]*A[2][k]
               for i,j,k in [(0,1,2),(1,2,0),(2,0,1),(0,2,1),(2,1,0),(1,0,2)])

W = [loop['winding'] for loop in data['free_loops']]
G = [[sum(F(w[i])*w[j] for w in W) - (F(1,3) if i==j else 0)
      for j in range(3)] for i in range(3)]
minors = [G[0][0], G[0][0]*G[1][1]-G[0][1]*G[1][0], det3(G)]
assert minors == [F(2,3),F(4,9),F(2,27)]
assert abs(det3(W)) == 1
assert [loop['winding'] for loop in data['original_free_loops']] == [[0,0,-1],[0,1,0]]

B = F(5671,4000)
probe_min_sq, probe_max_sq = F(399**2,80000), F(401**2,80000)
bounds = {}
for name,h,R2 in [('equal_radius',F(1,1000),F(2)),
                  ('full_probe_ratio_window',F(1,4000000),probe_max_sq)]:
    assert B > 2*h and (B-2*h)**2 > R2
    assert 3*h*h < 1
    one_d_kernel_norm_sq = 2*(F(1)-F(1)+F(1,3))/h
    norm_sq = one_d_kernel_norm_sq**3
    lower_o = F(4)/(2*81*20**2*norm_sq)
    lower_v = F(4,3)/(3*81*20**2*norm_sq)
    assert lower_o == h**3/4800
    assert lower_v == h**3/21600
    case = claimed['tube_and_tensor_bounds'][name]
    assert lower_o == F(case['original_K_over_sigma_lower_on_unblocked_plane'])
    assert lower_v == F(case['skew_K_over_sigma_global_lower'])
    bounds[name] = {'K_original_over_sigma':str(lower_o), 'K_skew_over_sigma':str(lower_v)}

Q = F(claimed['inherited_exact_checks']['plane']['max_squared_distance_upper_bound'])
assert Q < F(176,125)**2
assert F(176,125)+F(1,500) == F(141,100)
assert F(141,100)**2 < probe_min_sq
assert F(1,4000)*8 == F(1,500)
strip_half_width = F(1,4000)
b_amplitude = 1/(2*strip_half_width)
assert 2*strip_half_width*b_amplitude == 1
assert (b_amplitude-1)*2*strip_half_width - (1-2*strip_half_width) == 0
# A half-width-a strip can carry b=1/(2a), giving a periodic H(t)-t.

output = {'source_commit':PIN, 'source_head_observed':subprocess.check_output(
    ['git','rev-parse','HEAD'],cwd=SOURCE,text=True).strip(), 'frozen_files':FROZEN,
    'source_hashes_match_pinned_commit':True, 'lattice_interval_bounds':geometry,
    'skew_gram_minors':list(map(str,minors)), 'tensor_bounds':bounds,
    'periodic_strip_corrector_mean_check':True,
    'scope':'Independent rational derivation of constants and hash binding; functional-analysis arguments reviewed in INDEPENDENT_REVIEW.md.'}
for name, digest in FROZEN.items():
    assert hashlib.sha256((FOLLOWUP/name).read_bytes()).hexdigest() == digest
args.output.write_text(json.dumps(output,indent=2)+'\n')
print('PASS: independent interval, length, kernel, Gram, corrector-mean, tensor and source-hash checks.')
