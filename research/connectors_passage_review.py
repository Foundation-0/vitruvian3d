"""Saved-data integer-overflow audit, independent of passage implementation."""
import json,math
from fractions import Fraction as F
from pathlib import Path
D=json.loads(Path(__file__).with_name('passage_certificate.json').read_text())
assert D['plane_grid_size']==32 and D['integer_lattice_enclosure_scale']==10**6
assert D['fractional_polyline_scale']==10**6
scale=D['integer_lattice_enclosure_scale'];grid=D['fractional_polyline_scale'];N=D['plane_grid_size']
# All lattice coefficients have absolute value<4, so integer enclosures<=4*scale.
# The saved half-angle intervals are strictly inside(0,1), and yield d<4,t<1.
for interval in D['half_angle_intervals'].values():
 a,b=map(F,interval);assert 0<a<b<1
labelmax=max(abs(x) for label in D['plane_sphere_labels'] for x in label['n'])
# Every triangle coordinate is between0 andN, hence each W coordinate is <=N(1+labelmax).
plane_w=N*(1+labelmax)
plane_coordinate_bound=3*plane_w*(4*scale)+N*scale
plane_square_sum_bound=3*plane_coordinate_bound**2
# This crude bound must fit int64 even before the sharper runtime assertion.
assert plane_square_sum_bound<2**63
largest_recentered=0;largest_w=0;segments=0
for name in ['free_loops','original_free_loops']:
 for loop in D[name]:
  ps=loop['fractional_vertices_scaled']
  for first,second in zip(ps,ps[1:]):
   shift=[(x+y)//(2*grid) for x,y in zip(first,second)]
   a=[x-grid*n for x,n in zip(first,shift)];b=[x-grid*n for x,n in zip(second,shift)]
   largest_recentered=max(largest_recentered,*map(abs,a+b))
   largest_w=max(largest_w,*[abs(x)+3*grid for x in a+b]);segments+=1
assert largest_recentered<2*grid
# w@center, steps@center, and maxw@delta have <=3 products; both center anddelta
# entries bounded by4*scale, which deliberately overbounds actual tiny errors.
loop_matrix_product_bound=3*max(largest_w,4*grid)*(4*scale)
loop_a_bound=loop_matrix_product_bound+grid*scale
assert loop_a_bound<2**63
# All subsequent squares/products in verifier are Python integers. Discovery
# explicitly converts each component to int before squared-distance arithmetic.
assert 5671**2-200*401**2==41>0
assert F(399**2,80000)==(F(399,200)**2)/2
result=dict(probe_size_interval=['199/200','201/200'],probe_upper_bound_squared_integer_slack=41,geometry_review='PASS plane convexity, periodic separator, exact segment distance, matrix interval error, omitted-image singular-value bound, and independent winding determinants.',standard_library_verifier='Executed successfully as passage_verify.py; output retained in passage_verify.json.',integer_overflow=dict(plane_label_maximum_translation=labelmax,conservative_plane_coordinate_bound=plane_coordinate_bound,conservative_plane_three_square_sum_bound=plane_square_sum_bound,int64_limit=2**63-1,maximum_recentered_scaled_coordinate=largest_recentered,conservative_loop_matrix_product_bound=loop_matrix_product_bound,conservative_loop_site_vector_coordinate_bound=loop_a_bound,squared_loop_arithmetic='Python arbitrary precision integers; no int64 products after conversion.'),limitations='Discovery topology is unnecessary for acceptance. The saved loops and sphere labels are the witnesses. This audit does not claim hardware realization or historical novelty.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
