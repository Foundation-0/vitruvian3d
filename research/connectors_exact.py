"""Exact radius certificates for contact rods in a nonoverlapping equal-sphere framework."""
import json
from pathlib import Path
import sympy as s
r,q,a=s.symbols('r q a',positive=True)
# A segment has length2r. Distance to its nearest endpoint is <=sqrt(r²+d²)
# when d is distance to the segment and closest point lies inside.
probe_limit=s.sqrt(q*(2*r+q))-q
assert s.simplify(r*r+(q+probe_limit)**2-(r+q)**2)==0
# Equal-probe containment threshold.
assert s.simplify(probe_limit.subs(q,r)/r-(s.sqrt(3)-1))==0
# Nonincident rods: weighted mean of four endpoint separations gives
# (2r)² <= d²+4r²[s(1-s)+t(1-t)] <=d²+2r².
# Thus d>=sqrt2*r and equal rods of radius<=r/sqrt2 do not penetrate.
# Incident rods: endpoint separation>=2r ->angle>=pi/3. If a point belongs
# to both capsules outside their shared endpoint sphere, its radius must
# be<=a/sin(angle/2)<=2a. For a<=r/2 overlap stays in shared sphere.
result=dict(model='Capsules about all prescribed contact-center segments of length2r, including internal octahedral edges and intermotif bridges. Spheres remain full opaque obstacle balls.',capsule_contained_in_physical_endpoint_spheres='Only zero rod radius. Any a>0 protrudes near the touching point.',probe_exclusion_redundancy=dict(general_max_rod_radius=str(probe_limit),equal_probe_ratio=str(s.sqrt(3)-1),equal_probe_normalized_radius=str((s.sqrt(3)-1)/s.sqrt(2)),equal_probe_normalized_radius_decimal=float((s.sqrt(3)-1)/s.sqrt(2))),nonincident_rods=dict(minimum_centerline_separation='sqrt(2)*r',equal_rod_nonpenetration_radius='r/sqrt(2)',normalized_radius=.5),incident_rods=dict(minimum_angle='pi/3',overlap_region_contained_in_shared_junction_sphere_if='a<=r/2',normalized_radius=float(1/(2*s.sqrt(2)))),simultaneous_universal_certificate=dict(radius_ratio='a/r < 1/2',normalized_radius_upper_bound=float(1/(2*s.sqrt(2))),consequence='No extra equal-radius-probe obstruction; no nonincident rod intersection; incident rod overlaps confined to shared junction spheres. Applies to every globally nonoverlapping equal-sphere contact framework, including all exact separated nonuniform phases.'),limitations='Pure solid geometry. Junction articulation, finite wall thickness, bending, friction, loads, actuation, manufacturing and arbitrary solid motif faces are not modeled. Exact tangencies have no guaranteed tolerance against endpoint sphere growth or generic changes that destroy original contacts.')
h=r-s.sqrt(r*r-a*a)
neck=2*s.pi*(a*a*h-r*h*h+h*h*h/3)
neck_alt=2*s.pi*r*h*h-4*s.pi*h*h*h/3
assert s.simplify(neck-neck_alt)==0
assert s.simplify((s.Rational(1,2)+s.Rational(1,4))**2-s.Rational(1,4)*(2+s.Rational(1,4)))==0
result['composite_volume']=dict(contact_rods_per_cell=15,spheres_per_cell=6,neck_height='h=r-sqrt(r²-a²)',added_volume_per_rod='2*pi*r*h²-(4*pi/3)*h³',total_solid_volume='8*pi*r³+15*(2*pi*r*h²-(4*pi/3)*h³)',condition='a/r<1/2; exterior rod portions do not overlap one another or unrelated spheres',consequence='Matched cell volumes retain exactly matched composite solid fractions.',probe_ratio_band='For a/r<=1/2, every q/r>=1/4 satisfies rod exclusion redundancy.')
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
