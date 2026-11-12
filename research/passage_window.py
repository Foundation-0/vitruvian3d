"""Finite interval extension of the certified rank2/rank3 passage result.
No new geometric assumption: use the same saved free polylines and a common
64x64 plane-cover witness, verified over closed rational parameter tiles.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import passage_certificate as core

HERE=Path(__file__).parent
ORIGINAL=(F(355,1000),F(364,1000))
SKEW=(F(415,1000),F(426,1000))


def main():
    basepath=HERE/'passage_certificate.json';base=json.loads(basepath.read_text())
    core.ORIGINAL_INTERVAL=ORIGINAL;core.SKEW_INTERVAL=SKEW
    N,triangles,labels=core.plane_discovery(64)
    sphere_table=[];label_ids=[]
    for label in labels:
        if label not in sphere_table:sphere_table.append(label)
        label_ids.append(sphere_table.index(label))
    cells=[]
    for route,lo,hi in [('original',3550,3640),('skew',4150,4260)]:
        for k in range(lo,hi):
            interval=(F(k,10000),F(k+1,10000))
            C=core.C_ORIGINAL if route=='original' else core.C_SKEW
            loops=base['original_free_loops'] if route=='original' else base['free_loops']
            free=core.loops_check(loops,F(5671,4000),C,interval)
            assert free['accepted'],(route,k,free)
            result={'route':route,'interval':[str(x) for x in interval],'free_path_minimum_certified_distance_margin':free['minimum_certified_distance_margin_float']}
            if route=='original':
                core.ORIGINAL_INTERVAL=interval
                plane=core.plane_check(N,triangles,labels,F(399**2,80000));assert plane['accepted'],(k,plane)
                result['plane_distance_squared_upper_bound']=plane['largest_distance_squared_upper_bound_fraction']
            cells.append(result)
            if (k-lo+1)%10==0:print(route,k-lo+1,'/',hi-lo,flush=True)
    core.ORIGINAL_INTERVAL=ORIGINAL;core.SKEW_INTERVAL=SKEW
    volume=core.volume_bracket()
    data={'base_certificate_sha256':hashlib.sha256(basepath.read_bytes()).hexdigest(),'half_angle_intervals':{'original':[str(x) for x in ORIGINAL],'skew':[str(x) for x in SKEW]},'probe_ratio_interval':['199/200','201/200'],'integer_lattice_enclosure_scale':core.D,'fractional_polyline_scale':core.FGRID,'plane_grid_size':N,'plane_sphere_table':sphere_table,'plane_label_ids':label_ids,'parameter_tile_count':len(cells),'cells':cells,'matched_volume_bracket':volume,'angle_intervals_degrees_for_orientation_only':{'original':[2*math.degrees(math.atan(float(s))) for s in ORIGINAL],'skew':[2*math.degrees(math.atan(float(s))) for s in SKEW]},'scope':'Every real parameter and every real segment/triangle point, via exact enclosure on each closed rational tile. Angle decimals are explanatory only.'}
    (HERE/'passage_window.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS',len(cells),'closed rational tiles')


if __name__=='__main__':main()
