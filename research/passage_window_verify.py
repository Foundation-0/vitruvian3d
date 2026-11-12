"""Stdlib-only replay of the broadened continuum parameter window.
Uses the independent passage verifier, never the discovery implementation.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import passage_verify as independent

HERE=Path(__file__).parent


def main():
    window=json.loads((HERE/'passage_window.json').read_text())
    raw=(HERE/'passage_certificate.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==window['base_certificate_sha256']
    base=json.loads(raw)
    data={**base,'half_angle_intervals':window['half_angle_intervals'],'plane_grid_size':window['plane_grid_size'],'plane_sphere_labels':[window['plane_sphere_table'][i] for i in window['plane_label_ids']]}
    assert window['probe_ratio_interval']==['199/200','201/200']
    assert window['plane_grid_size']==64
    intervals={route:tuple(map(F,value)) for route,value in window['half_angle_intervals'].items()}
    totals={'original':{'tiles':0,'segment_site_checks':0,'plane_vertex_checks':0},'skew':{'tiles':0,'segment_site_checks':0}}
    for route in ['original','skew']:
        cells=[x for x in window['cells'] if x['route']==route]
        assert cells and F(cells[0]['interval'][0])==intervals[route][0]
        assert F(cells[-1]['interval'][1])==intervals[route][1]
        previous=intervals[route][0]
        for cell in cells:
            lo,hi=map(F,cell['interval']);assert lo==previous and hi>lo;previous=hi
            selected={**data,'half_angle_intervals':{**data['half_angle_intervals'],route:[str(lo),str(hi)]}}
            free=independent.verify_loops(selected,route)
            totals[route]['tiles']+=1;totals[route]['segment_site_checks']+=free['exact_segment_site_checks']
            if route=='original':
                plane=independent.verify_plane(selected)
                totals[route]['plane_vertex_checks']+=plane['vertex_checks']
    volume=independent.verify_volume_and_angle(data)
    assert sum(x['tiles'] for x in totals.values())==window['parameter_tile_count']==200
    result={'totals':totals,'volume_and_angle':volume,'probe_ratio_interval':window['probe_ratio_interval'],'half_angle_intervals':window['half_angle_intervals'],'conclusion':'Same rank2/rank3 and matched-volume theorem holds throughout the widened closed parameter intervals; original window spans over0.9degrees and includes exact40deg.'}
    (HERE/'passage_window_verify.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
