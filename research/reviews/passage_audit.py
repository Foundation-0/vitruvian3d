"""Read-only audit of the pinned published passage source and exact witnesses.

Only this script's sibling JSON is written. The reviewed checkout is not mutated.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import types
from fractions import Fraction as F

PIN='3dbc5c1a8ca42b0cc843546f18eece0a487f118e'
REPO=Path(__file__).resolve().parents[2]
CONNECTOR_PIN='0bb354ee74f48e921a1350c3d6d7dd4063397639'
HERE=Path(__file__).parent


def ints(values):
    if isinstance(values,list):return all(ints(v) for v in values)
    return type(values) is int


def main():
    research=REPO/'research'
    def pinned(name):
        return subprocess.check_output(['git','show',PIN+':research/'+name],cwd=REPO)
    names=['passage_proof.md','passage_verify.py','passage_window_verify.py','passage_certificate.json','passage_window.json']
    hashes={name:hashlib.sha256(pinned(name)).hexdigest() for name in names}
    base=json.loads(pinned('passage_certificate.json'))
    window=json.loads(pinned('passage_window.json'))
    assert type(base['integer_lattice_enclosure_scale']) is int and base['integer_lattice_enclosure_scale']>0
    assert type(base['fractional_polyline_scale']) is int and base['fractional_polyline_scale']>0
    for name in ['original_free_loops','free_loops']:
        for loop in base[name]:
            assert ints(loop['fractional_vertices_scaled']) and ints(loop['winding'])
    for label in base['plane_sphere_labels']+window['plane_sphere_table']:
        assert ints(label['n']) and ints(label['p'])
    for data in [base,window]:
        for interval in data['half_angle_intervals'].values():
            a,b=map(F,interval);assert 0<a<b<1
    assert hashes['passage_certificate.json']==window['base_certificate_sha256']
    module=types.ModuleType('pinned_passage_verifier')
    module.__file__=str(research/'passage_verify.py')
    exec(compile(pinned('passage_verify.py'),module.__file__,'exec'),module.__dict__)
    # Functions do not write outputs; do not invoke published main().
    narrow={'plane':module.verify_plane(base),'original':module.verify_loops(base,'original'),
            'skew':module.verify_loops(base,'skew'),'volume':module.verify_volume_and_angle(base)}
    assert module.SQ[0]**2<2<module.SQ[1]**2
    assert F(5671,4000)**2>F(401**2,80000)
    connector_hash=hashlib.sha256(subprocess.check_output(['git','show',CONNECTOR_PIN+':research/connectors_proof.md'],cwd=REPO)).hexdigest()
    assert connector_hash=='cc328d5b84f3bd0b02c1cc85dd51d50e77835acc2592f8af30ac14ce0594d441'
    output={'reviewed_passage_commit':PIN,'reviewed_repository':'https://github.com/Foundation-0/vitruvian3d','reviewed_connector_commit':CONNECTOR_PIN,
            'source_read_method':'git show at exact pinned commit, independent of moving checkout HEAD',
            'published_file_sha256':hashes,
            'published_clone_mutated':False,'integer_witness_types_and_positive_scales_checked':True,
            'narrow_exact_function_replay':narrow,
            'wide_replay_evidence':'Identical-hash staging verifier/witnesses replayed all 200 tiles successfully; reviewed interval-cover and hash-binding logic.',
            'pinned_connector_proof_sha256':connector_hash,
            'mathematical_review':'PASS continuum plane obstruction, connected-loop winding ranks, all-image control, matched-volume monotonicity, shape sandwich, rod exclusion and disjoint added material.',
            'dependency_note':'The pinned passage checkout references connectors_proof.md but does not contain it; the separate connector publication supplies that corollary dependency.',
            'scope':'No novelty, optimal throat, arbitrary probe shape, material mechanics, or hardware claim established.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
