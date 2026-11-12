"""Small independent incenter-extraction check; numerical examples, not a novelty claim.

Uses the tetrahedron face-area formula, independently of the axial closed form.
Requires numpy, already present in the research environment.
"""
import itertools
import json
from pathlib import Path

import numpy as np


def tetra_incenter(vertices):
    weights = []
    for i in range(4):
        face = np.delete(vertices, i, axis=0)
        weights.append(np.linalg.norm(np.cross(face[1] - face[0], face[2] - face[0])) / 2)
    return np.average(vertices, axis=0, weights=weights)


SIGNS = np.array(list(itertools.product([-1., 1.], repeat=3)))


def extract(basis):
    return np.array([tetra_incenter(np.vstack([np.zeros(3), s[:, None] * basis]))
                     for s in SIGNS])


def normalized_pair_squared(points):
    points = points - points.mean(axis=0)
    points = points / np.sqrt(np.mean(np.sum(points ** 2, axis=1)))
    return np.sort([float(np.sum((points[i] - points[j]) ** 2))
                    for i in range(len(points)) for j in range(i)])


def main():
    axial = []
    for axes in [[1, 1, 1], [2, 1, 3], [0.25, 4, 1.5]]:
        axes = np.array(axes, dtype=float)
        centers = extract(np.diag(axes))
        inverse = 1 / axes
        radius = 1 / (inverse.sum() + np.linalg.norm(inverse))
        error = float(np.max(np.abs(centers - SIGNS * radius)))
        assert error < 1e-12
        axial.append({'axes': axes.tolist(), 'closed_form_radius': float(radius),
                      'max_face_area_formula_error': error, 'centers': centers.tolist()})
    basis = np.array([[2., 0., 0.], [.4, 1., 0.], [0., 0., 3.]])
    centers = extract(basis)
    pairs = normalized_pair_squared(centers)
    cube_pairs = normalized_pair_squared(SIGNS)
    difference = float(np.max(np.abs(pairs - cube_pairs)))
    assert difference > .01
    result = {
        'status': 'Numerical verification of exact axial deduction; illustrative sheared counterexample.',
        'construction': 'Eight tetrahedra conv(0, +/-v1, +/-v2, +/-v3); one incenter per sign triple.',
        'axial_examples': axial,
        'shear_example': {'basis_rows': basis.tolist(), 'centers': centers.tolist(),
                          'rms_normalized_squared_pair_distances': pairs.tolist(),
                          'cube_squared_pair_distances': cube_pairs.tolist(),
                          'max_sorted_pair_distance_squared_difference_from_cube': difference},
        'limitations': 'No theorem for general polyhedra, no recursive stability test, no literature novelty search.'}
    destination = Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'axial_max_errors': [x['max_face_area_formula_error'] for x in axial],
                      'shear_distance_spectrum_difference': difference,
                      'output': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
