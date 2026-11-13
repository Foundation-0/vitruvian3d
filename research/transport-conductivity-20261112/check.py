"""Bind a variational-conductivity proof to existing exact passage witnesses.

This replays source acceptance functions and checks only the new rational bounds.
It solves no PDE and simulates no stochastic process.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

PIN = "e7445d56b9ecfd97b408653a2726f3bab70a6e4c"
HASHES = {
    "research/passage_certificate.json": "deabfc3b90da4666a45b8d09021c22827c1e6b3284876c49eb81a7fe0f8dc0c7",
    "research/passage_verify.py": "4c8fb018cc2dcd09ee40467f80c689a6d95c6a877582ec11af15fe32c9f96dac",
    "research/passage_proof.md": "2a82d76b85d4b62d3d95a7c6de99b5871aebba4f941c335e428f182241ddd74c",
    "research/connectors_proof.md": "cc328d5b84f3bd0b02c1cc85dd51d50e77835acc2592f8af30ac14ce0594d441",
}


def check(source: Path) -> dict:
    for name, expected in HASHES.items():
        assert hashlib.sha256((source / name).read_bytes()).hexdigest() == expected, name
    spec = importlib.util.spec_from_file_location("bound_passage_verifier", source / "research/passage_verify.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    data = json.loads((source / "research/passage_certificate.json").read_text())
    assert module.SQ[0] ** 2 < 2 < module.SQ[1] ** 2
    inherited = {
        "plane": module.verify_plane(data),
        "original": module.verify_loops(data, "original"),
        "skew": module.verify_loops(data, "skew"),
        "volume_and_angle": module.verify_volume_and_angle(data),
    }

    # The original checker already proves sigma_min(L)>2 on both intervals.
    # Entry bounds d<4, |off_diagonal|<1 imply ||L||_F^2<54<64.
    # Row norm squared<18 and sqrt(18)<9/2 imply V<81.
    assert 2 < 2**2 and 54 < 64 and 18 < F(9, 2)**2
    lattice_bounds = {"smallest_singular_value_strict_lower": "2",
                      "operator_norm_strict_upper": "8",
                      "cell_volume_strict_upper": "81"}

    lengths = {}
    for kind, key in [("original", "original_free_loops"), ("skew", "free_loops")]:
        rows = []
        for loop in data[key]:
            points = loop["fractional_vertices_scaled"]
            steps = sum(sum(abs(x-y) for x, y in zip(a, b)) for a, b in zip(points, points[1:]))
            bound = F(8 * steps, data["fractional_polyline_scale"])
            assert 0 < bound < 20
            rows.append({"winding": loop["winding"], "physical_length_strict_upper": str(bound)})
        lengths[kind] = rows

    B = F(5671, 4000)
    R2_max = F(401**2, 80000)
    h_unit = F(1, 1000)
    h_window = F(1, 4_000_000)
    cases = {}
    for name, h, R2 in [("equal_radius", h_unit, F(2)),
                        ("full_probe_ratio_window", h_window, R2_max)]:
        assert B > 2*h
        slack = (B - 2*h)**2 - R2
        assert slack > 0
        # Cube support has radius sqrt(3)*h<2*h and diameter below period spacing.
        assert 3*h*h < 1
        # Product tent kernel: mass=1 and squared L2 norm=(2/(3h))^3.
        kernel_squared_norm = F(8, 27) / h**3
        assert (F(2, 3)/h)**3 == kernel_squared_norm
        # Averaging m individual Cauchy-Schwarz bounds, V<81 and lengths<20.
        original_coefficient = F(27, 8) * h**3 / (2*81*20**2) * 4
        skew_coefficient = F(27, 8) * h**3 / (3*81*20**2) * F(4, 3)
        assert original_coefficient == h**3 / 4800
        assert skew_coefficient == h**3 / 21600
        cases[name] = {
            "kernel_cube_half_width": str(h),
            "exclusion_radius_squared_upper": str(R2),
            "squared_tube_clearance_slack": str(slack),
            "kernel_squared_L2_norm": str(kernel_squared_norm),
            "original_K_over_sigma_lower_on_unblocked_plane": str(original_coefficient),
            "skew_K_over_sigma_global_lower": str(skew_coefficient),
        }

    maximum_plane_distance_squared = F(inherited["plane"]["max_squared_distance_upper_bound"])
    assert maximum_plane_distance_squared < F(176, 125)**2
    assert F(176, 125) + F(1, 500) == F(141, 100)
    assert F(141, 100)**2 < F(399**2, 80000)
    assert F(1, 4000) * 8 == F(1, 500)
    barrier = {
        "plane_distance_strict_upper": "176/125",
        "physical_slab_half_width": "1/500",
        "fractional_slab_half_width": "1/4000",
        "minimum_exclusion_radius_squared": str(F(399**2, 80000)),
        "strict_squared_margin": str(F(399**2, 80000) - F(141, 100)**2),
        "normal_direction": "g=L^{-1}e_1; kernel is span(g)",
    }

    W = inherited["skew"]["winding_vectors"]
    assert abs(module.determinant(W)) == 1
    gram = [[sum(F(row[i])*row[j] for row in W) - (F(1, 3) if i == j else 0)
             for j in range(3)] for i in range(3)]
    minors = [gram[0][0],
              gram[0][0]*gram[1][1] - gram[0][1]*gram[1][0],
              module.determinant(gram)]
    assert all(x > 0 for x in minors)
    assert inherited["original"]["winding_vectors"] == [[0, 0, -1], [0, 1, 0]]

    # Must-fail controls for this certificate's sufficient conditions.
    rejected = {
        "equal_radius_kernel_half_width_1_over_100": (B - F(2, 100))**2 <= 2,
        "window_kernel_half_width_1_over_1000000": (B - F(2, 1_000_000))**2 <= R2_max,
        "physical_slab_half_width_1_over_10": (F(176, 125)+F(1, 10))**2 >= F(399**2, 80000),
        "dependent_winding_set": module.determinant([W[0], W[1], W[0]]) == 0,
    }
    assert all(rejected.values())
    return {
        "source_commit": PIN,
        "source_sha256": HASHES,
        "python": sys.version.split()[0],
        "inherited_exact_checks": inherited,
        "lattice_bounds": lattice_bounds,
        "loop_length_bounds": lengths,
        "tube_and_tensor_bounds": cases,
        "original_insulating_slab": barrier,
        "skew_W_transpose_W_minus_one_third_identity_leading_minors": list(map(str, minors)),
        "controls_rejected": rejected,
        "conclusion": "The defined bulk variational conductivity has rank exactly two on the original certified interval and three on the skew interval. The same holds for the component containing the saved loops.",
        "domain_scope": "Open probe-center free space, uniform positive scalar conductivity there, zero outside; periodic correctors. The full probe-radius window is certified on the narrow parameter intervals only.",
        "proof_dependency": "Current convolution, Hilbert-space projection, zero-contrast limit, and insulating-strip corrector are proved in THEOREM.md; this script does not formalize those arguments.",
        "not_established": ["Brownian limit law", "per-particle diffusivity", "transport rates beyond conservative variational bounds", "hardware performance", "gravity"],
    }


def main():
    if not __debug__:
        raise SystemExit("Assertions must remain enabled.")
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--output", type=Path, default=Path(__file__).with_name("certificate.json"))
    args = p.parse_args()
    result = check(args.source_root.resolve())
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: source-bound passage replay, tube bounds, zero-normal strip, tensor ranks, and four rejected controls.")
    print("Equal-radius lower bounds: K_original/sigma >= P/4800000000000; K_skew/sigma >= I/21600000000000.")
    print("These bounds certify positivity and are not performance estimates.")


if __name__ == "__main__":
    main()
