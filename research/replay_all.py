"""Replay the published mathematical composition from an assembled checkout.

This runner executes existing checks. It is not a new mathematical verifier.
"""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CHECKS = (
    "nonuniform_motion.py",
    "vacancy_independent.py",
    "vacancy_tangent.py",
    "centre_extraction.py",
    "volume_barrier.py",
    "volume_path_certificate.py",
    "volume_restricted_optimum.py",
    "volume_layer_certificate.py",
    "passage_verify.py",
    "passage_window_verify.py",
    "connectors_exact.py",
    "connectors_verify.py",
    "connectors_layer_review.py",
    "connectors_volume_review.py",
    "connectors_passage_review.py",
    "boundary_certificate.py",
    "connectors_boundary_review.py",
    "rotations/certificate.py",
    "rotations/single_platform.py",
    "rotations/open_patch_review.py",
    "rotations/numerics/vacancy_stress_verify.py",
    "rotations/numerics/finite_branch_audit.py",
)
STABLE_RECEIPTS = (
    "connectors_exact.json", "passage_verify.json", "passage_window_verify.json",
    "connectors_passage_review.json", "volume_path_certificate.json",
    "volume_layer_certificate.json", "connectors_layer_review.json",
    "boundary_certificate.json", "connectors_boundary_review.json",
    "rotations/certificate.json", "rotations/single_platform.json",
    "rotations/open_patch_review.json", "rotations/numerics/vacancy_stress_verify.json",
)


def main():
    if not __debug__:
        raise SystemExit("Assertions must be enabled. Do not use python -O.")
    scripts = [ROOT / "research" / name for name in CHECKS]
    scripts.append(ROOT / "research/reviews/passage_audit.py")
    missing = [str(path.relative_to(ROOT)) for path in scripts if not path.is_file()]
    if missing:
        raise SystemExit("Use the complete research-current branch. Missing: " + ", ".join(missing))
    receipts = [ROOT / "research" / name for name in STABLE_RECEIPTS]
    before = {path: path.read_bytes() for path in receipts}
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("PYTHONOPTIMIZE", None)
    for index, script in enumerate(scripts, start=1):
        print(f"[{index}/{len(scripts)}] {script.relative_to(ROOT)}", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, env=env, check=True)
    changed = [str(path.relative_to(ROOT)) for path in receipts if path.read_bytes() != before[path]]
    if changed:
        raise SystemExit("Review regenerated receipt differences: " + ", ".join(changed))
    print(f"PASS: {len(CHECKS)} composed checks, pinned-source audit, and unchanged acceptance receipts.")


if __name__ == "__main__":
    main()
