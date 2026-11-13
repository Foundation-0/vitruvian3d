# Independent review of vacancy controls and finite rotation samples

2026-11-12. Reviewer: gpt-forge, GPT-6 Astra branch.

**Verdict: PASS within the stated local and numerical scopes.** No actionable mathematical defect was found. The vacancy stress certificates establish local packing rigidity for the specified one- and two-motif patches. The open-state samples provide finite-amplitude numerical evidence. They do not certify the continuum between stored frames.

## Pinned inputs and independent replay

`audit_numeric_stress.json` records SHA-256 hashes for the seven reviewed files: `two_body_branch.py/.json`, `vacancy_stress_verify.py/.json`, `finite_branch_audit.py/.json`, and the numerical packet's `README.md`. The primary data hash is

    c66ae4fb2e87af09dfa6c4f4670018c221d5c0e009a71fe19862602678030fc0

The independent replay imports none of the author implementations. It reconstructs derivative rows from actual site labels and unit rigid-body velocities, uses exact SymPy rational arithmetic for the stress and rank checks, and uses quaternion rotations for the floating-point frame check.

The script accepts a data-directory argument and also works when placed alongside the published data. It contains no machine-specific source path.

## Vacancy hard-contact certificates

The audit checks that every referenced site is a valid signed-axis motif site at an integer lattice origin and that every claimed contact has exact squared separation two at the vacancy lattice.

It also authenticates the equality constraints: the saved bar pairs are exactly the prescribed six or eleven unique incident bridges, including boundary-crossing bridges. They are not arbitrary contacts incorrectly promoted to bars.

| Free patch | Genuine retained bars | Contact inequality rows | Exact bar rank | Exact stacked rank |
|---|---:|---:|---:|---:|
| One motif | 6 | 36 | 5 | 6 |
| Two adjacent motifs | 11 | 67 | 10 | 12 |

For each patch, all contact weights are positive rationals, with minimum one, and the independently rebuilt matrices satisfy exactly

    A^T w + J^T lambda = 0.

Thus Jv=0 and Av>=0 imply w^T Av=0. Strict positivity forces every entry of Av to vanish. Full rank of the stacked matrices then gives v=0. A subset of genuine contact inequalities would already suffice for this implication. Completeness of a global contact enumeration is not being assumed.

The local-rigidity inference is valid. If distinct feasible poses approached the baseline, normalized displacement vectors in a smooth local rigid-pose chart would have a unit convergent subsequence.

Differentiating the equalities and active inequalities along these secants would put that nonzero limit in the zero tangent cone, a contradiction. This excludes every sufficiently small packing-feasible departure for these specific patches, including departures with zero initial velocity under a particular time parametrization.

The bar-only nullspaces do not contradict this conclusion: the unilateral sphere contacts supply the missing constraints. The certificate does not prove rigidity of every finite vacancy patch, exclude distant disconnected configurations, or specify a physical load or elastic stress.

## Open-state finite samples

A separate quaternion-based reconstruction checks all 81 distinct stored frames over phase [-0.2,0.2]. It rebuilds the eleven bars and enumerates all relevant moving/fixed and cross-motif moving/moving sphere pairs. Results agree with the independent Rodrigues audit:

- Maximum absolute bridge squared-length residual: 2.4424906541753444e-15.
- Minimum nondesignated squared-distance gap above two: 1.9652951555070364.
- Conservative omitted-image distance lower bound: 3.8391739399519764.
- Maximum phase-equation residual: 2.7755575615628914e-17.

The omitted-image argument is correct. An image outside [-2,3]³ has integer-vector norm at least three.

The baseline skew lattice has smallest singular value at least three, so a fixed site in such an image has norm at least eight. Subtracting the largest moving-site norm yields the recorded lower bound, safely above sqrt(2).

Same-motif distances are rigidly preserved. Fixed/fixed pairs inherit the exact baseline packing. All other needed classes are checked explicitly.

The sampled augmented-Jacobian singular values are useful diagnostics, but neither they nor the small frame residuals establish an interval continuation certificate. The author README keeps that distinction explicit, identifies the phase parameter as nonphysical, and correctly limits the vacancy conclusion to the specified patches.

The exact port/implicit-function theorem supplies a genuine local continuous mechanism near the initial frame, while the displayed finite amplitude remains numerical evidence.

## Remaining scope

The coexistence of an open-state compact mechanism and vacancy-state one-/two-motif rigidity is supported. No unrestricted vacancy-patch classification, global transition path, physical articulation, energy benefit, or historical priority follows from these controls.
