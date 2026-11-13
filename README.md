# Vitruvian 3D

This project studies a three-dimensional Vitruvian idea: understand a body through the space it occupies, the routes available to it, and the collective changes of posture it can join. It extends [Foundation-0/incentre-lattices](https://github.com/Foundation-0/incentre-lattices) using spheres, octahedral motifs, and fixed contact bars.

**The investigation supports further research into accessible space and coordinated switching.** Equal occupied fraction can coexist with different directions of unbounded access. Allowing layers to switch separately reduces the relative expansion cost. From a homogeneous starting state, fixing the exterior prevents that motion in the translation-only model.

These are geometric results with written proofs and internal AI review. Human mathematical review, historical priority, mechanical performance, and anatomical significance remain unestablished. [Credit and scope](NOTICE.md), [research decision](research/status.md), and [related work](research/related-work.md) record those boundaries.

## Start here

Merged contributions are on `main`. The complete research snapshot is on **[gpt-forge/research-current](https://github.com/Foundation-0/vitruvian3d/tree/gpt-forge/research-current)**, assembled from the separately reviewable contributions below.

Assembly makes the work usable together and does not indicate maintainer acceptance. Exact source hashes and replay results are in [the composed receipt](research/reviews/composed.json).

| Result | What the evidence establishes | Owning contribution |
|---|---|---|
| Equal-density access | Exactly two versus three independent periodic directions for near-spherical probes, throughout a certified parameter window | [Certified passage](https://github.com/Foundation-0/vitruvian3d/pull/4) |
| One-cell switching | A positive expansion lower bound and an exact optimum with two bridge directions held fixed | [Transition barrier](https://github.com/Foundation-0/vitruvian3d/pull/2) |
| Distributed switching | Arbitrary layer compositions remain feasible, with relative peak expansion proportional to `1/N` | [Layer construction](https://github.com/Foundation-0/vitruvian3d/pull/5) |
| Rotating motifs | At the open skew posture, a finite patch has one local configuration degree of freedom per internal bridge, with its exterior fixed | [Proof and certificate](research/rotations/proof.md) |
| Fixed exterior | From a homogeneous state, no bounded patch can translate nontrivially while preserving all bridges and fixed motif orientations | [Boundary obstruction](https://github.com/Foundation-0/vitruvian3d/pull/6) |
| Finite rods | Capsules of radius less than half the obstacle radius preserve the relevant clearances and add constant material volume | [Connector geometry](https://github.com/Foundation-0/vitruvian3d/pull/3) |
| Figures and initial illustrations | Reproducible scientific plots and the original interactive fragments with source hashes | [Visual record](https://github.com/Foundation-0/vitruvian3d/pull/7) |

The earlier [distributed-fold proof](research/nonuniform_motion_proof.md) covers arbitrary separated phase sequences from 40 to 50 degrees. The [centre-extraction result](research/centre_dynamics_report.md) gives an alternative cycle and discontinuity, ruling out universal regularization. The [vacancy identification](research/vacancy_independent.txt) matches a known packing. These results remain part of the record.

## Reproduce the complete snapshot

```sh
git clone https://github.com/Foundation-0/vitruvian3d.git
cd vitruvian3d
git switch gpt-forge/research-current
python3.13 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python research/replay_all.py
```

The runner executes 22 certificate, control, and stored-frame audit commands, including the rotational follow-up. It also replays the separate pinned-source passage audit. Assertions must remain enabled. Scripts write JSON receipts beside their source.

The numerical frame audit records runtime metadata, so its output may differ between environments. The complete window check takes about a minute on the review machine, with hardware-dependent runtime.

The numerical discovery scripts remain available, but are not acceptance tests for the exact theorems. Figure reproduction has its own pinned requirements in `visualizations/`. The archived HTML fragments require their original host styles and predate the later proofs.

## Continue from a concrete question

For a physical mechanism, define its joints, boundary freedom, and actuation before simulation. For a mathematical extension, state which assumptions change and what result would count as progress. The research decision distinguishes proved local rotational motion from the unresolved dense-state transition, elongated probes, and unrestricted optimum.
