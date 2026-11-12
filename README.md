# Vitruvian 3D

This project explores a three-dimensional Vitruvian idea: understand a body through the space it occupies, the movements available to it, and the collective changes of posture it can join. It builds on [Foundation-0/incentre-lattices](https://github.com/Foundation-0/incentre-lattices).

The current models use spheres, octahedra, and contact bars. They establish no law of human proportions or biological mechanism. The useful results so far concern accessible space and constrained motion.

## Research record

| Investigation | Current evidence | Entry point |
|---|---|---|
| Distributed folding | Written infinite-family proof and 488 exact nearby-pair comparisons | [Construction](research/nonuniform_motion_proof.md) |
| Movement at equal density | Independent numerical routes give two periodic directions versus three | [Transport](research/transport_probe.txt) |
| Vacancy endpoint | Exact match to a known packing with one seventh of FCC sites vacant | [Prior-art identification](research/vacancy_independent.txt) |
| Required expansion | Strict local volume minimum in the one-cell fixed-motif model | [Volume obstruction](research/vacancy_tangent_report.md) |
| Centre extraction | Exact alternative cycle and discontinuity rule out universal regularization | [Centre dynamics](research/centre_dynamics_report.md) |
| Failed or limited ideas | Independent phases can overlap, a naive spatial phase fails closure, and fitting differs from escape | [Exploration notes](research/abstract-directions.md) |

These are research records with explicit scopes. Exact arguments have internal AI review, numerical calculations have finite tolerances, and historical priority remains unresolved except where a known source is identified. Later contributions should update the table as the evidence changes.

## Reproduce

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python research/nonuniform_motion.py
python research/vacancy_independent.py
python research/vacancy_tangent.py
python research/centre_extraction.py
```

The scripts write their JSON evidence beside the source. Run with assertions enabled. Full numerical explorations are `transport_probe.py`, `transport_independent.py`, `centre_dynamics.py`, `void_probe.py`, and `chirality_probe.py`. The transport grid reaches 96 cubed samples and takes longer than the exact arithmetic checks.

## Continuing the work

The next questions are a continuum certificate for the transport difference, a quantified minimum-expansion transition, and finite connector clearance. Keep numerical discoveries separate from accepted certificates, preserve failed approaches, and state the model assumptions in every result.

[NOTICE.md](NOTICE.md) records credit, AI assistance, the pinned source model, and evidence limits. Code is MIT. Text and original figures are CC BY 4.0.
