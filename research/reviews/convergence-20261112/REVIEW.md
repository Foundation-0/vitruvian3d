# Matched density, access and switching in Vitruvian 3D

Vitruvian 3D gives an exact example in which the same prescribed bar graph and occupied fraction leave different directions open to a moving probe. This is its clearest addition to the Foundation-0 convergence study. Its layer construction also shows why a barrier found in one repeating cell need not remain a bulk barrier when more cells can move independently.

The reviewed main commit is `e7445d56b9ecfd97b408653a2726f3bab70a6e4c`. All eight contributions are merged. Exact contribution heads, source bindings and replay results are in `evidence.json`. The original incentre-lattices contribution was closed without merging after migration. Its clearance checker is byte-identical to the new repository's checker.

| Contribution | Result supported by the proof and certificates | Limit |
|---|---|---|
| [Baseline](https://github.com/Foundation-0/vitruvian3d/pull/1) | Arbitrary separated angle sequences from 40 to 50 degrees integrate into a compatible infinite sphere-and-bar framework. | The phase in each direction depends only on that coordinate. This is prescribed motion. |
| [Passage](https://github.com/Foundation-0/vitruvian3d/pull/4) | Matched-volume configurations have maximum periodic winding rank two versus three for probes within 0.5% of the obstacle-sphere radius. | Rank counts independent directions of unbounded travel within a component. It is not local, spectral or fractal dimension. |
| [Finite rods](https://github.com/Foundation-0/vitruvian3d/pull/3) | Contact capsules with radius below half the sphere radius preserve the passage result and matched occupied fraction. | Hinge housings, articulation and forces are absent. |
| [One-cell transition](https://github.com/Foundation-0/vitruvian3d/pull/2) | Every feasible vacancy-to-FCC transition exceeds volume `14 + 270/6859`. The exact optimum with two bridge rows fixed is `9 + 3sqrt(22)/2`. | The unrestricted one-cell optimum remains between these bounds. |
| [Independent layers](https://github.com/Foundation-0/vitruvian3d/pull/5) | Sequential conversion in an N-layer cell has peak volume `14N + E`, where `E = 3sqrt(22)/2 - 5`. | The period vector moves. Absolute excess remains finite. |
| [Fixed exterior](https://github.com/Foundation-0/vitruvian3d/pull/6) | A homogeneous framework cannot undergo a nonzero translation of finitely many fixed-orientation motifs while all bridges and exterior motifs remain fixed. | Rotations, deformation, altered bridges and moving boundaries change the assumptions. |

The passage proof supplies both sides of the result. Entire periodic planes covered by inflated obstacles block one direction in the original configuration. Explicit connected polygonal loops establish two accessible directions there and three in the skew configuration. Exact interval arithmetic covers every point along these paths and whole parameter intervals. Analytic bounds handle all omitted periodic images.

The layer result is stronger than a favorable numerical path. Its relative peak expansion is `E/(14N)`. At 52 layers this certified upper bound is below the certified lower bound for every one-cell transition under the compared assumptions. Yet the active layer still needs extra normal height `E/sqrt(35)`, approximately 0.34408. The conversion front spans an infinite plane.

The fixed-exterior proof explains the boundary dependence. Summing the exact squared bridge-length changes cancels the linear terms and leaves the nonnegative sum of squared neighboring displacement differences. Preserving every bridge then forces all translations to vanish. At the vacancy state, some infinitesimal bar motions exist but cannot extend to finite motions under these conditions.

This extends [physics0's synthesis](https://github.com/Foundation-0/physics0/blob/3007774b935495335a59f54c5e9b5d5454eab635/notes/convergences-20261112.md): graph observables can stay unchanged while geometric observables change. The new example controls density while changing access. It supplies no constant-volume path between the matched-density configurations.

The mechanics assessment remains unchanged. [The fixed-angle study](https://github.com/Foundation-0/incentre-lattices/pull/100) includes angular stiffness, distributed holding loads and an endpoint instability. Vitruvian's theorems omit that energy law and restrict motif orientations. Their feasible paths and boundary obstruction neither remove that instability nor establish mechanical stability.

No gravitational field equation, spacetime dynamics or gravitational observation is derived. Winding rank concerns probe travel through obstacles and supplies no identification with the dimensions used in the gravity or fractals models.

The background theory is established. [Connelly and Whiteley](https://pi.math.cornell.edu/~connelly/pdf/10.1137_S0895480192229236.pdf) distinguish infinitesimal flexibility, finite rigidity and prestress stability. [Borcea and Streinu](https://scholarworks.smith.edu/csc_facpubs/261/) study related periodic deployment geometry.

[Chen, Upadhyaya and Vitelli](https://arxiv.org/abs/1404.2263) establish nonlinear moving fronts in a different mechanical linkage with prototypes. These comparisons support the scope boundaries. They do not settle historical priority for the particular certificates here.

Eight targeted replays passed on Python 3.14.6, including the independent passage verifier, layer composition, boundary identity and pinned-source audit. Generated receipts were unchanged. Of 88 files bound by the earlier full composition receipt, 87 match current main byte-for-byte. Only README differs.

The earlier 17-command replay and wide passage certificate are retained evidence. This review reuses their source bindings and does not claim to have rerun those full checks.

One further question is justified if transport becomes the target: put a spherical probe of radius r in each exact matched-volume configuration. Use isotropic diffusion with the same positive microscopic diffusivity and reflecting obstacle boundaries, within a chosen connected periodic component.

First establish the boundary regularity needed for reflected diffusion, the invariant measure in the periodic quotient, and existence of a long-time covariance limit. Then test whether the resulting diffusion tensor has rank two in the original component and rank three in the skew component. These assumptions and positive diffusion bounds have not been proved here.

The plane separator and certified open route tubes provide independent controls. A positive tensor eigenvalue needs an analytic lower bound or a numerical error bound, not merely a positive mesh estimate. Failure of a required regularity or homogenization assumption would block this inference. This proposed test connects accessibility to specified dynamics without a gravity claim.

Human mathematical review, a complete priority comparison, practical joints, autonomous propagation and material performance remain unestablished.
