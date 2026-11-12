# Centre dynamics: bounded investigation and decision

2026-11-12. GPT-6 Astra developed the exact family. GPT-6.1 Sol independently checked the formulas and replay.

## Decision

**No-go for treating unmodified face-pyramid centre extraction as a robust universal shape regularizer.

Go for a narrowly scoped mathematical investigation of its information loss and symmetry-dependent cycles, if that is independently useful.** This branch now has exact nontrivial structure, explicitly defined geometry: an arbitrarily small tetrahedrally symmetric perturbation of a cube selects a different exact cycle.

The distributed-folding investigation is a separate research direction..

This is a decision about this defined operation, not a proof that all Vitruvian-inspired investigations fail or succeed. Novelty has not been established.

## Reproducible definition

`centre_dynamics.py` and `centre_dynamics.json` contain bounded experiments.

Each step uses the solid's volume centroid as origin, merges genuinely coplanar hull triangles into polygonal faces, joins those faces to the origin, extracts each pyramid's largest-inscribed-ball centre, takes their convex hull, and recentres/rescales to unit RMS hull-vertex radius.

Ties in largest-ball centres use the unique minimum-Euclidean-norm centre in the optimal centre set. Tetrahedra use the independent opposite-face-area formula. Other pyramids use a linear program followed by a minimum-norm active-set calculation..

Main face-merging tolerance is 1e-8 in normalized coordinates. sensitivity runs use 1e-6 and 1e-10. Largest-ball active-set feasibility uses a 2e-8 numerical tolerance. These are floating-point explorations, not arbitrary-precision general dynamics proofs. Exact families below do not rely on long-run numerical convergence.

## Exact family: a singular cube and an alternative cycle

Write S={-1,+1}^3 and p(s)=s1*s2*s3. Consider

    P_epsilon = conv{(1+epsilon*p(s))*s : s in S},  0<epsilon<1/2.

These are two oppositely oriented regular tetrahedra on the cube's vertex rays, with slightly different radii. The volume centroid remains zero by tetrahedral symmetry. At epsilon=0 the hull is the cube with six square faces. For 0<epsilon<1/2 the hull has twelve actual triangular faces. The endpoint 1/2 degenerates to a tetrahedron, so it is excluded.

Each triangular face pyramid is a tetrahedron whose three origin-containing planes depend only on the fixed cube-vertex rays. Its incenter direction therefore does not depend on epsilon. The twelve pyramids are congruent by tetrahedral symmetry, so their incentres also share a radius.

For one representative pyramid, the three side-plane distance equations give an incenter ray proportional to (3,1,1), with suitable signs. Consequently normalized T(P_epsilon) is exactly the regular truncated tetrahedron

    Q_- = conv{ permutations of (3*s1,s2,s3)/sqrt(11) : p(s)=-1 }.

There are twelve distinct vertices. The sign of epsilon selects the opposite tetrahedral orientation. The amplitude of nonzero epsilon is lost in one step.

### The following step and cycle

Work at the unnormalized truncated-tetrahedron scale with vertices permutations of signed (3,1,1). Its four triangular faces have origin distance h_tri=5/sqrt(3) and apothem rho_tri=sqrt(2/3). Its four hexagonal faces have h_hex=sqrt(3), rho_hex=sqrt(6). Each face is regular and its pyramid is rotationally symmetric. The incenter lies on the face-normal axis at distance

    c = h / (1 + rho/sqrt(h^2+rho^2)).

Thus

    c_tri = (5/sqrt(3)) / (1+sqrt(6)/9),
    c_hex = sqrt(3) / (1+sqrt(6)/3),
    c_hex/c_tri = (sqrt(6)-1)/(sqrt(6)+1).

The eight centres form two opposite tetrahedra, with the triangular-face sector farther from the origin. This is P_{-1/sqrt(6)} up to scale. Since 1/sqrt(6)<1/2, the preceding exact result applies again:

    Q_- -> P_{-1/sqrt(6)} -> Q_+ -> P_{+1/sqrt(6)} -> Q_-.

Therefore the normalized map has a four-cycle in a fixed coordinate frame, or a two-cycle of shape classes if overall rigid rotations are identified. This is distinct from the cube–octahedron cycle. Verification gives the predicted radius ratio 0.4202041028867287, observed 0.4202041028867283, and four-step fixed-frame vertex Hausdorff residual 1.22e-15.

The exact formulas received internal independent AI review. Outside mathematical review remains open. It is not claimed to be previously unknown.

## Why universal regularization fails here

As epsilon tends to zero, P_epsilon approaches the cube, but its extracted hull stays Q_-. At the exact cube, the extracted normalized hull is the regular octahedron O=conv{±e_i}.

This is a difference of convex bodies, not merely a mismatch in their vertex counts. In direction (-1,-1,-1)/sqrt(3), Q_- has support 5/sqrt(33), while O has support 1/sqrt(3). Their Hausdorff distance is therefore at least

    5/sqrt(33) - 1/sqrt(3) = 0.2930380105888635 > 0.

The operation is discontinuous at the cube when its input space allows faces to split under arbitrarily small perturbations. Genuine coplanar-face merging does not remove this: the perturbed triangles are actual different faces. Treating every triangulation of the exact cube as separate faces would define a different operation, and would erase the original cube-to-octahedron result.

This is a decisive counterexample to a continuous or perturbation-robust universal regularizer, although not to every possible eventual convergence assertion. The alternative exact cycle separately excludes universal eventual convergence to the cube–octahedron cycle.

## General numerical probes and their limits

Seven initial cases were run for ten iterations: exact cube, unequal-axis octahedron, small and large shear octahedra, sheared cube, random near-octahedron, and random near-cube. Exact cube and unequal-axis octahedron reach the known cycle, as expected.

The random near-octahedron grows from 6 vertices/8 faces to 85 vertices/166 faces by iteration ten. random near-cube reaches 38 vertices/72 faces. Their two-step vertex-set distances remain substantial (about 0.364 and 0.696). These short runs do not establish chaos or nonconvergence.

Symmetry-preserving shear runs can appear to approach another eight-vertex/eight-face cycle. However, the sequence changes under a tighter coplanarity tolerance as numerical asymmetry grows. **No convergence claim is made from those runs.** This is precisely why the exact tetrahedral family is the stronger evidence.

A random cube perturbation sequence from amplitude 1e-2 to 1e-6 approaches the original cube while the sampled support-gap lower bound between the outputs approaches about 0.320918. The exact parity family gives the cleaner proof above, with an invariant gap 0.293038.

## Related prior work and novelty boundary

Joseph O'Rourke's 2015 primary research question, [Tetrahedron incenter iteration tree](https://mathoverflow.net/questions/215068/tetrahedron-incenter-iteration-tree), recursively subdivides tetrahedra through incentres and asks about the limiting hull. It is a related but different operation: it retains a subdivision tree rather than replacing a polyhedron with the hull of one centre per current face pyramid.

The 2013 primary discussion [Tetrahedron insphere iteration](https://mathoverflow.net/questions/144498/tetrahedron-insphere-iteration) instead joins insphere contact points. Its answer develops nonregular two-cycles and a convergence proposition for that map. The result cannot be transferred to the present T.

It does show that incenter-related iterative geometry and nonregular cycles already have substantial prior exploration..

Targeted searches did not identify a source establishing the exact T used here. That is a limited search outcome, not evidence of novelty.

## What remains worth asking

A bounded mathematical question is: can the cycles and discontinuity loci of T be classified on low-dimensional symmetry families, separating metric information lost within a fixed cone fan from changes caused by face combinatorics? The parity family answers one special case and provides exact regression data.

Its payoff would be understanding what the construction actually computes, not assuming it reveals universal proportions..

For a robust practical geometric process, first choose whether face combinatorics are meaningful input. If small faceting changes should be ignored, an explicit coarsening or alternate face measure is necessary. that would be a new model. No evidence here warrants expanding this into an application project yet.

