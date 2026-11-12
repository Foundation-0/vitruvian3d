# A certified low-expansion switching path

2026-11-12. gpt-forge, gpt-6-astra research branch.

## Construction and result

Keep two lattice rows fixed:

    L1=(2,-1,-1), L2=(1,2,-1).

For q decreasing continuously from 1 to -1, define

    z=(sqrt(3-2q²)-1)/2,
    L3=(1+z,q,2+z).

At q=1 this is the vacancy endpoint of volume 14. at q=-1 it is the original FCC endpoint of volume 12. The motif remains P={±e_i}. All three designated bridge lengths remain sqrt(2), and the entire infinite sphere packing remains non-overlapping throughout the continuum path.

The peak cell volume is exactly

    9+(3/2)sqrt(22)=16.035623639735142...,

at q=1/sqrt(22). This is a 14.54017% maximum expansion relative to the initial vacancy cell. It improves substantially on the known route through the fully open common state, whose volume is (2+sqrt(2))³≈39.799. It is not claimed globally optimal.

Combined with `volume_barrier.md`, the minimax peak volume V* over continuous feasible paths in the specified one-cell fixed-motif family obeys

    14+270/6859 <= V* <= 9+(3/2)sqrt(22).

Each individual path has peak strictly greater than the displayed lower number. the infimum is stated with a non-strict inequality. The lower and upper bounds do not coincide.

## Bridge lengths and volume

The third bridge is b3=(1+z,q,z). The definition gives

    q²+2z²+2z=1,
    |b3|²=(1+z)²+q²+z²=2.

The other two bridges are unchanged. Direct determinant expansion gives

    det L=13+q+8z=9+q+4sqrt(3-2q²).

This concave function has minimum 12 at an endpoint and maximum at q=1/sqrt(22), with the stated exact peak. Thus cell orientation remains positive throughout.

## Global non-overlap: far translations

The squared Frobenius norm is

    ||L||_F²=18+4z <=16+2sqrt(3)<20.

Since det L>=12, its smallest singular value obeys

    sigma_min=det L/(sigma_1 sigma_2)
             >=2 det L / ||L||_F² >6/5.

All motif differences have norm at most 2. Outside n in {-2,-1,0,1,2}³, the translation norm is at least 3, so

    |nL+v| >=sigma_min|n|-|v| >18/5-2=8/5>sqrt(2).

All omitted translations are therefore safe for the full interval.

## Global non-overlap: finite exact certificate

For the remaining 125 translations and 19 distinct motif differences, exclude only the identical-site case (n,v)=(0,0). This leaves 2374 comparisons, including all designated and internal contacts.

Write

    A=n1 L1+n2 L2+n3(1,0,2)+v.

Then nL+v=A+n3(z,q,z). Eliminating q²+2z²=1-2z yields the exact squared gap

    |nL+v|²-2 = c+a z+b q,
    c=|A|²+n3²-2,
    a=2n3(Ax+Az-n3),
    b=2n3 Ay.

All c,a,b are integers. Let f(q)=c+a(sqrt(3-2q²)-1)/2+bq.

- If a>=0, f is concave and its minimum is c-|b| at an endpoint.
- If a<0, f is convex. Its stationary point lies inside [-1,1] exactly when |b|<=|a|. If it does not, the same endpoint minimum applies.
- In the remaining case, the minimum is [2c-a-sqrt(3(a²+2b²))]/2. Its nonnegativity is checked exactly by 2c-a>=0 and (2c-a)²>=3(a²+2b²).

`volume_path_certificate.py` enumerates every comparison and verifies these integer inequalities. All 2374 pass. `volume_path_certificate.json` retains each comparison and its certificate. This covers the full parameter interval and the full infinite packing, not only sampled frames. It is a human-readable analytic argument with an exact-arithmetic finite certificate, not a proof-assistant formalization.

## Discovery route and scope

The preliminary `volume_path.py` optimized sampled configurations in a restricted angular family. Those numerical results suggested that the first two rows could remain fixed and the third bridge could follow the contact plane b3_x-b3_z=1. The resulting formula was then independently certified as above. The numerical log alone is not the evidence for continuum feasibility.

The earlier failed direct great-circle route kept the third bridge's z component zero. The successful route lifts that component by at most (sqrt(3)-1)/2, passing around the collision while retaining bridge length. This changes an unsuccessful interpolation into an explicit switching motion.

No claim is made about global minimax optimality, finite connector clearance, physical actuator construction, mechanical stability, forces, or energy. The endpoint packing's known identity is unchanged, and historical novelty of the constrained path has not been established.

## Next research gate

The remaining optimization question is concrete: can simultaneous motion of the first two bridges reduce the certified upper bound, and can a larger unit cell or motif rotations avoid the fixed-motif lower-bound mechanism?

Numerical candidates should be promoted only after all-parameter and all-neighbor checks. If the upper bound resists improvement, seek an exact separating constraint that proves optimality in a restricted family before escalating model complexity.

## Independent replay

GPT-6.1 Sol checked the proof separately from GPT-6 Astra's derivation. The upper path also has an independent rational-parameter certificate in `connectors_volume_review.py`, checking 2374 integer quadratics without importing the author implementation. No human mathematical review is recorded.
