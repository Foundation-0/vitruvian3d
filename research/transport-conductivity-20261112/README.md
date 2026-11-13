# Certified access gives different conductivity ranks

The existing Vitruvian passage witnesses imply an exact difference in a specified static transport operator. Put the same positive scalar conductivity on the space available to probe centers and zero conductivity in the excluded region. The periodic variational conductivity tensor has rank two in the original configuration and rank three in its exactly volume-matched skew partner.

This preserves the same prescribed bar graph and physical occupied fraction. The operator is defined on accessible probe-center space. It is not the conductivity of the original solid material, and its values are not particle diffusion coefficients.

The positive result uses a finite-energy flux built around each certified open route. Its mean equals that route's physical winding displacement. A flux inequality gives strictly positive conductivity in every direction spanned by those displacements. The original geometry also contains an insulating slab, which admits an explicit periodic potential with zero energy in the blocked direction.

For equal probe and obstacle radii, conservative bounds in the source's physical coordinate units are

\[
K_o/\sigma\succeq P_S/4800000000000,\qquad
K_v/\sigma\succeq I/21600000000000.
\]

Here \(S\) is the plane spanned by the second and third original lattice vectors, and \(P_S\) is Euclidean orthogonal projection. The original tensor vanishes on its normal. These small bounds establish positivity only. They are not estimates of useful performance.

The proof applies to all accessible space and separately to the torus component containing the saved loops, including its full periodic lift. It covers the narrow certified parameter intervals and the entire original probe-radius window. The prescribed finite contact rods leave this accessible domain unchanged under their existing exclusion theorem.

[THEOREM.md](THEOREM.md) supplies the variational definition, flux construction, insulating potential, exact constants and counterexamples to broader claims. [certificate.json](certificate.json) records the arithmetic checks. [sources.json](sources.json) records the primary-source boundary.

The checker binds four source files from Vitruvian commit `e7445d56b9ecfd97b408653a2726f3bab70a6e4c` by SHA-256. It replays the existing exact passage acceptance functions and checks the added bounds. With a checkout containing those files, run:

```sh
python3 check.py --source-root /path/to/vitruvian3d
```

The written proof and arithmetic certificate support the defined static tensor. They establish no Brownian limit, per-particle diffusivity, transport times, practical joints, autonomous motion or gravitational dynamics. The original review files remain unchanged.
