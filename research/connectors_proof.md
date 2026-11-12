# Finite contact rods and robust positive transport

2026-11-12. gpt-forge, gpt-6.1-sol.

The positive transport routes survive finite contact rods. A uniform geometric bound applies to every nonoverlapping configuration in the exact separated-phase motion family. The bound includes the twelve internal octahedron edges and the intermotif bridges. It is a theorem about solid geometry, with mechanical behavior excluded.

## Model

Put a closed ball of radius r at each site. Distinct site centers are at least 2r apart. For every prescribed contact edge, whose endpoints have separation exactly 2r, add the closed capsule of radius a around the straight center segment. The spheres remain physical solids. Overlap of a rod with its own endpoint spheres is intentional. Incident rods may merge inside their shared junction sphere.

A probe is a ball of radius q. A route describes its center. Collision-free here permits exact boundary tangency unless a strict margin is stated.

A rod of any positive radius protrudes outside the union of its two physical endpoint spheres near their touching point. It would be incorrect to claim that the physical spheres contain a finite rod.

## Exact probe exclusion bound

Let K be the capsule around a segment of length 2r. The probe-center exclusion set of this rod is a capsule of radius q+a. The endpoint spheres exclude probe centers from balls of radius r+q.

The rod exclusion is contained in the union of these two endpoint exclusion balls if and only if

    a <= sqrt(q*(2r+q)) - q.

To prove sufficiency, a point projecting onto the segment interior has distance at most sqrt(r^2+(q+a)^2) from the nearer endpoint. This is at most r+q under the stated inequality. A point whose closest segment point is an endpoint is covered directly, since the same inequality implies a<r. Necessity follows by taking a point at perpendicular distance q+a from the segment midpoint.

For q=r, the exact threshold is

    a <= (sqrt(3)-1)*r.

In the normalized model r=1/sqrt(2), this is a<=0.5176380902050415. Below this threshold, every sphere-only equal-size-probe route remains valid with all contact rods added. This statement holds for every possible route, rather than depending on one small numerical route margin.

## Rod clearance against unrelated spheres

The preceding containment argument with q=r also applies to any unrelated physical sphere center. Nonoverlap with the rod's endpoint spheres therefore prevents rod penetration into any unrelated sphere whenever a<=(sqrt(3)-1)r. Exact boundary tangency can occur at equality. The rod can enter its own endpoint spheres as intended.

## Nonincident rods

Consider two center segments with four distinct endpoints. All four cross-endpoint distances are at least 2r. Let x and y be closest points on the two segments, expressed as convex combinations with weights s and t in [0,1]. Independently choose random endpoints A and B with means x and y. Then

    4r^2 <= E|A-B|^2
           = |x-y|^2 + 4r^2*s*(1-s) + 4r^2*t*(1-t)
           <= |x-y|^2 + 2r^2.

Thus their centerline separation is at least sqrt(2)r. Equal rods cannot penetrate each other when a<=r/sqrt(2). In normalized units this gives a<=0.5. A strict inequality supplies a positive inter-rod clearance.

## Incident rods

Two contact segments sharing a center have lengths 2r. Their other endpoints have separation at least 2r, so their included angle is at least pi/3.

For two rays at angle alpha, a point lying within distance a of both rays has distance at most a/sin(alpha/2) from their common endpoint. One proof is to minimize the larger of the two ray distances on each sphere centered at that endpoint. The minimum occurs on the internal angular bisector and is the sphere radius times sin(alpha/2). Projection behind a ray endpoint only increases its distance.

Consequently, incident rod overlap is confined within radius 2a of the shared endpoint. If a<=r/2, the rods merge only inside their shared physical junction sphere.

## Simultaneous certificate

For every globally nonoverlapping equal-sphere contact framework, choosing

    0 < a/r < 1/2

ensures all of the following:

- Equal-size-probe routes acquire no extra obstruction from the rods.
- Rods do not penetrate unrelated spheres.
- Rods with no common endpoint remain disjoint.
- Overlap between incident rods stays inside their shared junction sphere.

For r=1/sqrt(2), the common radius bound is a<0.3535533905932738. For example, a=r/4=0.1767766952966369 satisfies the certificate with slack.

All prescribed bars in `nonuniform_motion.py` retain length sqrt(2)=2r, and its certificate establishes global sphere nonoverlap for every permitted separated sequence in [40,50] degrees. The connector theorem therefore extends the entire infinite configuration family, including its time-parameterized kinematic examples. No extra finite spatial or temporal sampling is required.

## Composite solid fraction remains matched

For the periodic controlled comparison with a/r<1/2, the added rod material outside its endpoint spheres is disjoint from unrelated spheres and from all other added rod material. Each periodic motif cell contains six spheres, twelve internal octahedron rods, and three periodic bridges. Every rod has the same length 2r and radius a. Its added volume is therefore the same in every configuration.

More explicitly, put h=r-sqrt(r^2-a^2). Integrating the cylindrical neck outside the endpoint spheres gives

    V_neck=2*pi*(a^2*h-r*h^2+h^3/3)
          =2*pi*r*h^2-(4*pi/3)*h^3.

Thus the composite solid volume per cell is exactly

    V_solid=8*pi*r^3+15*V_neck.

Consequently, matching cell volume also matches the solid fraction after these finite rods are added. The controlled comparison is preserved for the composite obstacle, rather than merely for its spherical part.

For a/r<=1/2, probe exclusion redundancy holds for every q/r>=1/4. This follows by squaring the nonnegative condition q+a<=sqrt(q(2r+q)). At a=r/2 it reduces to q/r>=1/4. In particular, every probe ratio in [199/200,201/200] satisfies this condition.

## Improved numerical transport witnesses

`connectors_routes.py` starts from the independently constructed skew routes at theta=46.1113417712 degrees. It simplifies them with line-of-sight shortcuts and then optimizes interior waypoints under explicit segment-to-site distance constraints. A final segment projection check supplies the reported minimum distance. The optimization is local and does not establish maximal possible clearance.

| Winding | Segments | Minimum site distance | Margin above sqrt(2) |
|---|---:|---:|---:|
| (1,0,0) | 1 | 1.4179630639981604 | 0.0037495016250653 |
| (0,1,0) | 9 | 1.4188776328096897 | 0.0046640704365946 |
| (0,0,1) | 8 | 1.4188776328095287 | 0.0046640704364336 |

All three routes share a basepoint and supply a connected rank-three periodic network. The worst route margin exceeds 0.00374, improving the prior grid witness minimum of approximately 0.000002684. The fixed common basepoint limits the other two route margins.

For a fixed route, nearest-site distance is 1-Lipschitz under a uniform bound on site displacement. Thus a conservative total budget of 0.0035 is available for the sum of maximum obstacle-center displacement, obstacle-radius increase, probe-radius increase, and route-point displacement, provided the periodic route and the omitted-image control are preserved.

This is a numerical robustness bound for the saved paths. It is not a claim that arbitrary perturbations preserve contact bars or the original transport upper bound.

## Scope and evidence

The exact rod theorem applies to contact-center capsules. It excludes solid triangular faces, shells spanning an octahedron, external joint housings, and offset connectors. A solid octahedron would fill regions that these spheres and rods leave open.

The rod theorem does not establish hinge articulation, wall thickness, stress, bending, friction, actuation, stability, manufacture, or autonomous wave propagation. Full solid spheres touching along an edge also do not guarantee that an embedded physical hinge can articulate.

The geometric bounds are deductions here, with elementary written proofs and symbolic identity checks in `connectors_exact.py`. No priority claim is made.

The route coordinates and projections use floating-point arithmetic. They have not been converted to interval certificates. The positive result supports continued investigation of physical realization, rather than establishing a finished material design.

Reproduce:

    python3 research/connectors_exact.py
    python3 research/connectors_routes.py

Evidence: `connectors_exact.json`, `connectors_routes.json`.
