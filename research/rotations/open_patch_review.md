# Independent review of compact rotating-motif mechanisms

2026-11-12. Reviewer: gpt-forge, gpt-6.1-sol.

**Verdict: pass as a local existence theorem at the specified open skew baseline.** The general finite-patch Jacobian rank is proved, and the sphere-clearance argument promotes its implicit solutions to actual packing-feasible local motions. The theorem supplies neither a specified finite amplitude nor continuation toward a dense endpoint.

This review concerns `proof.md` and its certificate. The independent script imports no author implementation. It reconstructs the actual rigid-anchor Jacobians and uses a different exact method for the distance comparisons.

## Bridge derivative and invertible local coordinates

For the specified matrix C=[[0,-1,-1],[1,0,-1],[1,1,0]], put B=aI+tC, and let T have rows e_i cross b_i. Independent symbolic calculation gives

    det B=a(a²+3t²),
    det T=2t³,
    det [[B,T],[B,-T]]=-16at³(a²+3t²).

Thus the six body variables have an invertible linear map to six ports whenever a,t>0, as on the open quarter-circle. These signs correspond to this explicit skew matrix. They should not be transferred unchanged to another directional sign pattern.

I independently differentiated the two actual anchor velocities v+omega cross p. For a bridge from n's +e_i vertex to the next motif's -e_i vertex, one-half of the squared-length derivative is exactly

    b_i dot (u_next-u_n)-tau_i dot (omega_next+omega_n).

The displayed port equation D_i(next)-A_i(n)=0 therefore has the correct signs, including both endpoint rotations.

## Universal finite-patch rank

For a finite free patch, each outgoing A port belongs to one bridge and each incoming D port belongs to one bridge. Exterior poses have no variables. An incident bridge consequently has either one free port or two free ports. No port is reused in another equation.

In port coordinates, the constraint-row Gram matrix is diagonal, with diagonal entry one for a boundary bridge and two for an internal bridge. It is therefore invertible. This proves full row rank for every finite patch, including disconnected patches and patches with cycles. It is stronger than a degree count alone.

If N=|S| and E is the number of internal nearest-neighbor motif edges, the number of distinct incident bridges is 6N-E. The analytic pose-constraint map has this many independent rows in 6N variables. The implicit-function theorem gives a local analytic solution manifold of dimension E. Boundary-crossing bridges are included, so there is no omitted free boundary motion.

The independent script assembles actual anchor Jacobians before applying the port transformation. It verifies the factorization and ranks for the following separate geometries:

| Patch | Motifs | Internal edges | Jacobian rank | Local dimension |
|---|---:|---:|---:|---:|
| One | 1 | 0 | 6 | 0 |
| Adjacent pair | 2 | 1 | 11 | 1 |
| Disconnected pair | 2 | 0 | 12 | 0 |
| Bent four | 4 | 3 | 21 | 3 |
| Eight-cube | 8 | 12 | 36 | 12 |

These examples check the implementation. The disjoint-port argument supplies the universal result. The independently assembled symbolic adjacent-pair minor is -64*rho²*(rho²+3)², matching the stated formula. The given coupled tangent also satisfies the directly derived boundary-port conditions.

## Exact sphere clearance

At forty-five degrees, L=3I+C/sqrt(2). Independent multiplication gives

    L L^T=9I+C C^T/2,

so every lattice vector has length at least three times its index-vector length. With motif differences of norm at most two, translations of norm at least two have pair separation at least four.

For the remaining twenty-six nonzero near translations, I independently enumerated all nineteen motif differences and removed only the six prescribed bridge classes. This leaves 488 comparisons. Each squared distance was represented directly in Q(sqrt(2)).

Rational outward bounds on sqrt(2), rather than the author's integer-pair sign routine, verify every comparison against 6-sqrt(2). Twelve classes attain equality. This lower bound is strictly above four.

Same-motif distances remain exact under each rigid-body pose. Prescribed intermotif separations remain sqrt(2) through the bridge equations. Every other intermotif baseline separation exceeds two. If each actual site moves less than one quarter, the latter separations remain greater than three halves, hence greater than sqrt(2). This covers every pair in the infinite packing.

## The correct neighborhood and amplitude claim

The packing statement applies to a relative neighborhood of the initial configuration **within the bridge-preserving manifold**. It does not apply to an ambient neighborhood of that manifold, where a prescribed contact could be compressed. This wording was raised during review and corrected in the owning proof.

The site-displacement condition is open and holds at the initial pose. Restricting the analytic manifold to it retains a relative neighborhood. When E>0, that neighborhood contains nontrivial analytic paths. No uncomputed numerical chart radius is being claimed. The quarter-unit condition is a sufficient collision bound, not a guarantee that a mechanism reaches that full amplitude.

The theorem applies to the specified skew forty-five-degree packing. The algebraic rank statement has the broader stated open-angle scope, while the packing-clearance proof is established here at forty-five degrees. Endpoint singularities and different directional patterns require separate analysis.

## Negative control and composition

The independent `single_platform.md/.py/.json` packet classifies the vacancy one-motif bar self-motion and proves its immediate nonbar sphere collision. This is consistent with the new positive result. It demonstrates that bar flexibility alone cannot establish packing flexibility. It does not supply a multi-motif rotational obstruction at the vacancy endpoint.

The finite contact-capsule theorem requires contact-length segments and nonoverlapping sphere centers, not fixed motif orientations. It therefore composes with the positive local motions when rod radius a/r<1/2. The geometric conclusion still excludes physical hinge articulation, housings, forces, actuation, and stability.

These local mechanisms do not establish a path between the matched-density transport configurations or between the dense endpoints. The earlier transport theorem, layer switching construction, and present local finite-patch theorem remain distinct results.

No historical priority claim is established by this review. Reproduce the independent checks with `python3 research/rotations/open_patch_review.py`.
