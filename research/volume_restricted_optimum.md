# Exact minimax optimum when only the third bridge moves

2026-11-12. gpt-forge, gpt-6-astra research branch.

## Statement

Fix L1=(2,-1,-1) and L2=(1,2,-1), retain motif P={±e_i} and sphere diameter sqrt(2), and let L3=(p,q,2+z) vary continuously while its prescribed bridge b=(p,q,z) has length sqrt(2). Among all non-overlapping continuous paths from b=(1,1,0) to b=(1,-1,0), the minimum possible maximum cell volume is exactly

    Vrestricted = 9+(3/2)sqrt(22) = 16.035623639735142... .

The path in `volume_path_proof.md` attains this value. This result allows arbitrary third-bridge routes, not just the earlier q-parametrization or a chosen plane. It does not prove optimality when the first two rows also move.

## A necessary spherical region

Single-row sphere contacts imply z>=|p|-1 and z>=|q|-1. These force z>=0: first z>=-1, and if -1<=z<0 then p²+q²<=2(1+z)² together with p²+q²+z²=2 would require 0<=z(4+3z), a contradiction.

Five actual pair constraints further require

    p-z<=1, q-z<=1, q-p<=1, p+q>=-1, q>=-1.

The exact contact witnesses and their squared gaps are checked by `volume_restricted_optimum.py`. Define Omega as the sphere p²+q²+z²=2 satisfying these five inequalities and z>=0. Every fully feasible path is contained in Omega. Extra contacts may shrink the feasible set. they cannot invalidate a lower bound proved for the larger Omega.

Omega is a spherical region with boundary cycle

    A -> C -> E -> D -> B -> A,
    A=(1,1,0), B=(1,-1,0), C=(1/3,4/3,1/3),
    D=(0,-1,1), E=(-1,0,1).

Here and below these are curved circular boundary arcs, not straight chords. Exact pair-plane/sphere intersections confirm these are all boundary vertices. The five arcs can also be derived directly:

| Arc | Equation and parameter interval | Volume function |
| --- | --- | --- |
| AB | p=1+z, q in [-1,1], z=(sqrt(3-2q²)-1)/2 | 9+q+4sqrt(3-2q²) |
| AC | q=1+z, p in [1/3,1], z=(sqrt(3-2p²)-1)/2 | 8+3p+3sqrt(3-2p²) |
| CE | q=p+1, p in [-1,1/3], z=sqrt(1-2p-2p²) | 11+4p+5sqrt(1-2p-2p²) |
| ED | q=-1-p, p in [-1,0], z=sqrt(1-2p-2p²) | 9+2p+5sqrt(1-2p-2p²) |
| DB | q=-1, p in [0,1], z=sqrt(1-p²) | 9+3p+5sqrt(1-p²) |

Substitution into the defining inequalities gives exactly these intervals. The z=0 boundary has only A and B, so it contributes no additional arc.

## Why an interior route cannot beat the boundary routes

Volume is the linear function V(b)=10+h dot b, with h=(3,1,5). Each boundary-arc volume above is concave in its displayed parameter, so its minimum is at an endpoint and is at least 12.

On the sphere, the only unconstrained critical points of V are ±sqrt(2)h/|h|. The negative one has z<0 and is excluded. the positive one is a maximum. Therefore the global minimum on Omega is on the boundary, and V>=12 throughout Omega. In particular h dot b>=2>0, so Omega lies strictly in the hemisphere facing h.

Let u=h/|h|. The point sqrt(2)u lies strictly inside Omega. Express any other point as

    b=sqrt(2)(cos(theta)u+sin(theta)w),
    w perpendicular to u, |w|=1, 0<theta<pi/2.

Every defining inequality has the form c+a dot b>=0 with c>=0 and a dot u>0. Along a fixed w, if a dot w>=0 the inequality stays positive throughout this angular interval. If a dot w<0 its left side decreases strictly with theta. Thus the allowed ray is an interval from u to a unique outer boundary angle theta_max(w). These boundary angles vary continuously with w.

Sending each point outward along its ray to theta_max defines a continuous retraction of Omega minus sqrt(2)u onto its boundary. It fixes boundary points and never increases V, because V=10+sqrt(70)cos(theta).

A path whose maximum is below the claimed optimum cannot pass through sqrt(2)u, where V=10+sqrt(70)>Vrestricted. It can therefore be retracted to a boundary path without increasing its peak. Consequently it suffices to compare the two boundary routes from A to B.

## Both boundary routes have the same bottleneck

Elementary one-variable maximization of the five concave expressions gives

    max_AB V = 9+(3/2)sqrt(22),
    max_AC V = 8+(9/2)sqrt(2),
    max_CE V = 9+(3/2)sqrt(22),
    max_ED V = 8+(9/2)sqrt(2),
    max_DB V = 9+sqrt(34).

The last three unequal values are smaller than Vrestricted. The direct boundary route AB therefore has peak Vrestricted, and the other route ACEDB also has peak Vrestricted because of CE. Every boundary path contains one of these two routes (possibly with backtracking). No path in Omega can have a smaller peak.

Finally, the explicit AB path is globally non-overlapping by the independently certified all-neighbor argument in `volume_path_proof.md`. Hence the lower bound for the necessary region is achieved by a truly feasible packing path, proving restricted minimax optimality.

## Bounded search beyond the restriction

`volume_general_search.py/.json` performs an exploratory sampled optimization allowing both first bridges to change their polar and azimuthal angles, in their initial off-diagonal sign sectors, while the third bridge moves on its upper sphere. At 41 fixed q values, with the reference path and continuation seeds, no sampled improvement was found: the optimizer returned the same configuration with the first two rows fixed to numerical precision.

This is not a global search and not an optimality argument for the full six-degree-of-freedom bridge family. It can miss disconnected routes, alternative sign sectors, and routes that do not progress through the chosen q interval in this way. Its limited result does not upgrade the restricted theorem.

## Remaining questions and tractability

- **All three bridges move, one cell:** six continuous shape degrees of freedom with polynomial length and non-overlap constraints. The exact restricted optimum supplies a strong candidate and test case. A rigorous global minimax lower bound would need a separating invariant or a certified connectivity calculation. a larger local optimization alone is unlikely to settle it.
- **Larger cells or nonperiodic changes:** substantially larger configuration spaces may permit a moving transition interface instead of a uniform switch. The one-cell lower bound does not apply directly. Start with a specified small supercell and retained contact graph. avoid treating arbitrary local phases as automatically compatible.
- **Motif rotations or deformation:** these change the model and may remove the obstruction. Useful only if the intended physical or geometric application permits them.

The restricted problem is now solved analytically. Further work should target one of these explicit model changes rather than continue resampling the same third-bridge family. No historical novelty, energy barrier, or material-performance claim is made.

## Independent review

GPT-6.1 Sol independently reviewed the contact constraints, every boundary interval, the absence of other boundary components, and the volume-decreasing retraction. The review found no issue. This is internal AI review, with outside mathematical review still open.
