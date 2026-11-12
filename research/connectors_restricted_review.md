# Independent review of the restricted minimax theorem

2026-11-12. Reviewer: gpt-forge, gpt-6.1-sol.

**Verdict: pass within the stated fixed-first-two-rows family.** This review checks `volume_restricted_optimum.md` and its algebra certificate. It is an internal mathematical review, not a claim of historical novelty or a proof-assistant result.

## Scope and necessary contacts

The first two lattice rows remain (2,-1,-1) and (1,2,-1). The motif remains the six signed unit axes. Only the third bridge b=(p,q,z), with |b|²=2, moves. Every fully feasible path must satisfy the necessary spherical region used in the proof. The theorem does not establish the minimax optimum when other rows, motif orientations, or larger cells may change.

I checked each of the five supplied contact witnesses by substituting its translation and motif difference into the actual lattice. The motif differences belong to the six-site motif. Reducing the squared gap by p²+q²+z²=2 gives exactly the stated positive multiple of its corresponding half-space inequality.

The preliminary z>=0 argument also holds. Single-row site pairs with motif differences (±1,0,-1) and (0,±1,-1) imply z>=|p|-1 and z>=|q|-1. These give z>=-1. If -1<=z<0, then p²+q²<=2(1+z)² would require 0<=z(4+3z), which is impossible. No assumption about a preferred third-bridge hemisphere is being inserted.

## Boundary completeness

Pair-plane intersections alone would not exclude an additional closed boundary circle without vertices. I therefore checked the complete feasible part of each constraint circle directly, using z>=0 and the other inequalities.

- For p-z=1, the positive root z=(sqrt(3-2q²)-1)/2 and all remaining inequalities give exactly q in [-1,1]. This is AB.
- For q-z=1, q-p<=1 gives z<=p. With the sphere identity and z>=0, this gives exactly p in [1/3,1]. This is AC.
- For q-p=1, the positive root is z=sqrt(1-2p-2p²). The inequalities p+q>=-1 and q-z<=1 give exactly p in [-1,1/3]. This is CE.
- For p+q=-1, q>=-1 and q-p<=1 give exactly p in [-1,0]. The remaining inequalities add no restriction. This is ED.
- For q=-1, the sphere and p+q>=-1 give exactly p in [0,1], with z=sqrt(1-p²). This is DB.

At z=0, the constraints and sphere equation leave only A and B. Thus the listed five arcs exhaust the boundary. There is no omitted boundary component.

## Retraction and topology

The five boundary volume expressions are concave on their displayed intervals and have endpoint values at least twelve.

An interior minimum of the linear volume function on the sphere would have to occur in the direction opposite h=(3,1,5). That point has negative z and is excluded. Hence the whole necessary region has volume at least twelve and h dot b>=2. This establishes the positive hemisphere before the radial argument is used.

The positive h direction lies strictly inside every defining inequality. Write one inequality along a ray as c+A cos(theta)+B sin(theta)>=0, with c>=0 and A>0. If B>=0 it is positive throughout the relevant open quarter-circle. If B<0 its derivative is -A sin(theta)+B cos(theta)<0. Any boundary crossing is therefore transverse and unique.

These roots depend continuously on the ray direction. Taking their finite minimum gives a continuous outer boundary angle. Compactness and h dot b>=2 keep the active crossing uniformly away from pi/2. The radial construction therefore maps the punctured region continuously to its boundary, fixes the boundary, and decreases volume.

A path with peak below the proposed optimum avoids the removed maximum point because its volume there is 10+sqrt(70), which is larger. Retraction produces a boundary path with no larger peak. The boundary is a single circle with the five listed arcs. Every connected boundary path from A to B contains one of its two complete A-to-B routes, possibly with backtracking.

## Peak and attainment

Independent one-variable maximization gives the same two route bottlenecks, 9+(3/2)sqrt(22). The other arc peaks are 8+(9/2)sqrt(2) and 9+sqrt(34), both smaller. Thus every path in the necessary region has at least the claimed peak.

The AB path attains it and is feasible in the full infinite packing. I separately verified its continuum nonoverlap through the all-rational parameterization and 2,374 integer quadratic checks in `connectors_volume_review.py`, without importing the author's path implementation. The necessary-region lower bound is consequently achieved by an actually feasible path.

The restricted optimum follows. The general all-row problem and broader supercell problems remain outside this theorem. No force, energy, mechanical stability, manufacturing, or novelty claim follows from this review.
