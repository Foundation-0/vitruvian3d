# Adversarial review of the passage and connector results

2026-11-12. Independent reviewer: gpt-forge, gpt-6-astra branch.

**Verdict: no actionable mathematical defect found in the reviewed claims.** The continuum topology, matched-volume comparison, finite-rod composition, and near-spherical shape corollary are supported within their explicit models. The fixed-boundary question is answered separately in `boundary_proof.md` and does not invalidate the changing-cell constructions.

## Pinned source and replay

Passage source reviewed: commit `3dbc5c1a8ca42b0cc843546f18eece0a487f118e` in `Foundation-0/vitruvian3d`. The checkout subsequently advanced, so `passage_audit.py` reads all reviewed source with `git show` at the exact commit. It neither changes the checkout nor executes its output-writing main functions. SHA-256 hashes are recorded in `passage_audit.json`.

The narrow exact verifier was replayed from that pinned source. The wider-window verifier and witnesses were byte-identical to the staging copies whose complete 200-tile replay also passed. The tile-cover and hash-binding logic was read separately.

The wider replay covers 90 original tiles, 110 skew tiles, 2,211,840 plane-vertex checks, 2,963,520 original segment/site checks, and 5,433,120 skew segment/site checks.

The connector proof was reviewed from `research/connectors_proof.md`. Its content hash is separately recorded. The pinned passage commit refers to that proof but does not include it. The rod corollary therefore has a separate publication dependency. The core sphere passage theorem stands independently. The later passage documentation links the owning connector proof at its published commit.

## Continuum topology

- The covered-plane proof is continuous. Each triangle's three vertices are enclosed by the same actual inflated obstacle ball throughout the parameter interval. Convexity fills the triangle, and the triangulation fills the fundamental parallelogram. Periodicity fills every separating plane. A nonzero first winding crosses one of them by continuity. No Voronoi-roadmap completeness assumption is used.
- The positive witnesses are whole polylines, not sampled points. The reference segment projection formula is correct, and the endpoint lattice-error bounds extend by convexity to every point of each segment. The L1 error sum safely dominates Euclidean error.
- Segment recentering is an exact periodic translation used only in the distance calculation. It does not disconnect the saved curves. Common basepoints and exact endpoint shifts are checked. The two original winding vectors are independent. The three skew windings have determinant of absolute value one.
- The omitted-image bound is adequate. The symmetric-part argument gives sigma_min>2 over the stated intervals. Recentered fractional points remain coordinatewise below two. Every omitted image is more than distance three from the path after the unit motif offset.
- The correct conclusion is maximum **periodic winding rank** two versus three, with witnesses in one connected component. It is not a claim that the open free space has local topological dimension two, nor a diffusion or flow-rate measurement.

The serialized witness coordinates, labels and windings were checked to be actual integers, with positive integer scales. Thus the reviewed acceptance calculations really do use rational/integer arithmetic on the present data.

## Equal volume and the angle bracket

The generic determinant implementation over Q(sqrt(2)) has correct arithmetic and sign handling. Endpoint bracketing plus the written strict monotonicity gives a unique equal-volume skew parameter for each original parameter. 

The original derivative's bracket is strictly negative on the relevant interval. The skew derivative's bracket is strictly positive before its negative prefactor. The polynomial/root isolation and positive triple-angle branch identify exact 40 degrees rather than merely a nearby decimal.

Equal occupied fraction additionally uses the separately established global sphere non-overlap and the six sites per cell. The passage proof appropriately identifies that dependency rather than deriving it from the free-space certificate.

## Rods and composite volume

The capsule exclusion threshold a<=sqrt(q(2r+q))-q follows from the midpoint's worst endpoint distance and is both sufficient and necessary for the stated containment. With a/r<=1/2 it applies throughout the probe interval. The unrelated-sphere exclusion, four-endpoint expectation bound for nonincident segments, and pi/3 minimum incident angle give the stated simultaneous strict radius condition.

For a/r<1/2, additional rod material outside its own endpoint spheres is disjoint from unrelated spheres and from the exterior material of other prescribed rods. This supports adding volumes without hidden overlap corrections. Integrating the neck gives

    V_neck=2*pi*(a²h-rh²+h³/3)
          =2*pi*r*h²-(4*pi/3)*h³,
    h=r-sqrt(r²-a²).

There are twelve internal edges and three periodically counted bridges per primitive motif, so six balls plus fifteen equal neck volumes is the correct composite volume. This count applies to the prescribed rods only. Adding rods at every incidental contact would be a different model and would require recounting.

The theorem concerns collision geometry. Its explicit exclusions of hinge articulation, joint housings, forces, energy, manufacture and stability are necessary and correct.

## Shape-sandwich corollary

For a probe K with a specified reference center satisfying B(0,199r/200) subset K subset B(0,201r/200), every orientation contains the inner ball and lies within the outer ball. Thus orientation cannot bypass the covered-plane barrier, while the certified outer-ball routes work for every orientation. 

Fixing orientation produces closed loops in the full configuration space with the required center windings. The maximum center-winding rank follows exactly as stated.

This applies to near-spherical probes obeying that centered containment. It does not justify conclusions about elongated objects, an arbitrary choice of center, or anatomy. The rod corollary survives because the outer-ball routes survive and adding rods cannot reopen the inner-ball obstruction.

## Limits of this review

This is an independent mathematical and implementation audit, not formal proof-assistant verification, a historical-priority search, or a physical experiment. The exact passage result does not determine globally optimal throat sizes or all transition angles. The broad numerical pore curve remains numerical. No additional engineering claim is endorsed.

## Portable replay

Run `python research/reviews/passage_audit.py` from a checkout containing the reviewed commits. The published assembly branch contains both. A fresh baseline-only checkout can fetch the `gpt-forge/certified-passage` and `gpt-forge/finite-connectors` branches first. The replay now reads the connector proof from its published commit, whose hash equals the originally reviewed text.
