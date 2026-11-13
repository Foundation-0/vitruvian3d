# Compact rotating-motif mechanisms in an open skew packing

2026-11-12. gpt-forge, gpt-6-astra branch. General port-variable argument developed jointly with the lead. Exact algebra and clearance replay supplied here.

## Main result

The translation-only fixed-boundary obstruction does **not** extend to rigid motifs that may also rotate.

At the homogeneous open skew configuration with angle 45 degrees, let S be any finite set of octahedral motifs that may translate and rotate, while every exterior motif remains exactly fixed. Retain all internal rigid motif geometry and all prescribed bridges, including those crossing the boundary of S. Let E(S) count nearest-neighbor lattice edges whose two endpoints both lie in S.

Then, near the initial configuration, the bridge-preserving rigid-body configurations form a real-analytic manifold of dimension **E(S)**. A relative neighborhood of the initial configuration within this manifold is globally sphere-nonoverlapping.

Consequently any finite patch containing an adjacent pair admits nontrivial, compactly supported, packing-feasible local motions with fixed exterior and unchanged bridge lengths.

For two adjacent motifs the local dimension is one. A single free motif, or any set with no internal nearest-neighbor edges, is locally rigid at this open configuration.

This is a local existence theorem, not a quantified large-motion continuation or a path between the dense endpoints. It applies to the explicit skew geometry below, not automatically to the original non-skew sign pattern. No physical hinge, actuation, energy, or historical novelty claim is made.

## Geometry and constraints

Let P={±e1,±e2,±e3}, sphere radius r=1/sqrt(2), and

    C=[[0,-1,-1],[1,0,-1],[1,1,0]],
    a=sqrt(2)cos(theta), t=sin(theta),
    B=a I+t C, L=2I+B.

Each row b_i of B has length sqrt(2). The baseline motif origin at n in Z³ is nL. A free motif has translation u(n) and orientation R(n) in SO(3), initially u=0 and R=I. Its six sites are nL+u(n)+R(n)P. Exterior motifs keep u=0 and R=I.

The bridge in direction i joins the +e_i site of n to the -e_i site of n+e_i. Its vector is

    g_i(n)=L_i+u(n+e_i)-u(n)-R(n+e_i)e_i-R(n)e_i.

The retained constraint is |g_i(n)|²=2. This definition includes every boundary-crossing bridge. Internal distances need no additional equations because each motif is represented by a rigid-body pose.

## Exact full rank for any finite patch

Use local angular-velocity coordinates omega(n), so the derivative of R(n)e_i is omega(n) cross e_i. The derivative of one-half the squared bridge constraint is

    b_i dot [u(n+e_i)-u(n)]
      - tau_i dot [omega(n+e_i)+omega(n)],
    tau_i=e_i cross b_i.

Let T have rows tau_i. Direct calculation gives

    det B=a(a²+3t²),  det T=2t³.

Both are nonzero for 0<theta<90 degrees. Thus the following local change of linear coordinates is invertible at every free motif:

    v_i(n)=b_i dot u(n),  w_i(n)=tau_i dot omega(n),
    A_i(n)=v_i(n)+w_i(n),  D_i(n)=v_i(n)-w_i(n).

In these port coordinates, each bridge equation becomes simply

    D_i(n+e_i)-A_i(n)=0.

Every port of a free motif appears in exactly one incident bridge equation. An internal edge pairs two free ports. A boundary edge sets one free port to zero because the other motif is fixed. Distinct rows therefore have disjoint nonempty column supports and are linearly independent.

There are 6|S| ports. Each internal edge reduces the number of distinct incident edges by one, giving 6|S|-E(S) independent constraints. The Jacobian has full row rank. The real-analytic implicit-function theorem then gives an actual local solution manifold of dimension

    6|S|-[6|S|-E(S)]=E(S).

This is not a Maxwell count extrapolated to a special geometry: the port transformation proves exact rank at the specified configuration. Since the exterior is fixed, the positive-dimensional motions cannot be dismissed as global rigid motions.

This bar-only rank result holds throughout the open skew interval 0<theta<90 degrees. The following clearance argument establishes packing feasibility at the concrete 45-degree baseline. It makes no automatic endpoint claim, where one of B or T becomes singular.

## A useful explicit two-body check

For S={0,e1}, scale the Jacobian by 1/t and write rho=a/t. Deleting its final angular-coordinate column gives an 11-by-11 minor

    -64 rho²(rho²+3)²,

which is nonzero for rho>0. One tangent normalization has

    u(0)=u(e1)
      =[2(rho²+1), -2(rho+1), -2(rho-1)]/[rho(rho²+3)],
    omega(0)=(1,1,-1), omega(e1)=(-1,-1,1).

The coupled rotation is essential. Setting all rotations to zero returns the prior fixed-exterior translation obstruction.

## Exact global clearance at 45 degrees

At theta=45 degrees, a=1, t=1/sqrt(2), and L=3I+C/sqrt(2). Since C is skew-symmetric,

    |nL|²=9|n|²+|nC|²/2 >=9|n|².

Every motif difference v has |v|<=2. For |n|>=2,

    |nL+v|>=3|n|-2>=4.

It remains to inspect the 26 nonzero n in {-1,0,1}³ and the nineteen distinct motif differences. Excluding the six designated bridge classes leaves 488 comparisons. Set u=3n+v and z=nC. Then

    |nL+v|² = |u|²+|z|²/2+sqrt(2)(u dot z).

The exact integer-pair certificate checks that every such value is at least 6-sqrt(2)>4. The bound is attained, for example, by n=e1 and v=(-1,1,0). Therefore every nondesignated intermotif pair has baseline distance strictly greater than 2. Same-motif distances are preserved by rigidity.

Restrict the local manifold to poses for which every actual sphere centre moves less than 1/4. This is an open neighborhood of the initial pose, and any nonconstant local analytic branch has nonzero points in it. A nondesignated pair then loses less than 1/2 of separation, leaving distance greater than 3/2>sqrt(2).

Internal contacts are unchanged. All designated intermotif contacts are held at sqrt(2) by the exact equations. Every pair in the infinite packing is therefore non-overlapping.

The 1/4 number is a geometric clearance condition on actual site motion, not a claimed parameter radius of the implicit-function chart. The theorem guarantees some nonzero local branch inside it. It does not certify a particular finite rotation amplitude.

## Vacancy negative control: a bar branch can still be blocked by spheres

At the vacancy endpoint a=0,t=1, B=C is singular and the port argument does not apply. A single motif with all neighbors fixed has an explicit bar-only branch. Let w=(1,-1,1), translate by u=tau*w, and rotate about w/|w| through angle -alpha. The six incident bridges retain their exact lengths when

    tau²=(4/9)(cos(alpha)-1)+(4/(3sqrt(3)))sin(alpha).

Near alpha=0 this gives alpha=O(tau²), so rotation first appears at second order in translation. It is a direct way for rotations to evade the translation-only sum-of-squares obstruction.

It is not a packing-feasible branch near the starting state. The moving +e2 sphere has a nondesignated contact with the fixed sphere at L3-e3=(1,1,1), whose initial separation is d=(1,0,1). Since d dot w=2, its squared distance along the branch is

    2-4*tau+O(tau²),

which overlaps for small tau>0. The opposite pair overlaps for small tau<0. This is a negative control demonstrating why a bar-only mechanism is not enough. No full classification of multi-motif vacancy mechanisms is claimed here.

## Connector composition and limits

The established contact-capsule theorem applies to these sufficiently small packing-feasible motions: every designated segment still has length 2r and all distinct sphere centres remain at least 2r apart. Capsules with radius below r/2 inherit the geometric nonintersection bounds. This does not establish that physical joints embedded in the spheres can articulate or that a chosen actuation realizes the motion.

The body-and-bar viewpoint is established. Relevant primary prior art is [Tiong-Seng Tay, *Rigidity of multi-graphs. I. Linking rigid bodies in n-space*, JCTB 36 (1984), 95–112](https://doi.org/10.1016/0095-8956(84)90016-9).

The result here is the explicit nonsingular port transformation and clearance certificate for this lattice, not the invention of rigid-body linkage counting. Its historical priority has not been established.

## Independent review and finite samples

[The independent positive review](open_patch_review.md) rebuilds the anchor Jacobians and verifies all 488 distance classes by a separate exact method. [The single-motif control](single_platform.md) supplies a complete local classification at the vacancy state.

[The numerical packet](numerics/README.md) records 81 distinct frames of the open-state pair, reaching about 7.11 degrees of rotation. Its exact positive-stress controls also establish local isolation of the specified one- and two-motif vacancy patches. [A separate Astra audit](numerics/audit_numeric_stress.md) authenticates those controls and reconstructs the frames with quaternion rotations.

These results distinguish a proved local continuous mechanism, sampled finite amplitudes, and exact dense-state obstructions. They do not establish a transition between the open and dense postures.

## Reproduction and stopping point

Run `python3 research/rotations/certificate.py`. It verifies det B, det T, the explicit two-body minor, port incidence for representative finite patches, all 488 exact clearance comparisons, and the single-vacancy bar-only branch identity. The universal rank and existence statements rest on the written port and implicit-function argument, not on a finite patch sample.

The specific open gap is resolved constructively: allowing rigid motif rotations permits compact local motion with fixed exterior in an open packing, despite translation-only rigidity. Further global continuation toward the dense vacancy endpoint would be a separate question. No unbounded search is needed to establish the present conclusion.
