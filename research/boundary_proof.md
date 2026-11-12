# Fixed-exterior rigidity and the limits of localized translation

2026-11-12. gpt-forge, gpt-6-astra independent derivation, with the same telescoping argument independently proposed by the lead.

## Exact no-go result

Start from a homogeneous lattice of identically oriented, undeformed motifs. Retain every designated bridge between neighboring motifs. If the exterior motifs are fixed, **no nonzero finitely supported motif translation can preserve all those bridge lengths**. The result holds for finite displacements of any size, not only infinitesimal motions, and does not require sphere non-overlap.

Likewise, for a periodic displacement field with the period cell held fixed, the only length-preserving motif translations are common translations of the whole framework. Fixing one motif removes that freedom.

These statements concern a homogeneous starting state and fixed motif orientations/anchor offsets. They do not establish rigidity of every nonuniform background or of models with motif rotations.

## Model and telescoping identity

Let motif origins initially be X(n)=nL, n in Z³. Let b_i be the constant initial vector of the designated bridge in lattice direction i. For the octahedral motif, b_i=L_i-2e_i. A motif translation u(n) changes that bridge vector to

    b_i + D_i u(n),  where D_i u(n)=u(n+e_i)-u(n).

Its squared-length change is exactly

    s_i(n)=2 b_i dot D_i u(n)+|D_i u(n)|².

If u has finite support, all nonzero summands lie on finitely many edges. For each i,

    sum_n D_i u(n)=0

by index translation. Therefore

    sum_{n,i} s_i(n) = sum_{n,i}|D_i u(n)|².             (1)

The right side is a discrete Dirichlet quadratic form. This name does not specify a physical elastic energy.

If every retained bridge length is unchanged, every s_i(n)=0. Equation (1) forces every D_i u(n)=0. Connectivity of Z³ makes u constant. Finite support makes it zero. More strongly, a nontrivial finitely supported field cannot keep every bridge no longer than its initial length: at least one must lengthen.

For a periodic field on a fixed finite quotient lattice, the same sum telescopes with periodic indexing. Equation (1) again implies a constant field. The baseline bridge vectors include the fixed period translations across the quotient boundary. Their displacement differences are periodic and hence telescope.

One may interpret the cancellation as a positive equilibrium-stress argument: opposite constant bridge directions cancel, while the sum of squared displacement differences is positive. The proof above is elementary and self-contained.

Related established rigidity theory includes [Connelly and Whiteley, *Second-Order Rigidity and Prestress Stability for Tensegrity Frameworks* (1996)](https://pi.math.cornell.edu/~connelly/pdf/10.1137_S0895480192229236.pdf). No novelty claim is made for this general mechanism or its application here.

## Approximate bridge preservation

Identity (1) remains exact when lengths are allowed to change. It quantifies an unavoidable total squared-length excess rather than merely detecting failure.

For a field supported in [-R,R]³, let U=max_n|u(n)|. Along either direction of each coordinate line through a maximizer, telescope from that point to the fixed zero exterior. If the two path lengths are ell_plus and ell_minus, their sum is 2R+2. Cauchy-Schwarz gives a contribution at least

    U²(1/ell_plus+1/ell_minus) >= 2U²/(R+1)

from that coordinate line. The edges in the three coordinate directions are disjoint, so

    sum_edges s = sum_edges|D u|² >= 6U²/(R+1).         (2)

If M>0 edges have nonzero displacement differences, at least one has squared-length increase at least 6U²/((R+1)M). This is a conservative bound, not an optimal strain-distribution calculation. Positive and negative individual squared-length changes may coexist. Their signed sum is the nonnegative quantity in (1).

## Why infinitesimal analysis gives a different answer at the vacancy state

The bar-only linearized equations are

    b_i dot D_i u(n)=0.

For finite-support u, following an i-coordinate line to the zero exterior gives b_i dot u(n)=0 at every site. Thus all finitely supported infinitesimal translation flexes are exactly the fields with values in the common kernel of the three b_i rows. If the bridge matrix is invertible, this kernel is zero.

At the vacancy endpoint,

    B=[[0,-1,-1],[1,0,-1],[1,1,0]],
    ker B=span{w}, w=(1,-1,1).

Every field u(n)=f(n)w with finitely supported scalar f is therefore a bar-only infinitesimal flex. For example, move only the motif at n=0 with infinitesimal velocity w. All six incident retained bridges have zero first-order length change.

Nevertheless, displacing that motif by t w increases each of those six squared lengths from 2 to 2+3t². Their total squared-length excess is 18t², exactly the Dirichlet sum. There is no finite length-preserving motion integrating this localized linearized flex. Nonlinear rigidity must not be mislabeled infinitesimal rigidity.

### Optional ideal axial-spring calculation

For the hypothetical bar-only field u=t f w, an affected bridge has length

    ell=sqrt(2+3t²(D_i f)²).

If one additionally posits independent ideal axial springs with energy (k/2)(ell-sqrt(2))², Taylor expansion gives the leading total energy

    (9k/16)t^4 sum_edges(D_i f)^4 + O(t^6)

for a fixed finite field f. This is a consequence of that explicitly added spring law, not a material-performance claim.

Such a field can violate nonbar sphere-contact inequalities already at first order.

It is **not established as a feasible hard-sphere packing motion or soft mode**. The exact fixed-boundary theorem needs no spring assumption.

## Consequence for the conversion-front construction

The layer-switching theorem does not contradict this result. Its periodic cell vector changes, and its infinite planar-front version changes far-field displacements. Those fields are not finitely supported translations with a fixed exterior. The affine or boundary contribution prevents the cancellation used in (1).

A finite block with all external motifs fixed and all crossing bridges retained is covered by the no-go theorem, by extending its displacement field by zero. Removing the crossing bridges would change the problem.

To seek a localized conversion from a homogeneous state one must relax at least one assumption, for example allow length changes, motif rotations/deformation, a changing boundary, or displacement extending into the far field. This proof does not guarantee that any particular relaxation succeeds.

## Replay and limits

`boundary_certificate.py` uses exact rational arithmetic to verify (1) on nontrivial finite fields, a fixed-period field, the bound (2), and the single-site infinitesimal counterexample. `boundary_certificate.json` records the outcomes. These examples replay the identities. The universal result rests on the written summation proof, not a finite sample extrapolation.

The theorem deliberately closes the specified fixed-boundary translation avenue. Further unspecified engineering models are not pursued here.

## Independent review

The independent review in `connectors_boundary_review.md` checks the summation argument and its scope. Its separate standard-library replay imports no author implementation.
