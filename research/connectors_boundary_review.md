# Independent review of the fixed-exterior theorem

2026-11-12. Reviewer: gpt-forge, gpt-6.1-sol.

**Verdict: pass within the stated translation-only model.** The exact theorem, approximate bound, periodic version, and infinitesimal counterexample are consistent. No author implementation was imported by the independent replay.

## Exact identity and scope

Each designated bridge changes from a constant vector b_i to b_i+D_i u. Its squared-length change is 2b_i dot D_i u+|D_i u|². For a finitely supported field, the sum of D_i u over all lattice sites is zero by index translation. The linear terms therefore vanish exactly, leaving the nonnegative sum of squared displacement differences.

If every bridge retains its length, that sum vanishes. Every difference is consequently zero. Lattice connectivity makes the displacement constant, and the fixed exterior makes it zero. The argument also rules out a nonzero finite field in which every retained bridge stays no longer than its initial length.

Periodic indexing gives the same cancellation when the cell is fixed. The conclusion is a common translation, with that freedom removed by pinning one motif. Bridges across the cell boundary still use their fixed baseline vectors, so there is no missing boundary term in this case.

This argument requires the baseline bridge vector to be constant in each lattice direction. A varying background, a changing cell or exterior, changing anchor offsets, or motif rotations are different models. It uses no sphere-separation assumption. The theorem concerns finite displacements, rather than infinitesimal rigidity.

## Approximate preservation bound

Choose a displacement maximizer in the support cube [-R,R]³. Along each of its six coordinate rays, telescope to the zero exterior. The two lengths in one coordinate direction are R+1-x and R+1+x. Their reciprocal sum obeys the exact identity

    1/(R+1-x)+1/(R+1+x)-2/(R+1)
      =2x²/[(R+1)((R+1)²-x²)] >=0.

Cauchy–Schwarz therefore gives at least 2U²/(R+1) on that coordinate line. The three pairs of rays have disjoint edges. Their sum gives 6U²/(R+1), as stated. For a nonzero field, the number M of affected edges is positive and the displayed per-edge increase bound follows by averaging the signed squared-length excesses.

## Linearized versus finite behavior

The bar-only linearized equations imply b_i dot u is constant along each i-direction line. Finite support makes each such constant zero. Thus the complete finite-support linearized translation space is exactly the pointwise common kernel of the three b_i rows.

For the vacancy bridge matrix, that kernel is spanned by (1,-1,1), whose squared norm is three. Moving one motif by t times this vector changes all six incident squared bridge lengths by exactly 3t². The total is 18t². The first-order changes vanish, while finite length preservation fails. The proof correctly avoids claiming infinitesimal rigidity.

The optional spring calculation is also algebraically correct: the leading extension is 3t²(Df)²/(2sqrt(2)). Squaring and multiplying by k/2 gives 9k t^4(Df)^4/16. This calculation requires its separately stated spring law. It does not make the bar-only field a feasible hard-sphere motion, and the manuscript explicitly preserves that distinction.

## Relation to layer conversion and replay

The layer construction changes the periodic vector or far-field displacement. Its affine boundary contribution prevents the telescoping cancellation used in the fixed-exterior proof. There is no contradiction with the layer front. No fixed-box front follows from the existing construction.

`connectors_boundary_review.py` independently evaluates a sparse rational displacement field with arbitrary rational baseline bridge vectors on a surrounding box. The signed squared-length sum equals its Dirichlet sum exactly, and the three directional linear sums vanish. It also checks the harmonic identity and vacancy-kernel arithmetic. These examples supplement the universal written summation argument.

The result closes the specified finite-support translation avenue. It establishes no energy barrier, physical actuator, rotational rigidity, or historical novelty.
