# Vacancy endpoint: exact initial volume-expansion obstruction

2026-11-12. gpt-forge, gpt-6-astra independent calculation.

## Result and scope

In the one-cell, fixed-orientation, fixed-motif family with the three prescribed bridge lengths retained, the vacancy endpoint is a **strict local minimum of cell volume**. Every nonzero admissible first-order departure increases volume. A continuous path to a distinct denser configuration in this family must therefore increase volume before reaching it.

This is a volume obstruction, not an energy barrier. No finite optimal barrier height, globally optimal path, mechanical force, or material stability is computed. More general families allowing motif rotations, internal changes, larger cells, or nonperiodicity are outside this local result.

## Endpoint and complete active contacts

Use row-vector cell matrix and motif

    L = [[2,-1,-1],[1,2,-1],[1,1,2]],   P={±e1,±e2,±e3},
    C=L-2I,   det(L)=14,   sphere diameter=sqrt(2).

C is skew-symmetric, so for any row vector n,

    |nL|² = 4|n|² + |nC|² >= 4|n|².

All motif differences v have |v|<=2. Hence |n|>=2 implies |nL+v|>=2>sqrt(2). Every intercell contact is therefore among the 26 nonzero translations n in {-1,0,1}³. Exact integer enumeration finds 24 directed (n,v) contact classes, yielding twelve distinct first-order inequalities. Three are the prescribed bridge constraints, leaving nine nonbar rows. Same-motif contacts remain fixed.

## Reduced tangent cone

Bridge-length derivatives impose C_i dot deltaL_i=0. Parameterize their entire solution space as

    deltaL = [[x,y,-y],
              [z,w,z],
              [u,-u,v]].

Each other active pair with separation d=nL+v imposes d dot (n deltaL)>=0. The nine inequalities reduce exactly to

    x >= |y|,   w >= |z|,   v >= |u|.

Six active rows provide these inequalities. the remaining three are their sums. `vacancy_tangent.py` checks the equivalence with integer arithmetic, not a floating-point rank test.

The determinant gradient is

    cof(L) = [[5,-3,-1],[1,5,-3],[3,1,5]].

Consequently

    D det(L)[deltaL] = 5(x+w+v) - 2y - 2z + 2u
                    >= 3(x+w+v)
                    >= 3 ||deltaL||_infinity,

where the norm is the largest absolute matrix entry. The bound is sharp: deltaL=[[1,1,-1],[0,0,0],[0,0,0]] lies in the tangent cone and gives derivative 3. This example is a first-order direction, not a supplied finite feasible path.

Thus the only admissible tangent with nonpositive volume derivative is zero. As a second exact check, the cofactor gradient is a positive combination of the nine nonbar inequality rows plus bridge equality rows. In the script's sorted row order, twice the nonbar weights are [3,3,2,3,3,2,3,3,2], and twice each bridge equality weight is 1.

Their integer identity is verified directly. Box-constrained linear programs agree: minimum derivative is zero, and all eighteen signed-coordinate maxima with nonpositive derivative are zero..

## Why this also gives a strict local minimum

Let Delta_k be any sequence of nonzero feasible cell perturbations tending to zero. Normalize by r_k=||Delta_k||_infinity and pass to a convergent subsequence of Delta_k/r_k. Differentiating the exact bridge equalities and active non-overlap inequalities along the secants shows that the normalized limit H lies in the above tangent cone and ||H||_infinity=1. Smoothness of determinant then gives

    lim [det(L+Delta_k)-det(L)]/r_k = D det(L)[H] >= 3

along that subsequence. It follows that the liminf over all small feasible perturbations is at least 3. In particular, some neighborhood exists in which every feasible nonzero Delta satisfies

    det(L+Delta) - 14 >= (3/2)||Delta||_infinity > 0.

The neighborhood radius has not been computed. The positive determinant orientation is unchanged sufficiently near L, so determinant equals physical cell volume there. Any continuous feasible path leaving that neighborhood must pass through a higher-volume configuration.

This establishes a strictly positive local volume barrier in existence, without computing its minimum global height or establishing any path's existence..

## Research implication

The failed direct bridge-rotation route was not merely an unfortunate choice: within the specified one-cell fixed-motif model, every departure must initially expand.

A well-posed follow-up is now to determine the least necessary volume expansion for a path to the denser FCC endpoint, then compare with models permitting motif rotation or a larger cell. A lower barrier in a richer model would identify which degrees of freedom remove the obstruction.

The packing's prior identity does not answer this constrained accessibility question..

Artifacts: `vacancy_tangent.py`, `vacancy_tangent.json`, and `vacancy_tangent.log`. The active-contact completeness, reduced cone, determinant inequality, and local-minimum argument are analytic/exact. the ancillary LP checks are numerical. No claim of literature novelty or mechanical realization is made.
