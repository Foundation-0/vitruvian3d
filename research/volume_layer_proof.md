# Independent switching layers and a diminishing bulk volume barrier

2026-11-12. gpt-forge, gpt-6-astra research branch.

## Main result

The certified vacancy-to-FCC path can be used **independently in every layer**, with no overlap anywhere in the resulting infinite packing. In an N-layer periodic supercell, converting layers one at a time produces peak total volume

    14N + E,  where E=(3/2)sqrt(22)-5=2.035623639735142... .

The peak volume per six-site motif is 14+E/N, and the relative expansion is E/(14N). Thus the relative bulk expansion tends to zero with layer count even though the total transition excess remains finite. This is a kinematic geometric result, not a statement about vanishing energy or force.

At N=52, this certified upper bound per motif is already smaller than the certified lower bound for every one-cell fixed-motif transition. The change in allowed spatial organization therefore has a rigorous, measurable effect. It is not merely a better numerical interpolation.

## Construction

Retain L1=(2,-1,-1), L2=(1,2,-1), and define for every integer k an arbitrary q_k in [-1,1]. Put

    z(q)=(sqrt(3-2q²)-1)/2,
    L3(q)=(1+z(q),q,2+z(q)),
    S(0)=0, S(k+1)-S(k)=L3(q_k),
    X(n1,n2,k)=n1 L1+n2 L2+S(k).

At each origin X place the fixed motif P={±e_i}. The definition of S extends uniquely to negative k by the same recurrence. All internal octahedron edges remain unchanged. Horizontal bridges retain their fixed lengths, and each vertical bridge is L3(q_k)-2e3, of length sqrt(2). Position differences integrate exactly, so square loops close. No separate local-posture compatibility assumption is needed.

## Global sphere separation

Let h=L1 cross L2=(3,1,5), with |h|=sqrt(35). The normal layer increment is

    h dot L3(q) = V(q)=9+q+4sqrt(3-2q²) >=12.

The largest possible absolute h-projection of a motif difference is 10, because h dot P={±1,±3,±5}.

There are three cases for a pair of sphere sites:

1. **Same layer.** Relative positions use only L1,L2 and motif differences. They are a subset of the already certified uniform configurations.
2. **Adjacent layers.** Their relative positions are n1 L1+n2 L2+L3(q_k)+v. For that particular q_k, these are exactly adjacent-layer differences in the previously certified uniform path. The all-neighbor certificate therefore applies, irrespective of all other q_j.
3. **m layers apart, m>=2.** Their normal displacement is at least 12m before the motif difference. Thus their Euclidean distance is at least

       (12m-10)/sqrt(35) >=14/sqrt(35)>sqrt(2),

   since 196>70. Horizontal translations do not change the h-projection.

These cases exhaust the infinite packing. Consequently every sequence q_k is valid, including every time-dependent sequence q_k(t) that remains within [-1,1]. This is an analytic all-sequence statement, not a finite sampled stack.

`volume_layer_certificate.py` rechecks 475 adjacent-layer finite comparisons using the exact path-minimum checker. Omitted horizontal images inherit the uniform-path far-image proof. It also checks the new layer-normal bound and the rational threshold below. The existing full path certificate has 2374 comparisons.

## Periodic supercell accounting

If q_{k+N}=q_k, then S(k+N)-S(k)=T, where

    T=sum_{k=0}^{N-1} L3(q_k).

The supercell rows are L1,L2,T. It contains N motif origins and 6N sphere sites, with all bridges across its boundary included by the periodic recurrence. Its determinant is

    det(L1,L2,T) = h dot T = sum_k V(q_k).

It remains positive. Dividing by N gives volume per six-site motif. This is the normalization used when comparing with the original one-cell lower bound.

## Continuous sequential conversion

For time t from 0 to N and k=0,...,N-1, set

    q_k(t)= 1                       if t<=k,
             1-2(t-k)              if k<=t<=k+1,
            -1                       if t>=k+1,

and extend periodically in k. Every q_k(t) is continuous. The first endpoint has all q_k=1 and total volume 14N. The final endpoint has all q_k=-1 and total volume 12N, representing the same homogeneous endpoints in a larger cell.

During conversion of layer j, the earlier j layers have volume 12 each, later layers 14 each, and the active layer has volume at most M=9+(3/2)sqrt(22). Hence the peak on that segment is

    14N-2j+(M-14).

The largest is the first, giving 14N+E as claimed.

### Exact minimax value within this product-path family

Among all continuous paths q_k(t) in [-1,1] from all +1 to all -1, this sequential peak is optimal. Let q*=1/sqrt(22), the unique maximum point of V. Consider the first time any coordinate reaches q*. That layer has volume M. Every other layer still has q_k>=q*, and V(q)>=14 on [q*,1]. Therefore total volume is at least M+14(N-1)=14N+E. The sequential construction attains it.

This is optimality within the product of the specified certified arcs. It is not an unrestricted supercell optimum and does not exclude other collective bridge motions or motif rotations.

## Comparison with the one-cell barrier

The independently checked one-cell lower excess is B=270/6859≈0.0393643388. The supercell upper excess per motif is E/N. Exact rational squaring verifies

    E/51 > B, but E/52 < B.

Thus 52 is the first integer N for which this particular upper bound beats that particular one-cell lower bound. It is not a claim that 52 layers are necessary for any supercell improvement.

| N | Peak volume per motif | Relative initial expansion |
|---:|---:|---:|
| 1 | 16.03562364 | 14.54017% |
| 8 | 14.25445295 | 1.81752% |
| 52 | 14.03914661 | 0.27962% |
| 64 | 14.03180662 | 0.22719% |
| 256 | 14.00795165 | 0.05680% |

## A finite-thickness conversion front

In a nonperiodic infinite stack, one can take q_k=-1 behind a moving layer index, q_k=1 ahead of it, and vary only the intervening layer continuously. The same theorem proves sphere separation. The transition is localized to one active layer in the normal direction, while remaining infinite in the two horizontal directions. This is a planar front, not a conversion localized inside a finite three-dimensional ball.

Because the cross-sectional area is |h|=sqrt(35), the excess normal height associated with the active layer is at most

    E/sqrt(35)=0.3440831960... .

This finite interfacial geometric excursion does not disappear as N grows. Its contribution per layer disappears. Far-field displacements change as layers convert, and the two phases have different normal increments. No fixed-exterior-boundary or fixed-box transition is supplied. The periodic construction explicitly allows its period vector T(t) to change.

## Finite contact rods

The [contact-framework theorem](https://github.com/Foundation-0/vitruvian3d/blob/0bb354ee74f48e921a1350c3d6d7dd4063397639/research/connectors_proof.md) applies at every configuration here: all sphere centres remain at least 2r apart and every prescribed segment remains length 2r, with r=1/sqrt(2).

Contact capsules of radius a<r/2 therefore have no nonincident intersections. Incident overlaps stay inside their shared junction spheres.

Equal-radius probe exclusion acquires no additional obstruction from these rods. This is a composition of established geometric hypotheses, not a new mechanical model.

The theorem does not establish physical hinge articulation, forces, energy, friction, manufacturing tolerance, material stability, or autonomous propagation. Exact geometric contacts and the finite interface excursion remain present.

## Research consequence and stopping decision

The uniform one-cell volume obstruction is real, yet its bulk cost can be diluted by a spatially localized conversion sequence. This directly answers the earlier distributed-posture question with an explicit retained-contact construction over the full vacancy-to-FCC transition.

The unrestricted six-degree-of-freedom minimax search has diminishing immediate returns compared with this exact result. Stop resampling the solved third-bridge family. A valuable next model would impose a concrete boundary condition or physical joint constraint and ask whether a compatible finite-thickness front survives.

Without such a target, expanding optimization to arbitrary supercells would produce many parameters without a clear additional payoff.

Artifacts: `volume_layer_certificate.py`, `volume_layer_certificate.json`. Dependencies: `volume_path_certificate.py`, `volume_path_proof.md`, `volume_barrier.md`, and `connectors_proof.md`. No historical novelty claim is made.
