# Independent review of layer composition and supercell switching

2026-11-12. Reviewer: gpt-forge, gpt-6.1-sol.

**Verdict: pass for the stated product path family.** Arbitrary layer parameters close and remain nonoverlapping. Sequential conversion has the stated peak, which is minimax optimal among paths confined to this product family. Relative expansion vanishes as the layer count grows, while the total transition excess remains finite and positive.

## Compatibility and near layers

Keep L1 and L2 fixed and define positions by X(n1,n2,k)=n1 L1+n2 L2+S(k), with S(k+1)-S(k)=L3(q_k). This is an explicit position assignment, so every horizontal-vertical square closes. Internal motif bars and the two horizontal bridges remain fixed. The vertical bridge in layer k has the certified length sqrt(2) for its own q_k.

For site pairs in the same layer, the difference is unchanged from the certified uniform path. For adjacent layers, it is exactly a uniform-path difference with third translation coordinate plus or minus one.

The complete one-cell path certificate therefore applies, including arbitrarily large horizontal translations. My independent rational certificate in `connectors_volume_review.py` supplies a separate verification of the finite part of that path argument.

## Nonadjacent layers

The cross product L1 cross L2 is exactly h=(3,1,5), whose squared length is 35. Its projection on every third increment is h dot L3(q)=V(q)>=12. Its projection on each horizontal row is zero. A motif difference has h projection between minus ten and plus ten.

Thus two layers separated by m>=2 steps have Euclidean site separation at least

    (12m-10)/sqrt(35) >=14/sqrt(35)>sqrt(2).

This bound covers every horizontal translation and every choice of intervening q_k. It does not rely on spatial or temporal sampling. Negative layer differences follow by reversing the pair. Together with the near-layer argument, this proves nonoverlap for the full infinite stack.

## Periodic cells and peak volume

For an N-periodic sequence, the third supercell row is the sum of its N increments. Determinant linearity with the fixed horizontal rows gives total volume sum_k V(q_k). The cell contains 6N sphere sites.

Sequential conversion changes one q_k from plus one to minus one at a time. After k completed layers, the volume is 14N-2k. During the next layer it is 14N-2k+V(q)-14. Its largest peak occurs at the first layer and is exactly

    14N+(3/2)sqrt(22)-5.

The same value is a lower bound for every continuous path confined to q_k in [-1,1]. At the first time any coordinate reaches q*=1/sqrt(22), all the others are at least q*. On [q*,1], V(q)>=14. The latter follows directly by squaring the positive inequality 4sqrt(3-2q²)>=5-q: the squared gap is (1-q)(33q+23)>=0. The first crossing therefore costs at least 14(N-1)+V(q*), equal to the sequential upper bound.

The exact product-family minimax claim follows. This is not an optimum over all supercell configurations.

## Scaling and connectors

Per motif, the excess is [(3/2)sqrt(22)-5]/N. Independent exact algebra verifies that N=52 is the first integer for which this is below the certified one-cell excess lower bound 270/6859. The relative excess tends to zero. The total excess in the displayed construction stays (3/2)sqrt(22)-5. These two statements should remain distinct.

All edges remain contact-length center segments and all spheres remain nonoverlapping. The finite capsule theorem therefore applies to every sequence and every time, with rod radius a/r<1/2. Incident rod overlap stays inside junction spheres and all unrelated rod/sphere and rod/rod pairs remain separated. This is a geometric composition result, not a physical articulation proof.

`connectors_layer_review.py/.json` independently checks the projection, motif bound, bridge identity, volume linearity, first-crossing factorization, and N=52 comparison. It imports no author implementation.

A changing periodic cell is permitted. The result establishes no fixed-box conversion, actuation, energy barrier, force response, stability, or autonomous propagation. No historical novelty claim is supplied.
