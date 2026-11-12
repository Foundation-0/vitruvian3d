# Matched-volume configurations with different dimensions of accessible motion

2026-11-12. gpt-forge, gpt-6.1-sol. Written argument with exact finite arithmetic certificates. Not formalized in a proof assistant. No claim of historical novelty.

## Statement

Let the motif be $P=\{\pm e_1,\pm e_2,\pm e_3\}$. Put physical spheres of radius $r=1/\sqrt2$ at every point of $P+\mathbb Z^3L$. Write

$$
L(s)=d(s)I+t(s)C,\qquad
d(s)=2+\sqrt2\frac{1-s^2}{1+s^2},\quad
t(s)=\frac{2s}{1+s^2}.
$$

Thus $s=\tan(\theta/2)$. Compare

$$
C_o=\begin{pmatrix}0&-1&-1\\1&0&-1\\1&-1&0\end{pmatrix},
\qquad C_v=\begin{pmatrix}0&1&1\\-1&0&1\\-1&-1&0\end{pmatrix}.
$$

The original interval is

$$
I_o=[363969/10^6,363972/10^6],
$$

and the skew interval is

$$
I_v=[425610/10^6,425635/10^6].
$$

For every spherical probe radius satisfying

$$
199/200\ \le q/r\le\ 201/200,
$$

the following statements hold in the continuum, not merely on a graph or a spatial sample:

1. For every $s\in I_o$, the maximum periodic winding rank of a connected free-space component is **exactly two**.
2. For every $s\in I_v$, there is a connected free-space component with winding rank **three**.
3. For every original parameter in $I_o$, there is a unique skew parameter in $I_v$ with exactly the same cell volume. Both configurations therefore have the same physical sphere packing fraction.
4. The original interval contains the exact angle $40^\circ$. Its equal-volume skew angle is numerically about $46.11134177^\circ$. the theorem uses the rigorous interval, not that decimal approximation.

Here a free path is a path of the probe's center outside the obstacle spheres inflated to radius $r+q$. Winding rank is the rank of the image in $\mathbb Z^3$ of loops within **one connected component** of the periodic free space. Directions in disconnected components are not pooled. The statement guarantees paths. it does not quantify flow rates, diffusion, mechanical work, or bodily motion.

The original and skew sphere frameworks are globally nonoverlapping along these paths. The existing exact folding and sign-pattern certificates establish that separate condition. The passage proof below also makes sense directly as a statement about the complement of the specified obstacle union.

## A covered plane rules out every bypass on the original route

Use fractional coordinates $x=fL(s)$. Consider the entire plane $f_1=0$, with fundamental parallelogram $0\le f_2,f_3\le1$. Divide each coordinate into 32 intervals, and split each resulting square along its increasing diagonal. This gives 2,048 triangles.

For each triangle, the saved witness chooses **one actual sphere center** $p+nL$. Only eight distinct obstacle spheres are needed for the whole parallelogram. The verifier proves that all three vertices of that triangle lie strictly inside that same inflated sphere, uniformly for every real $s\in I_o$. A Euclidean ball is convex, so the entire triangle lies inside that sphere.

Repetition by $L_2,L_3$ covers the whole plane. Translation by $L_1$ covers every plane $f_1=k\in\mathbb Z$.

The largest certified squared vertex-to-center distance is at most

$$
\frac{202766238443119}{102400000000000}
<\frac{399^2}{80000}=(r+(199/200)r)^2.
$$

The corresponding distance upper bound is about 1.40717414. Even the smallest permitted probe has exclusion distance about 1.41067803. The cover therefore has strict slack.

Any lifted continuous path with nonzero first winding must cross one of these planes by the intermediate value theorem applied to its first fractional coordinate. Every point of every such plane is forbidden. Consequently **no path anywhere in the original continuum can have nonzero first winding**. This is a direct obstruction. it does not assume that a finite Voronoi graph captures all possible bypasses.

## Explicit polylines prove positive access

The saved witness contains fractional-coordinate polylines with rational coordinates, denominator $10^6$. Multiplying these coordinates by (L(s)) defines a continuous family of actual spatial polylines throughout the selected parameter interval.

The original route has two loops with a common basepoint and windings

$$
(0,0,-1),\quad(0,1,0).
$$

They contain 16 nontrivial segments in total. Their winding vectors are independent and span the plane allowed by the obstruction above. This proves winding rank at least two, and hence exactly two.

The skew route has three loops with a common basepoint and windings

$$
(0,0,-1),\quad(-1,0,1),\quad(0,-1,0).
$$

They contain 24 nontrivial segments. The determinant of the three winding vectors has absolute value one. Thus they generate the full periodic lattice and prove rank three within one component.

All original and skew polylines are checked against the largest permitted probe. The calculation uses the rational upper exclusion distance

$$
B=5671/4000=1.41775,
\qquad B^2>\frac{401^2}{80000}=(r+(201/200)r)^2.
$$

The paths therefore remain free for all smaller probes in the stated interval.

### Exact distance checks over full segments and parameter intervals

The interval bounds use rational lower and upper bounds for $\sqrt2$:

$$
\frac{141421356237309504880}{10^{20}}
<\sqrt2<
\frac{141421356237309504881}{10^{20}}.
$$

Both squared inequalities are checked exactly. Monotonicity of (t(s)) and of $(1-s^2)/(1+s^2)$ gives rational interval bounds for every entry of (L(s)). The verifier rounds these bounds outward to multiples of $10^{-6}$.

For a fixed path segment and sphere center, subtract the center from both segment endpoints. Let $L_c$ be the midpoint of the entrywise interval matrix. The reference segment obtained using $L_c$ has rational endpoints. Its minimum squared distance to the origin is computed exactly by the clamped orthogonal projection formula:

$$
\min_{0\le u\le1}|a+uv|^2=
\begin{cases}
|a|^2,&a\cdot v\ge0,\\
|a+v|^2,&a\cdot v\le-|v|^2,\\
|a|^2-(a\cdot v)^2/|v|^2,&\text{otherwise}.
\end{cases}
$$

Entrywise lattice errors bound the displacement of both segment endpoints. Because a point of the segment is a convex combination of its endpoints, the same bound holds over the **entire** segment. Summing the three coordinate error bounds gives a rational upper bound $E$ on the Euclidean displacement. The verifier checks

$$
\text{reference minimum distance}>B+E
$$

by squaring and clearing positive denominators, using Python integers. No floating-point comparison enters acceptance. The reverse triangle inequality then proves that the actual segment stays outside every tested inflated sphere for every real $s$ in the interval and every real position along the segment.

There are 32,928 exact segment/site checks for the original loops and 49,392 for the skew loops. The baseline equal-size-probe calculations have conservative certified distance margins above $\sqrt2$ greater than 0.01142 and 0.00936 respectively. The stronger probe-size-window calculation is retained separately in the JSON.

### All omitted periodic images are excluded analytically

For each segment, apply an integer lattice translation so that all fractional endpoint coordinates have absolute value less than two. Convexity gives the same bound along the whole segment. This independent recentering is only a calculation device. The saved full polylines remain connected.

For $C_o$, the smallest eigenvalue of its symmetric part is $-1$. For $C_v$, it is zero. On both parameter intervals, the rational enclosure proves $d-t>2$. Cauchy--Schwarz therefore gives the uniform bound

$$
|wL|>2|w|\quad(w\ne0).
$$

If an integer image index lies outside $[-3,3]^3$, some coordinate has absolute value at least four. Since the path's corresponding fractional coordinate has absolute value less than two, $|f-n|>2$. As $|p|=1$,

$$
|(f-n)L-p|>4-1=3>B.
$$

Thus every omitted image is safe. The finite checks against $[-3,3]^3$ cover **all** remaining periodic spheres.

## The equal-volume pair is exact

The two determinants are

$$
D_o=d^3+dt^2+2t^3,\qquad D_v=d^3+3dt^2.
$$

Both decrease strictly on the relevant angle interval. For the original determinant, the canonical exact-motion proof gives, with $a=\sqrt2\cos\theta$ and $t=\sin\theta$,

$$
\frac{dD_o}{dt}=\frac{2t}{a}[-16-10a+3t^2+3at]<0.
$$

For the skew determinant, direct differentiation gives

$$
\frac{dD_v}{d\theta}=-3\sin\theta
[5\sqrt2+4\cos\theta-\sqrt2\cos^2\theta]<0.
$$

The independent verifier evaluates endpoint determinants using a generic $3\times3$ determinant over the exact field $\mathbb Q(\sqrt2)$. It proves

$$
D_v(425610/10^6)>D_o(363969/10^6),
$$

and

$$
D_o(363972/10^6)>D_v(425635/10^6).
$$

Strict monotonicity and continuity give a unique equal-volume skew parameter inside $I_v$ for every original parameter inside $I_o$. Both frameworks contain six equal spheres per cell, so equal cell volume gives equal sphere packing fraction.

For exact angle inclusion, let

$$
g(s)=s^6-33s^4+27s^2-3.
$$

The verifier proves $g(s_o^-)<0<g(s_o^+)$. Its derivative is $6s(s^4-22s^2+9)>0$ throughout a containing interval $0<s<0.37$. The triple-angle identity gives

$$
\left(\frac{3s-s^3}{1-3s^2}\right)^2=3.
$$

In the selected interval both numerator and denominator are positive, and $0<3\arctan s<\pi/2$. The isolated root therefore satisfies $3\arctan s=\pi/3$, hence $2\arctan s=40^\circ$. No decimal trigonometric evaluation is used to establish inclusion.

## Relation to finite rods

The separate elementary theorem in `connectors_proof.md` shows that a capsule of radius $a$ around a contact-center segment of length (2r) adds no probe-center exclusion when

$$
a\le\sqrt{q(2r+q)}-q.
$$

For $a/r\le1/2$, this condition holds whenever $q/r\ge1/4$: squaring reduces it to $(q/r+1/2)^2\le(q/r)(2+q/r)$. It therefore covers the entire probe interval above. The positive paths remain valid with these rods, and the existing covered plane remains a barrier after any rods are added. The stronger physical connector conditions in that separate proof apply when $a/r<1/2$.

This corollary concerns contact-center capsules. It does not establish working hinges, actuation, strength, friction, energy barriers, manufacture, or arbitrary additional joint housings.

For $a/r<1/2$, the separate connector bounds also make the added rod volumes disjoint outside the physical endpoint spheres. Each of the 15 rods per primitive cell has the same length and radius.

Its volume outside its two endpoint spheres is therefore constant. Six equal physical spheres plus 15 equal additional rod volumes give the same occupied solid volume per cell in either configuration. Equal cell volume consequently preserves equal occupied fraction for this sphere-and-capsule model too. This conclusion depends on the stated contact-rod geometry and the strict overlap bounds.

## Reproduction, controls, and evidence boundary

### A wider control window

`passage_window.py` extends the same theorem to the exact half-angle intervals

$$
I_o^{\rm wide}=[355/1000,364/1000],\qquad
I_v^{\rm wide}=[415/1000,426/1000].
$$

The original interval spans about $39.090^\circ$ to $40.003^\circ$. The skew interval spans about $45.077^\circ$ to $46.148^\circ$. These angle decimals are explanatory. The endpoints of the certified intervals are the rational half-angle parameters above.

The extension uses the same two original free loops and three skew free loops. A common $64\times64$ plane-cover witness replaces the smaller grid. Each whole interval is partitioned into closed parameter tiles of width $1/10000$: 90 original tiles and 110 skew tiles. Exact interval checks cover every real parameter in every tile. Adjacent tiles share their endpoints, so there are no unverified gaps.

The separate standard-library program `passage_window_verify.py` checks the saved witnesses through the independent verifier. It verifies interval coverage, the hash binding to the original loop witnesses, the covered plane on every original tile, the free loops on every tile, and the exact matched-volume bracket over the whole wider intervals. The earlier proof arguments then give rank exactly two versus rank three throughout this wider window, for every $q/r\in[199/200,201/200]$.

### A shape-sandwich corollary

Let a rigid probe $K$, referred to its chosen center, satisfy

$$
B(0,(199/200)r)\subseteq K\subseteq B(0,(201/200)r).
$$

Every orientation of $K$ contains the smaller concentric ball. The certified covered planes therefore block every nonzero first winding regardless of the probe's orientation. Every orientation of $K$ lies within the larger ball. The certified free polylines therefore work for every orientation, including arbitrary orientation changes while following those paths.

Holding orientation fixed gives closed loops in the full position-and-orientation configuration space with the already certified center windings. The maximum center-winding rank is consequently exactly two on the original route and three on the skew route for these probes too. This is a robustness result for shapes close to a sphere. It supplies no claim about elongated objects or anatomical bodies.

The contact-rod corollary also survives this sandwich. The outer-ball paths remain free with the specified rods, and adding rods cannot remove the inner-ball barrier.

### Commands and scope

Run:

```sh
python3 research/passage_certificate.py
python3 research/passage_verify.py
python3 research/passage_window.py
python3 research/passage_window_verify.py
```

The first program discovers convenient witnesses using NumPy/SciPy, then accepts them only through exact finite inequalities. It writes `passage_certificate.json`. The second program reads the saved witnesses and independently reconstructs every acceptance calculation with the Python standard library. It imports neither the discovery program nor NumPy/SciPy, and writes `passage_verify.json`.

The discovery program also rejects an insufficient plane-cover squared radius of 1.93 and rejects the proposed free loops when the exclusion radius is enlarged to 1.43. The independent verifier reconstructs the grid topology, requires the saved sphere labels to be real motif/image centers, checks common loop basepoints and exact endpoint shifts, and checks independent winding vectors.

The trusted mathematical ingredients are convexity of balls, continuity of lifted fractional coordinates, the projection formula, rational interval enclosure, the stated singular-value bound, and the determinant/angle arguments above. The arithmetic implementations have not been formalized. Reimplementation reduces shared-code risk. it does not replace review of the written argument.

The exact theorem establishes a robust change in the dimension of accessible motion at matched volume. It does not prove the numerically observed optimal pore thresholds or every transition angle. In particular, the numerical throat curve from `transport_probe.py` is not silently promoted to a theorem by this certificate. Novelty and practical value require separate literature and physical-design work.

Canonical original-motion proof: [continuous fold to fcc]$https://github.com/Foundation-0/incentre-lattices/blob/main/results/fold-continuous-fcc/PROOF.md$.
