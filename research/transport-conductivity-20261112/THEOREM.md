# From certified passages to a static conductivity tensor

The existing passage certificates imply rank two versus rank three for a specified variational conductivity tensor. The same prescribed bar graph and occupied fraction therefore coexist with different ranks of this static transport operator. No mesh calculation, reflected Brownian limit or material experiment is needed for this statement.

This note uses the narrow passage intervals at source commit `e7445d56b9ecfd97b408653a2726f3bab70a6e4c`. The original interval contains the exact 40-degree state and has an exactly volume-matched skew partner. The result also applies throughout the certified probe-radius interval. Historical priority and human mathematical review remain open.

**Geometry and operator.** Use physical coordinates in the source normalization. The sphere radius is \(r=1/\sqrt2\), the motif is \(P=\{\pm e_1,\pm e_2,\pm e_3\}\), and row vectors of \(L\) generate the lattice \(\Lambda=\mathbb Z^3L\). For \(s=\tan(\theta/2)\),

\[
L=dI+tC,\quad
d=2+\sqrt2(1-s^2)/(1+s^2),\quad t=2s/(1+s^2).
\]

The original and skew matrices are

\[
C_o=\begin{pmatrix}0&-1&-1\\1&0&-1\\1&-1&0\end{pmatrix},
\qquad
C_v=\begin{pmatrix}0&1&1\\-1&0&1\\-1&-1&0\end{pmatrix}.
\]

Use \(s_o\in[363969,363972]/10^6\) and \(s_v\in[425610,425635]/10^6\). Let \(q/r\in[199/200,201/200]\). In the flat torus \(T=\mathbb R^3/\Lambda\), let \(D\) be the open set whose distance from every site \(p+nL\) exceeds \(r+q\). Its full lift is the allowed space of probe centers. Boundary tangency is excluded from this open conducting set.

Put constant scalar conductivity \(\sigma>0\) on \(D\) and zero conductivity outside. With \(V=\det L>0\), define

\[
\xi^{\mathsf T}K_D\xi
=\frac{\sigma}{V}
\inf_{\phi\in H^1(T)}
\int_D|\xi+\nabla\phi|^2\,dx. \tag{1}
\]

Here \(\phi\) is a single-valued periodic potential with square-integrable first derivatives. Normalization is by the whole cell volume \(V\), not the conducting volume \(|D|\). The prescribed uniform field \(\xi\) is a physical Cartesian vector. Thus matrices from different lattice bases are not compared as though fractional coordinates were orthonormal.

Equation (1) defines the operator under study. It conducts on probe-center accessible space, not on the original solid sphere-and-rod material. It supplies no claim that particle diffusion, fluid flow or mechanical response follows this equation.

No boundary smoothness or extension theorem is required for this variational definition. In \(L^2(D,\mathbb R^3)\), take the closure \(G\) of the restricted periodic gradients. The infimum is the squared distance of the constant field \(-\xi\) from \(G\). Orthogonal projection therefore makes (1) a finite symmetric positive semidefinite quadratic form, whether or not a minimizing potential exists.

For context, the uniformly elliptic periodic cell formula is given by [Allaire, Section 1.2.3, equation 1.16](https://www.cmap.polytechnique.fr/~allaire/homog/lect1.pdf). Its coercivity assumption does not hold at zero conductivity. There is a direct cell-level limit here: assign conductivity \(\varepsilon>0\) outside \(D\) and keep \(\sigma\) inside. For each fixed periodic potential its energy tends to (1).

The infimum for positive \(\varepsilon\) is bounded below by (1). Testing with any fixed potential arbitrarily close to the infimum in (1) gives the reverse upper bound as \(\varepsilon\downarrow0\). Thus the positive-contrast cell tensors converge to \(K_D\). This does not interchange a macroscopic homogenization limit with a contrast limit or establish a stochastic limit law.

**Flux bound from one open route.** Let a piecewise smooth lifted route \(\gamma:[0,1]\to\mathbb R^3\) have endpoint shift \(v=\gamma(1)-\gamma(0)\in\Lambda\). Its projection is a closed loop in \(D\). Suppose the periodic cube neighborhood \(\gamma+[-h,h]^3\) stays inside \(D\), and the kernel support does not overlap its nonzero lattice translates.

Define a product tent kernel of mass one:

\[
\rho_h(x)=\prod_{j=1}^{3}\frac{(1-|x_j|/h)_+}{h},
\qquad
\|\rho_h\|_2^2=\frac{8}{27h^3}.
\]

Periodize it over \(\Lambda\), and convolve it with the oriented route current:

\[
J_\gamma(x)=\int_0^1
\dot\gamma(u)\rho_h^{\rm per}(x-\gamma(u))\,du. \tag{2}
\]

The current is supported inside \(D\). In the distributional sense its divergence is the difference of the endpoint kernels, which vanishes because the endpoint shift is a period. Integrating (2) over one cell gives \(\int_TJ_\gamma\,dx=v\). Self-intersections of the route cause no difficulty. The construction does not require an embedded tubular coordinate system.

If the route length is \(\ell\), the triangle inequality for the \(L^2\) norm gives

\[
\|J_\gamma\|_2\le\ell\|\rho_h\|_2.
\]

Periodic integration by parts gives \(\int_T\nabla\phi\cdot J_\gamma=0\) for every \(\phi\in H^1(T)\). Cauchy-Schwarz, followed by the infimum in (1), yields

\[
\xi^{\mathsf T}K_D\xi
\ge\frac{\sigma}{V}
\frac{(\xi\cdot v)^2}{\|J_\gamma\|_2^2}
\ge\frac{27\sigma h^3}{8V\ell^2}(\xi\cdot v)^2. \tag{3}
\]

A collection of routes with independent physical shifts therefore supplies strictly positive response on their span. This follows from finite-energy fluxes, not from a line integral of an arbitrary Sobolev gradient along a zero-width curve.

**Source-bound route widths and lengths.** The existing exact verifier proves a distance greater than \(B=5671/4000\) from every route point to every obstacle site. It covers each entire segment, each parameter interval and every periodic image. For an equal-radius probe, the exclusion radius squared is 2. The choice \(h=1/1000\) satisfies

\[
(B-2h)^2>2.
\]

The kernel support lies within distance \(\sqrt3h<2h\) of the route, so it is wholly inside \(D\). For the entire probe-radius window, \(h=1/4000000\) instead satisfies

\[
(B-2h)^2>401^2/80000.
\]

The source also gives the smallest singular value of \(L\) greater than 2. Every nonzero lattice vector consequently has length greater than 2. Since \(2\sqrt3h<2\), the periodized kernel copies do not overlap, preserving the displayed kernel norm.

On both intervals, \(d<4\) and \(|t|<1\). Therefore \(\|L\|_{\rm op}<8\) and \(V<81\). For a route with fractional vertices \(f_j\),

\[
\ell\le\|L\|_{\rm op}\sum_j|f_{j+1}-f_j|
<8\sum_j|f_{j+1}-f_j|_1<20.
\]

The final inequality is checked with integers for all five saved routes. It is uniform over the narrow parameter intervals. These deliberately loose bounds are sufficient to certify positivity.

**Original configuration: a zero direction.** The plane certificate bounds every point on every plane \(f_1=k\) by squared distance

\[
Q=202766238443119/102400000000000<(176/125)^2
\]

from a genuine obstacle center. The same center covers every point within physical distance \(1/500\) of that plane because

\[
176/125+1/500=141/100,\qquad
(141/100)^2<399^2/80000.
\]

Thus an entire periodic slab lies outside \(D\), uniformly over the probe-radius interval. Let \(g=L^{-1}e_1=\nabla f_1\). Its norm exceeds \(1/8\) because \(Lg=e_1\) and \(\|L\|_{\rm op}<8\). The fractional strip \(\operatorname{dist}(f_1,\mathbb Z)<1/4000\) is consequently inside the insulating slab.

Choose a bounded periodic function \(b(t)\), supported in these fractional strips, with integral one over a period. Let \(H'=b\) and \(H(t+1)=H(t)+1\). The function

\[
\phi(x)=H(f_1(x))-f_1(x)
\]

is periodic and belongs to \(H^1(T)\). For \(\xi=g\), the field \(\xi+\nabla\phi=b(f_1)g\) vanishes throughout \(D\). Hence \(K_Dg=0\). This construction supplies the zero-energy potential explicitly.

The saved original routes have physical shifts \(-L_3\) and \(L_2\). They span \(S=\operatorname{span}\{L_2,L_3\}=g^\perp\). The two nonzero singular values of the matrix with these rows exceed 2, because its row Gram matrix is a principal submatrix of \(LL^{\mathsf T}>4I\).

Averaging their two bounds (3), and using \(\ell<20\) and \(V<81\), gives

\[
K_o\succeq\frac{\sigma h^3}{4800}P_S,\qquad K_og=0. \tag{4}
\]

Here \(P_S\) is Euclidean orthogonal projection onto \(S\). Therefore \(K_o\) has rank exactly two.

**Skew configuration: positive in every direction.** The three saved winding rows form

\[
W=\begin{pmatrix}0&0&-1\\-1&0&1\\0&-1&0\end{pmatrix}.
\]

Their physical shifts are the rows of \(WL\). The leading principal minors of \(W^{\mathsf T}W-I/3\) are \(2/3,4/9,2/27\), all positive. Hence

\[
\sum_j(\xi\cdot v_j)^2
=|WL\xi|^2>\frac43|\xi|^2.
\]

Averaging the three bounds (3) gives

\[
K_v\succeq\frac{\sigma h^3}{21600}I. \tag{5}
\]

Thus \(K_v\) is positive definite. At \(q=r\), (4) and (5) become \(K_o/\sigma\succeq P_S/4800000000000\) and \(K_v/\sigma\succeq I/21600000000000\). These are conservative witness bounds. They estimate neither an optimal conductivity nor useful engineering performance.

**Components, rods and normalization.** All saved routes in each configuration share a basepoint. Their cube neighborhoods lie in the same connected component of \(D\) on the torus. The proof works for all of \(D\), or for conductivity supported only on that quotient component. The covered slabs obstruct every original component, while the lower bounds use one specified component.

The component model includes the quotient component's full periodic lift. It does not select just one disconnected slab of that lift in \(\mathbb R^3\), which would generally fail to retain the original three-dimensional period lattice.

The exact contact-rod exclusion theorem states that rods with radius less than \(r/2\) add no probe-center exclusion over this probe-radius interval. Consequently \(D\), and the tensor defined by (1), remain unchanged when those prescribed capsules are added. Extra hinge housings are outside this statement.

The cell energy is normalized by \(V\). Dividing it by the conducting fraction \(\vartheta=|D|/V>0\) would preserve its rank but change its numerical value. Identifying such a normalized tensor with per-particle diffusivity requires additional dynamical and homogenization results. Those results are not supplied here.

**Counterexamples to broader claims.** A periodic grid of zero-width conducting curves can have three winding directions while occupying zero volume. Equation (1) then gives the zero tensor. A winding statement that permits such paths, or permits touching through excluded boundaries, does not supply the open route neighborhoods used in (2).

Bulk conductivity also need not have the maximum winding rank of any one connected component. In the unit cubic torus, take disjoint thin tubes around the loops \((t,0,0)\), \((1/3,t,1/3)\), and \((2/3,2/3,t)\), with radius less than \(1/20\). Each component has winding rank one. Their combined variational tensor has rank three by (3).

Uniform positive conductivity in the allowed region matters too. On the full unit torus, let scalar conductivity vanish quadratically near \(x_1\in\mathbb Z\). Geometric routes exist in every direction. A periodic potential whose entire residual field lies in a strip of width \(\delta\) has energy of order \(\delta\) in direction \(e_1\), so conductivity in that direction is zero.

These examples explain the hypotheses. The present result does not assume that geometric connectivity alone determines every transport law. Irregular-domain homogenization has separate analytical requirements, as discussed in [Heida, Jahnel and Vu, Sections 1 and 8](https://aimspress.com/aimspress-data/nhm/2025/1/PDF/nhm-20-01-010.pdf). No result from that stochastic setting is applied here.

**Replay and stopping point.** Run `python3 check.py --source-root /path/to/vitruvian3d` against the reviewed source. The script binds four files by SHA-256, replays the existing exact narrow passage verifier, and checks the new rational widths, length bounds, matrix inequalities and four rejected controls. It writes `certificate.json`. The written flux and variational arguments remain mathematical proof obligations for review.

This resolves the static rank question. It does not compute a transport time, prove a Brownian limit, establish a material design or derive a gravitational law.
