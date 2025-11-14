//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= The Einstein field equations
So far we've seen how we can write physical laws in general relativity by enforcing the principle of equivalence using the principle of general covariance. However, this requires we know what the metric describing some gravitational field looks like. So we want to know what the metric (the gravitational field) $g_(mu nu) (x)$ looks like given some distribution of matter.

== The energy-momentum tensor
In the Newtonian limit we found $ nabla^2 g_(00) (arrow(x)) = - 8 pi G_N rho(arrow(x)) $ with $rho(arrow(x))$ being the matter density which is equal to the energy density for small velocities. We'd like to generalize this for relativistic velocties, and from the LHS it seems like $rho(arrow(x))$ should be generalized to a tensor. In the case of a flow of matter then $rho(arrow(x))$ should be the energy density as seen by an observer moving with the flow given by the four-velocity $U^mu (x)$. This implies that the energy-momentum tensor of some relativistic flow
$
  T^(mu nu) = rho(x) U^mu (x) U^nu (x)
$
has the properties we are after.

Consider the $T^(00)$ component
$
  T^(00) = rho(x) (dv(x^0, tau))^2 = rho/(1-v^2)
$
with $arrow(v) = dd(arrow(x))\/dd(x^0)$ being the typical velocity, in special relativity
$
  dd(tau)^2 = (dd(x^0))^2 (1- arrow(v)^2)
$
in the Newtonian limit $T^(00) tilde.eq rho (x)$, so in this limit
$
  nabla^2 g_(00) (x) = - 8 G_N T_(00) (x)
$
The other components of $T^(mu nu)$ are
$
  T^(i j) = (rho v_i v_j)/(1-v^2) " and " T^(i 0) = T^(0 i) = (rho v_i)/(1 - v^2)
$
so $T^(mu nu)$ is the energy-momentum tensor in special relativity of a system with no pressure in the rest-frame, with $T^(0 i)$ being the density of momentum seen by a fixed observer, and $T^(i j)$ being the current of momentum.

In a closed system we have energy-momentum conservation
$
  partial_nu T^(mu nu) = 0
$
the $nu = 0$ component gives local energy conservation:
$
  partial_nu T^(0 nu) = pdv(, t) (rho/(1-v^2)) + nabla ((rho arrow(v))/(1-v^2)) = 0
$
this is sometimes called the relativistic Euler equation. The spatial part gives
$
  partial_nu T^(i nu) &= pdv(, t) ((rho v_i)/(1-v^2)) + pdv(, x) ((rho v_i v_x)/(1-v^2)) + pdv(, y) ((rho v_i v_y)/(1-v^2)) + pdv(, z) ((rho v_i v_z)/(1-v^2)) \
  &= rho/(1-v^2) (pdv(v^i, t) + arrow(v) dot nabla v^i)
$
this vanishes if the four-velocity is constant (no forces). This implies
$
  partial_nu T^(i nu) = rho/(1- v^2) (pdv(v^i, t) + arrow(v) dot nabla v^i) = 0
$
which is the relativistic Navier-Stokes equation. So conservation of the energy-momentum tensor in special relativity is a direct consequence of relativistic fluid mechanics.

Consider the total energy-momentum $P^mu$ defined by
$
  P^mu = integral dd(x, 3) T^(0 mu) (arrow(x),t)
$
this is conserved aswell
$
  dv(P^mu, t) = integral dd(x, 3) pdv(T^(0 mu), t) = integral dd(x^3) partial_nu T^(nu mu) = 0
$
where we use
$
  integral dd(x, 3) pdv(T^(i mu), x^i) = 0
$
if $T^(i mu)$ vanishes at infinity which is typical.

=== Perfect fluid
Previously we ignored pressure, now we consider a perfect fluid. A perfect fluid is defined by a velocity field $arrow(v) (x)$, with the property that any observer moving with the fluid sees it as being isotropic at all points---due to the cosmological principle our universe can be treated as a perfect fluid. To find $T^(mu nu)$ then we note that in the rest frame we have
$
  "is isotropic  " cases((T^(00))_"rest" = rho, (T^(i 0))_"rest"=(T^(0i))_"rest" = 0, (T^(i j))_"rest"=(T^(j i))_"rest" = p delta_(i j))
$
if the fluid has a velocity $arrow(v) (arrow(x),t)$ then this becomes
$
  T^(mu nu) = p eta^(mu nu) + (p + rho) U^mu U^nu
$
this reduces to the previous $T^(mu nu)$ in the rest frame, and is quite obviously a tensor. Now by conservation of this guy
$
  partial_nu T^(mu nu) = 0
$
yields the equations of motion of a perfect fluid in special relativity.

The covariant form is then easily achieved by $eta^(mu nu) -> g^(mu nu)$
$
  T^(mu nu) = p g^(mu nu) + (p + rho) U^mu U^nu
$
since in the free-falling frame $g^(mu nu) -> eta^(mu nu)$. Similarly to get the fluid equations in general relativity we write
$
  D_nu T^(mu nu) = 0
$
and now the total energy-momentum is
$
  P^mu = integral dd(x, 3) sqrt(-g) T^(mu 0)
$
this is no longer conserved! This is due to energy exchange between matter and the gravitational field---since we can write
$
  D_nu T^(mu nu) = 1/sqrt(-g) partial_nu (sqrt(-g) T^(mu nu)) + tensor(Gamma, mu, -rho nu) T^(rho nu)
$
the first term will vanish upon integrating, but the second term doesn't vanish. $P^mu$ is also not a contravariant vector.

== The curvature tensor
We have mentioned that a non-trivial metric does not imply a gravitational field, since assuming Minkowski spacetime we have
$
  g_(mu nu) (x) = eta_(mu nu) => dd(tau)^2 = - eta_(mu nu) dd(x^mu) dd(x^nu)
$
taking $x -> x'$ we obtain
$
        dd(tau^2) & = - g'_(mu nu) (x') dd(x'^mu) dd(x'^nu) \
  g'_(mu nu) (x') & = eta_(rho sigma) pdv(x^rho, x'^mu) pdv(x^sigma, x'^nu)
$
the metric seems to deviate from flat space even though it doesn't. We'd like to know how we can clearly distinguish between true curved space, and weirdly parametrized flat space. Suppose we had some covariant object which vanishes in flat space, but not in curved space. Then it would vanish in flat space independent of how we parametrize our space---therefore we could use this object to characterize our space.

We know $[partial_mu,partial_nu] = 0$ in flat space so if $[D_mu,D_nu] eq.not 0$ in curved space then we'd be done, since $D_mu -> partial_mu$ in flat space.

#proof[
  We show $[D_mu,D_nu] eq.not 0$. We recall the definition of the covariant derivative
  $
    D_kappa T_(mu nu) = partial_kappa T_(mu nu) - tensor(Gamma, lambda, -nu kappa) T_(mu lambda) - tensor(Gamma, lambda, -mu kappa) T_(lambda nu)
  $
  and let $T_(mu nu) equiv D_mu V_nu$. Then we obtain
  $
    D_(kappa) D_mu V_nu = partial_kappa D_mu V_nu - tensor(Gamma, lambda, -nu kappa) D_mu V_lambda - tensor(Gamma, lambda, -mu kappa) D_lambda V_nu
  $
  and commuting the derivatives
  $
    D_mu D_kappa V_nu = partial_mu D_kappa V_nu - tensor(Gamma, lambda, -nu mu) D_kappa V_lambda - tensor(Gamma, lambda, -kappa mu) D_lambda V_nu
  $
  then by symmetry the commutator is
  $
    D_kappa D_mu V_nu - D_mu D_kappa V_nu &= partial_kappa D_mu V_nu - partial_mu D_kappa V_nu - tensor(Gamma, lambda, -nu kappa) D_mu V_lambda + tensor(Gamma, lambda, -nu mu) D_kappa V_lambda \
    &= partial_kappa (partial_mu V_nu - tensor(Gamma, sigma, -mu nu) V_sigma) - partial_mu (partial_kappa V_nu - tensor(Gamma, sigma, -kappa nu) V_sigma) \ &#h(10pt)- tensor(Gamma, lambda, -nu kappa) (partial_mu V_lambda - tensor(Gamma, sigma, -mu lambda) V_sigma) + tensor(Gamma, lambda, -nu mu) (partial_kappa V_lambda - tensor(Gamma, sigma, -kappa lambda) V_sigma) \
    &= - partial_kappa tensor(Gamma, sigma, -mu nu) V_sigma - tensor(Gamma, sigma, -mu nu) partial_kappa V_sigma + partial_mu tensor(Gamma, sigma, -kappa nu)V_sigma + tensor(Gamma, sigma, -kappa nu) partial_mu V_sigma \
    &#h(10pt)- tensor(Gamma, lambda, -nu kappa) partial_mu V_lambda + tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda) V_sigma + tensor(Gamma, lambda, -nu mu) partial_kappa V_lambda - tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda) V_sigma \
    &= (-partial_kappa tensor(Gamma, sigma, -mu nu) + partial_mu tensor(Gamma, sigma, -kappa nu) - tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda)+ tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda) ) V_sigma \
    &= - tensor(R, sigma, -nu mu kappa) V_sigma
  $
  where we define the Riemann-Christoffel curvature tensor
  $
    tensor(R, sigma, -nu mu kappa) = partial_kappa tensor(Gamma, sigma, -mu nu) - partial_mu tensor(Gamma, sigma, -kappa nu) + tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda)
  $
  in Minkowski space $tensor(R, sigma, -nu mu kappa)$ vanishes since $D_mu -> partial_mu$ and $[partial_mu, partial_nu]=0$, and since it's a tensor it must then vanish independent of coordinates. So in flat spacetime it always vanishes and we are done.

]

We can also construct a purely covariant tensor
$
  R_(lambda mu nu kappa) = g_(lambda sigma) tensor(R, sigma, -nu mu kappa)
$
which becomes
$
  R_(lambda mu nu kappa) = 1/2 (pdv(g_(lambda nu), x^kappa, x^mu)-pdv(g_(mu nu), x^kappa, x^lambda) - pdv(g_(lambda kappa), x^nu, x^mu)+pdv(g_(mu kappa), x^nu, x^lambda)) + g_(eta sigma) (tensor(Gamma, eta, -nu lambda) tensor(Gamma, sigma, -mu kappa) - tensor(Gamma, eta, -kappa lambda) tensor(Gamma, sigma, -mu nu))
$
#proof[ (sketch)][
  The second part is trivial so we just show the first part
  $
    g_(lambda sigma) (partial_kappa tensor(Gamma, sigma, -mu nu) - partial_mu tensor(Gamma, sigma, -kappa nu))
  $
  the two terms are very similar, consider
  $
    g_(lambda sigma) partial_kappa tensor(Gamma, sigma, -mu nu) &= 1/2 partial_kappa [ g^(sigma alpha) (pdv(g_(mu alpha), x^nu)+pdv(g_(nu alpha), x^mu) - pdv(g_(mu nu), x^alpha))] \
    &= 1/2 g_(lambda sigma) {partial_k g^(sigma alpha) [partial_nu g_(mu alpha) + partial_mu g_(nu alpha) - partial_alpha g_(mu nu)] \ &#h(45pt) + g^(sigma alpha) [partial_kappa partial_nu g_(mu alpha)+partial_kappa partial_mu g_(nu alpha) - partial_kappa partial_alpha g_(mu nu)]}
  $
  when taking the difference with the same but $kappa <--> mu$ then we get two kinds of terms, second derivatives and the rest. The second derivative part becomes
  $
    partial_kappa partial_nu g_(mu lambda) - partial_kappa partial_lambda g_(mu nu) - partial_mu partial_nu g_(kappa lambda) + partial_mu partial_lambda g_(kappa nu)
  $
  which is the first part.
]

This representation nicely shows how it only vanishes if $partial^2 g = 0$, and easily shows its symmetries (and the first Bianchi identity)
$
  R_(lambda mu nu kappa) &= R_(nu kappa lambda mu) "     " (lambda mu)(nu kappa) <-> (nu kappa) (lambda mu) \
  R_(lambda mu nu kappa) &= - R_(mu lambda nu kappa) "   " lambda <-> mu \
  R_(mu lambda nu kappa) &= R_(lambda mu kappa nu) "     " lambda <-> mu, kappa <-> nu \
  R_(lambda mu nu kappa) &= - R_(lambda mu kappa nu) "   " kappa <-> nu \
  0 &= R_(lambda mu nu kappa) + R_(lambda kappa mu nu) + R_(lambda nu kappa mu) "cyclic permutation of "(mu, nu, kappa)
$
these are easily seen to hold in the local elevator where the Christoffel symbols vanish, and therefore they hold everywhere.

We define the Ricci tensor
$
  R_(mu kappa) = g^(lambda nu) R_(lambda mu nu kappa) = tensor(R, lambda, -mu lambda kappa)
$
which is symmetric by the previous. Further we can define the Ricci scalar
$
  R = g^(mu kappa) R_(mu kappa) = tensor(R, mu, -mu)
$
Lastly we'll need the second Bianchi identity
$
  D_eta R_(lambda mu nu kappa) + D_kappa R_(lambda mu eta nu) + D_nu R_(lambda mu kappa eta) = 0
$
which is most easily verified by going to the local elevator where $D -> partial$ and the Christoffel symbols vanish. Now in a local elevator we know $partial_sigma g_(mu nu) = 0$ so
$
  D_sigma g_(mu nu) = 0
$
everywhere. This lets us rewrite the second Bianchi identity by contracting with $g^(lambda nu)$
$
  D_eta R_(mu kappa) - D_kappa R_(mu eta) + D_nu tensor(R, nu, -mu kappa eta)
$
and now contracting by $g^(mu kappa)$
$
  D_eta R - D_mu tensor(R, mu, -eta)-D_nu tensor(R, nu, -eta) &= 0 \
  D_mu (tensor(R, mu, -eta) - 1/2 tensor(delta, mu, -eta) R) &= 0 \
  D_mu underbrace((R^(mu nu) - 1/2 g^(mu nu) R), "trace-reversed Ricci tensor") &= 0
$
so the trace-reversed Ricci tensor is conserved.

== The field equations
Recall we have
$
  nabla^2 g_(00) (x) = - 8 pi G_N T_(0 0) (x)
$
naturally this generalizes as
$
  G_(mu nu) (x) = - 8 pi G_N T_(mu nu) (x)
$
where $G_(mu nu)$ depends on the metric and its derivatives. By dimensional analysis the LHS should have units $[L^(-2)]$ which matches the second derivative of $g_(mu nu)$. For higher derivative we'd need to multiply by some fundamental length $L_p$. We might guess
$
  G_(mu nu) eq^? R_(mu nu)
$
but we also require
$
  D_nu G^(mu nu) =^! 0 " since " D_nu T^(mu nu) = 0
$
this implies that
$
  G_(mu nu) = c (R_(mu nu) - 1/2 g_(mu nu) R)
$
where $c$ is some constant to be determined.

In the Newtonian limit $abs(T_(i j)) << abs(T_(00))$ so the zeroth component becomes
$
  G_(00) = c (R_(00) + 1/2 R) = - 8 pi G_N T_(00)
$
with $T_(i j) -> 0$ implying
$
  R_(i j) - 1/2 g_(i j) R -> 0 " so " R_(i j) tilde.eq 1/2 g_(i j) R
$
in the weak field limit $g_(i j) eq.not n_(i j)$ so
$
  R tilde.eq R_(i i) - R_(00) tilde.eq 3/2 R - R_(0 0) => 1/2 R tilde.eq R_(00)
$
so we obtain
$
  G_(00) = c(R_(00)+1/2 R) = 2 c R_(00)
$
now in the Newtonian limit
$
  R_(lambda mu nu kappa) tilde.eq 1/2 (pdv(g_(lambda nu), x^kappa, x^mu) - pdv(g_(mu nu), x^nu, x^lambda) - pdv(g_(lambda kappa), x^nu, x^mu) + pdv(g_(mu kappa), x^nu, x^lambda))
$
all time-derivatives vanish in the static limit so
$
  R_(0000) tilde.eq 0 " and " R_(i 0 j 0) tilde.eq 1/2 pdv(g_(00), x^i, x^j)
$
and we obtain
$
  G_(00) tilde.eq 2 c R_(0 0) tilde.eq 2 c(R_(i 0 i 0) - R_(0 0 0 0)) tilde.eq c nabla^2 g_(0 0)
$
implying $c = 1$. Finally we obtain the Einstein field equations
$
  G_(mu nu) = R_(mu nu) - 1/2 g_(mu nu) R = - 8 pi G_N T_(mu nu)
$
which can be written in trace-inverted form
$
  R_(mu nu) = - 8 pi G_N (T_(mu nu) - 1/2 g_(mu nu) tensor(T, alpha, -alpha))
$
which is easily derived by taking the trace of the Einstein equation, from this we see $T_(mu nu) = 0 => R_(mu nu) = 0$.

We can generalize the Einstein equation by including a cosmological constant
$
  R_(mu nu) - 1/2 g_(mu nu) R - Lambda g_(mu nu) = - 8 pi G_N T_(mu nu)
$
this is allowed since the covariant derivative of $g_(mu nu)$ vanishes. Moving $Lambda g_(mu nu)$ to the RHS it can be interpreted as an energy density present even if $T_(mu nu)$ vanishes---vacuum energy or dark energy.

#pagebreak()
= Einstein-Hilbert action
We briefly discuss an alternate derivation of the EFE since it is nice.

== Without matter
The simplest action we can create using only the metric is
$
  S = integral dd(x, 4) sqrt(-g) R
$
writing $R = g^(mu nu) R_(mu nu)$ we vary the metric to find
$
  dd(S, d: delta) = integral dd(x, 4) [(dd(sqrt(-g), d: delta))g^(mu nu) R_(mu nu) + sqrt(-g) (dd(g^(mu nu), d: delta))R_(mu nu) + sqrt(-g) g^(mu nu) dd(R_(mu nu), d: delta)]
$
note that
$
  g_(rho mu) g^(mu nu) = tensor(delta, -rho, nu) => dd(g^(mu nu), d: delta) = - g^(mu rho) g^(nu sigma) dd(g_(rho sigma), d: delta)
$
we claim
$
  dd(sqrt(-g), d: delta) = - 1/2 sqrt(-g) g_(mu nu) dd(g^(mu nu), d: delta)
$
this follows from the Jacobi formula $log det A = tr log A$. The second term already has a nice form now we show
$
  dd(R_(mu nu), d: delta) = nabla_rho dd(tensor(Gamma, rho, -mu nu), d: delta) - nabla_nu dd(tensor(Gamma, rho, -mu rho), d: delta)
$
(it is a total derivative) where
$
  dd(tensor(Gamma, rho, -mu nu), d: delta) = 1/2 g^(rho sigma) (nabla_mu dd(g_(sigma nu), d: delta) + nabla_nu dd(g_(sigma mu), d: delta)- nabla_sigma dd(g_(mu nu), d: delta))
$
note $dd(tensor(Gamma, rho, -mu nu), d: delta)$ is a tensor due to it being a difference between two Christoffel symbols. This enables us to work in a freely-falling elevator where $tensor(Gamma, rho, -mu nu) = 0$, giving to linear order
$
  dd(tensor(Gamma, rho, -mu nu), d: delta) &= 1/2 g^(rho sigma) (partial_mu dd(g_(sigma nu), d: delta) + partial_nu dd(g_(sigma mu), d: delta)-partial_sigma dd(g_(mu nu), d: delta)) \
  &= 1/2 g^(rho sigma) (nabla_mu dd(g_(sigma nu), d: delta) + nabla_nu dd(g_(sigma mu), d: delta)-nabla_sigma dd(g_(mu nu), d: delta))
$
and since the LHS is a tensor the RHS is also a tensor meaning this is valid everywhere not only in the freely-falling elevator.

Next consider in a freely-falling elevator:
$
  tensor(R, sigma, -rho mu nu) = partial_mu tensor(Gamma, sigma, -nu rho) - partial_nu tensor(Gamma, sigma, -mu rho)
$
giving
$
  dd(tensor(R, sigma, -rho mu nu), d: delta) &= partial_mu dd(tensor(Gamma, sigma, -nu rho), d: delta) - partial_nu dd(tensor(Gamma, sigma, -mu rho), d: delta) \
  &=^"everywhere" nabla_mu dd(tensor(Gamma, sigma, -nu rho), d: delta) - nabla_nu dd(tensor(Gamma, sigma, -mu rho), d: delta)
$
then to leading order
$
  dd(R_(rho nu), d: delta) = nabla_mu dd(tensor(Gamma, mu, -nu rho), d: delta)-nabla_nu dd(tensor(Gamma, mu, -rho mu), d: delta)
$
this is nice because
$
  g^(mu nu) dd(R_(mu nu), d: delta) = nabla_mu X^mu "with" X^mu = g^(rho nu) dd(tensor(Gamma, mu, -rho nu), d: delta) - g^(mu nu) dd(tensor(Gamma, rho, -nu rho), d: delta)
$
then combining results
$
  dd(S, d: delta) = integral dd(x, 4) sqrt(-g) [(R_(mu nu)-1/2 R g_(mu nu))dd(g^(mu nu), d: delta) + underbrace(nabla_mu X^mu, "total derivative")]
$
requiring $dd(S, d: delta) = 0$ then gives
$
  G_(mu nu) equiv R_(mu nu) - 1/2 R g_(mu nu) = 0
$
contracting with $g^(mu nu)$ gives $R = 0$ so the vacuum EFE are
$
  R_(mu nu) = 0
$
we can get a slightly more non-trivial action by multiplying the volume form by a constant
$
  S = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (R - 2 Lambda)
$
varying the action then gives an additional term
$
  R_(mu nu) - 1/2 R g_(mu nu) = - Lambda g_(mu nu)
$
now contracting with $g^(mu nu)$ gives $R = 4 Lambda$ meaning
$
  R_(mu nu) = Lambda g_(mu nu)
$
these are the vacuum EFE in the presence of a cosmological constant.

== With matter
Recall in Minkowski space a scalar field has the action
$
  S_"scalar" = integral dd(x, 4) (-1/2 eta^(mu nu) partial_mu phi.alt partial_nu phi.alt - V(phi.alt))
$
this trivially generalizes to curved space by
$
  S_"scalar" = integral dd(x, 4) sqrt(-g) (-1/2 g^(mu nu) nabla_mu phi.alt nabla_nu phi.alt - V(phi.alt))
$
of course here $partial_mu -> nabla_mu$ is redundant since $phi.alt$ is a scalar field. However, curved space enables us to add new terms e.g.
$
  S_"scalar" = integral dd(x, 4) sqrt(-g) (-1/2 g^(mu nu) nabla_mu phi.alt nabla_nu phi.alt - V(phi.alt) - 1/2 xi R phi.alt^2)
$
for a constant $xi$, and the new term vanishes for $g_(mu nu) = eta_(mu nu)$. We find
$
  dd(S_"scalar", d: delta) &= integral dd(x, 4) sqrt(-g) (-g^(mu nu) nabla_mu dd(phi.alt, d: delta) nabla_nu phi.alt - pdv(V, phi.alt) dd(phi.alt, d: delta)- xi R phi.alt dd(phi.alt, d: delta)) \
  &= integral dd(x, 4) sqrt(-g) [underbrace((g^(mu nu) nabla_mu nabla_nu phi.alt - pdv(V, phi.alt) - xi R phi.alt), = 0)dd(phi.alt, d: delta)-nabla_mu (dd(phi.alt, d: delta) nabla^mu phi.alt)]
$
the Maxwell action also generalizes
$
  S_"maxwell" = -1/4 integral dd(x, 4) sqrt(-g) F_(mu nu) F^(mu nu)
$
with
$
  F_(mu nu) = nabla_mu A_nu - nabla_nu A_mu
$
giving
$
  nabla_mu F^(mu nu) = 0
$
as expected.

The Einstein-Hilbert action becomes
$
  S = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (R-2 Lambda) + underbrace(S_M, "matter")
$
we define
$
  T_(mu nu) = - 2/sqrt(-g) dv(S_M, g^(mu nu), d: delta)
$
then
$
  dd(S, d: delta) = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (G_(mu nu)+Lambda g_(mu nu)) dd(g^(mu nu), d: delta) - 1/2 integral dd(x, 4) sqrt(-g) T_(mu nu) dd(g^(mu nu), d: delta)
$
immediately giving the full EFE
$
  G_(mu nu) + Lambda g_(mu nu) = 8 pi G T_(mu nu)
$
or just
$
  G_(mu nu) = 8 pi G T_(mu nu)
$


#pagebreak()
= A detour
== Geodesic as minimal curve
We'd like a geometric interpretation of the geodesic, as we shall see it is the shortest path between two points in spacetime. Consider the line element $dd(s^2) = g_(mu nu) dd(x^mu) dd(x^nu)$, then the total line length is
$
  S = integral dd(s) = integral L dd(tau)
$
where
$
  L = sqrt(g_(mu nu) dv(x^mu, tau) dv(x^nu, tau))
$
so $L$ is a function of $x$ and $dd(x)\/dd(tau) = dot(x)$, as we'd expect for a classical Lagrangian. We let a curve be parametrized by $x^mu (tau)$ then the shortest curve between $x_i^mu$ and $x_f^mu$ is obtained by minimizing the action $dd(S, d: delta) = 0$. This easily leads one to the Euler-Lagrange equations
$
  pdv(L, x^mu) - dv(, tau) pdv(L, dot(x)^mu) = 0
$
as should be obvious the resulting equations of motion from using $tilde(L) = L^2$ are the same as those obtained by using $L$, so we use $tilde(L)=g_(mu nu) dot(x)^mu dot(x)^nu$ and
$
  pdv(L, dot(x)^mu) & = 2 g_(mu nu) dot(x)^nu \
       pdv(L, x^mu) & = pdv(g_(nu lambda), x^mu) dot(x)^nu dot(x)^lambda
$
to obtain
$
  0 &= pdv(g_(mu nu), x^lambda) dv(x^lambda, tau) dv(x^nu, tau) + g_(mu nu) dv(x^nu, tau, 2) - 1/2 pdv(g_(nu lambda), x^mu) dv(x^nu, tau) dv(x^lambda, tau) \
  0 &=^(g^(sigma mu) (dots)) dv(x^sigma, tau, 2) + 1/2 g^(sigma mu) (2 pdv(g_(mu nu), x^lambda) - pdv(g_(nu lambda), x^mu)) dv(x^lambda, tau) dv(x^nu, tau) \
  0 &= dv(x^sigma, tau, 2) + tensor(Gamma, sigma, -lambda nu) dv(x^lambda, tau) dv(x^nu, tau)
$
which is the familiar geodesic equation.

=== A trick
Comparing the equation obtained by
$
  dd(integral g_(mu nu) dv(x^mu, tau) dv(x^nu, tau) dd(tau), d: delta) = 0
$
with the geodesic equation it should be clear that by comparison we can obtain the Christoffel symbols---this is typically much simpler.

== The time-dependent spherically symmetric metric
We'll apply the previous trick to work toward some non-trivial solutions of the field equations---for simplicity we assume spherical symmetry, anything else goes.

For the metric to be spherically symmetric $dd(tau^2)$ can only depend on rotationally invariant quantities, with $r^2 = x^2 + y^2 + z^2$ these are $ {t, dd(t), r, r dd(r) = bold(x) dd(bold(x)), dd(r^2) + r^2 (dd(theta^2) + sin^2 theta dd(phi^2)) = dd(bold(x)^2)} $ defining the solid angle $dd(Omega^2) = dd(theta^2) + sin^2 theta dd(phi^2)$ then the metric must take the form
$
  dd(tau^2) = A dd(t^2) - B dd(r^2) - C dd(r, t) - D r^2 dd(Omega^2)
$
with ${A,B,C,D}$ being arbitrary functions of $t$ and $r$. We'll now massage and rewrite this until we get something more workable. Since relativity is coordinate invariant by construction we start by absorbing $D$ into $r$ by redefining $r -> r' = r sqrt(D)$, then dropping primes for convenience we have
$
  dd(tau^2) = A dd(t^2) - B dd(r^2) - C dd(r, t) - r^2 dd(Omega^2)
$
where the $A, B$ and $C$ are changed, similarly by time-reparametrization we can absorb $C$ by
$
  dd(t') = eta(r, t) [A dd(t)-1/2 C dd(r)]
$
note that $eta$ is not independent since
$
  eta A = pdv(t', t)";  " 1/2 eta C = -pdv(t', r)
$
which can be solved knowing $A$ and $C$ to give $eta$. Then
$
  1/(A eta^2) dd(t'^2) = A dd(t^2)- C dd(t, r) + C^2/(4 A) dd(r^2)
$
implying
$
  dd(tau^2) = 1/(eta^2 A) dd(t^2) - (B+ C^2/(4 A)) dd(r^2) - r^2 dd(Omega^2)
$
we write this as
$
  dd(tau^2) = E(r, t) dd(t^2) - F(r,t) dd(r^2) - r^2 dd(Omega^2)
$
so we have simplified our problem immensily. As a side-effect the metric has become diagonal, meaning it is very simple to determine the components
$
  g_(r r) &= F";  " g_(theta theta) = r^2";  " g_(phi phi) = r^2 sin^2 theta";  " g_(t t) = - E \
  g^(r r) &= 1/F";  " g^(theta theta) = 1/r^2";  " g^(phi phi) = 1/(r^2 sin^2 theta)";  " g^(t t) = -1/E
$
and
$
  -g = - det g_(mu nu) = r^4 E F sin^2 theta
$
meaning
$
  sqrt(-g) dd(x, 4) = r^2 sin theta sqrt(E F) dd(r, theta, phi, t)
$
if $E = F = 1$ this reduces to the usual spherical Jacobian.

Now we compute the Christoffel symbols by $dd(integral dd(tau), d: delta) = 0$ in our case this becomes
$
  dd(integral dd(tau) [E dot(t)^2 - F dot(r)^2 - r^2 dot(theta)^2 - r^2 sin^2 theta dot(phi)^2], d: delta) = 0
$
with $dot(x) = dd(x)\/dd(tau)$. Consider the $mu = 0$ component of the Euler-Lagrange equation
$
  pdv(L, x^mu) = dv(, tau) pdv(L, dot(x)^mu)
$
giving
$
  dot(t)^2 pdv(E, t) - pdv(F, t) dot(r)^2 &= dv(, tau) (2 E dot(t)) \
  &= 2 E dot.double(t) + 2 dot(t) dv(E, tau) \
  &= 2 E dot.double(t) + 2 dot(t) (pdv(t, tau) pdv(E, t) + pdv(r, tau) pdv(E, r)) \
  &= 2 E dot.double(t) + 2 dot(t)^2 pdv(E, t) + 2 dot(t) dot(r) pdv(E, r)
$
simplifying
$
  0 &= pdv(F, t) dot(r)^2 + 2 E dot.double(t) + dot(t)^2 pdv(E, t) + 2 dot(t) dot(r) pdv(E, r) \
  0 &= dot.double(t) + 1/(2 E) pdv(E, t) dot(t)^2 + 1/E pdv(F, r) dot(t) dot(r) + 1/(2 E) pdv(F, t) dot(r)^2
$
we compare this with
$
  0 = dot.double(t) + tensor(Gamma, 0, -mu nu) dot(x)^mu dot(x)^nu
$
to obtain the Christoffel symbols
$
  tensor(Gamma, t, -t t) & = 1/(2 E) pdv(E, t) \
  tensor(Gamma, t, -r r) & = 1/(2 E) pdv(F, t) \
  tensor(Gamma, t, -r t) & = tensor(Gamma, t, -t r) = 1/(2 E) pdv(E, r)
$
similarly the $r$-component gives
$
          tensor(Gamma, r, -t r) & = tensor(Gamma, r, -r t) = 1/(2 F) pdv(F, t) \
          tensor(Gamma, r, -r r) & = 1/(2 F) pdv(F, r) \
          tensor(Gamma, r, -t t) & = 1/(2 F) pdv(E, r) \
  tensor(Gamma, r, -theta theta) & = - r/F \
      tensor(Gamma, r, -phi phi) & = - (r sin^2 theta)/F
$
and the $theta$-component gives
$
  tensor(Gamma, theta, -r theta) & = tensor(Gamma, theta, theta r) = 1/r \
  tensor(Gamma, theta, -phi phi) & = - sin theta cos theta
$
and the $phi$-component gives
$
  tensor(Gamma, phi, -r phi) &= tensor(Gamma, phi, -phi r) = 1/r \
  tensor(Gamma, phi, -theta phi) &= tensor(Gamma, phi, -phi theta) = (cos theta)/(sin theta)
$
=== Ricci tensor
We have by definition
$
  tensor(Gamma, mu, -mu nu) = 1/2 g^(mu sigma) pdv(g_(mu sigma), x^nu)
$
then using
$
  g^(mu sigma) pdv(g_(mu sigma), x^nu) = pdv(ln g, x^nu)
$
and $R_(mu kappa) = tensor(R, lambda, -mu lambda kappa)$ with
$
  tensor(R, sigma, -mu nu kappa) = partial_kappa tensor(Gamma, sigma, -mu nu) - partial_nu tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu nu) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -mu kappa) tensor(Gamma, sigma, -nu lambda)
$
we find
$
  R_(mu kappa) &= partial_kappa tensor(Gamma, sigma, -sigma mu) - partial_sigma tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu sigma) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -mu kappa) tensor(Gamma, sigma, -sigma lambda) \
  &= 1/2 partial_kappa partial_mu ln g - partial_sigma tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu sigma) tensor(Gamma, sigma, -kappa lambda) - 1/2 tensor(Gamma, lambda, -mu kappa) partial_lambda ln g
$
then the components are
$
  R_(r r) &= - 1/(2 E) pdv(E, r, 2) + 1/(4 E^2) (pdv(E, r))^2 + 1/(4 E F) pdv(E, r) pdv(F, r) \
  &#h(12pt)+ 1/(r F) pdv(F, r) + 1/(2 E) pdv(F, t, 2) - 1/(4 E^2) pdv(E, t) pdv(F, t) - 1/(4 E F) (pdv(F, t))^2 \
  R_(theta theta) &= -1 + 1/F - r/(2 F^2) pdv(F, r) + r/(2 E F) pdv(E, r) \
  R_(t t) &= 1/(2 F) pdv(E, r, 2) - 1/(4 F^2) pdv(E, r) pdv(F, r) + 1/(r F) pdv(E, r) \
  &#h(12pt)- 1/(4 E F) (pdv(E, r))^2 - 1/(2 F) pdv(F, t, 2) + 1/(4 F^2) (pdv(F, t))^2 + 1/(4 E F) pdv(E, t) pdv(F, t) \
  R_(t r) &= 1/(r F) pdv(F, t) \
  R_(phi phi) &= sin^2 theta R_(theta theta)
$
then the trace-inverted form of the field equations give a differential equation for $E$ and $F$ in terms of $T_(mu nu)$ by
$
  R_(mu nu) = - 8 pi G_N (T_(mu nu) - 1/2 g_(mu nu) tensor(T, lambda, -lambda))
$
in our case $g_(t r) = 0$ so e.g. $ 1/(r F) pdv(F, t) = 8 pi G T_(t r) $

#pagebreak()
= The Schwarzchild metric
We now try to get one of the most important solutions to the Einstein equations, namely the Schwarzchild solution.

== The metric
We consider a point mass with mass $M$ at the origin, and assume that $m_"obs" << M$ so we can ignore any backreaction. This is quite obviously spherically symmetric, meaning we can use the previous general metric,
$ dd(tau^2)=E dd(t^2) - F dd(r^2) - r^2 dd(Omega^2) $
and importantly $T_(mu nu) = 0$ for all $r eq.not 0$, so by the Einstein equations $R_(mu nu) = 0$. In particular $R_(r r) = R_(t t) = 0$ so
$
  0 &= R_(r r)/F + R_(t t)/E \
  &= - 1/(2 E F) pdv(E, r, 2) + 1/(4 E^2 F) (pdv(E, r))^2 + 1/(4 E F^2) pdv(E, r) pdv(F, r) \
  &#h(12pt)+ 1/(r F^2) pdv(F, r) + 1/(2 E F) pdv(F, t, 2) - 1/(4 E^2 F) pdv(E, t) pdv(F, t) - 1/(4 E F^2) (pdv(F, t))^2 \
  &#h(12pt)+1/(2 E F) pdv(E, r, 2) - 1/(4 E F^2) pdv(E, r) pdv(F, r) + 1/(r E F) pdv(E, r) \
  &#h(12pt)- 1/(4 E^2 F) (pdv(E, r))^2 - 1/(2 E F) pdv(F, t, 2) + 1/(4 E F^2) (pdv(F, t))^2 + 1/(4 E^2 F) pdv(E, t) pdv(F, t) \
  0 &= 1/(r F^2) pdv(F, r) + 1/(r E F) pdv(E, r) \
  0 &= 1/(r F) underbrace((1/F pdv(F, r) + 1/E pdv(E, r)), =^!0)
$
assuming $E$ and $F$ are time-independent we can then write
$
  pdv(ln F, r) = - pdv(ln E, r)
$
implying $E(r) F(r) = "const"$, and we've gone from four independent functions down to one. We can fix the constant with our boundary conditions. In the limit $r->oo$ the metric should reduce to $eta_(mu nu)$ where $E = F = 1$ so $"const" = 1$ meaning
$
  E = 1/F
$
for all $r$. Now consider
$
  0 & = R_(theta theta) \
  0 & = -1 + E + r pdv(E, r) \
    & => pdv(, r) (r E) = 1 \
$
which can be solved to give $r E = r + C$ or
$
  E = 1 + C/r
$
we can fix this constant by considering the Newtonian limit where
$
  g_(0 0) = - (1 + 2 Phi) =^(Phi = - G_N M r^(-1)) - 1 + (2 G_N M)/r
$
since $E = - g_(0 0)$ we immediately find $C = -2 G_N M$ so
$
  E = 1 - (2 M G_N)/r;"   "F = 1/(1-(2 M G_N)/r)
$
and we've found the Schwarzchild metric
$
  dd(tau^2) = (1 - (2 M G_N)/r) dd(t^2) - dd(r^2)/(1-(2 M G_N)/r) - r^2 (dd(theta^2) + sin^2 theta dd(phi.alt^2))
$
which is valid outside any spherically symmetric point mass.

== Birkhoff's theorem
Birkhoff's theorem is very important and useful since it allows us to extend the Schwarzchild metric to any spherically symmetric mass distribution, e.g. stars, planets, etc. In fact we could even have a mass distribution which expands or contracts in time, as long as the total mass is unchanged.

#thm[Birkhoff's theorem][
  A spherically symmetric gravitational field in empty space must be static with a metric given by the Schwarzchild metric.
]
#proof[
  Our main assumption was that $E$ and $F$ are time-independent, this can be shown to always be the case. From $R_(t r) = 0$ we have
  $
    1/(r F) pdv(F, t) = 0 => pdv(F, t) = 0
  $
  so this is nice. Now in principle we can write
  $
    E = f(t) (1-(2 M G_N)/r) eq.not E(r)
  $
  but we can just redefine our time-coordinate (reset our time) to absorb $f(t)$ by
  $
    t' = integral sqrt(f(t)) dd(t) => dd(t') = sqrt(f(t)) dd(t)
  $
  so $E$ and $F$ are always time-independent, meaning everything done before is valid as long as $T_(mu nu) = 0$. So we can have a weird oscillatory mass distribution, we just require it be spherically symmetric, then the Schwarzchild metric is valid in vacuum.

]

This also explains why it's relatively rare to achieve gravitational waves since we need to break spherical symmetry for the metric to be changing in time. This also makes it trivial to analyze cavities since these are spherically symmetric and $T_(mu nu) = 0$, but $M = 0$ so in cavities we have $eta_(mu nu)$!

== Precession of Mercury
The precession of Mercury's perihelion is famously one of the first testable predictions made by relativity---since this is not explained by Newtonian gravity.

We want to solve the geodesic equation, so we consider the Lagrangian (setting $theta = pi\/2$ for simplicity)
$
  L = - (1-R/r) dot(t)^2 + 1/(1-R/r) dot(r)^2 + r^2 dot(phi)^2
$
where we've defined the Schwarzchild radius $R = 2 M G_N$---we'll also use that $L = dd(s^2)\/dd(tau^2) = -1$. To solve this we introduce Killing vectors.

=== Killing vectors & constants of motion
In the Schwarzchild metric all $g_(mu nu)$ are time-independent. In general if $g_(mu nu)$ is independent of $x^0$ then
$
  pdv(L, x^0) = 0
$
and by the Euler-Lagrange equation we have a conserved quantity
$
  dv(, tau) underbrace((pdv(L, dot(x)^0)), "constant") =^! 0
$
so
$
  pdv(L, dot(x)^0) = pdv(L, dot(t)) = 2 g_(mu 0) dot(x)^mu
$
is a conserved quantity, but this guy is not covariant! To ensure covariance we define the Killing vector
$
  K^mu_((t)) = vecrow(1, 0, 0, 0)
$
then
$
  gamma_((t)) = g_(mu nu) dot(x)^mu K^nu_((t))
$
is covariant and conserved. In general
$
  gamma_((alpha)) = g_(mu nu) dot(x)^mu K^nu_((alpha))
$
is covariant and conserved if $g_(mu nu)$ is independent of $x^alpha$.

The Schwarzchild metric only depend on $r$, so we have a conserved quantity related to both $t$ and $phi$ given by
$
  kappa &equiv - g_(mu nu) dot(x)^mu K^nu_((t)) = - g_(00) dot(t) = (1- R/r) dot(t) \
  l &equiv m g_(mu nu) dot(x)^mu K^nu_((phi)) = m g_(phi phi) dot(phi) = m r^2 phi
$
we recognize $kappa$ as an energy and $l$ as the angular momentum---the inclusion of $m$ and $-1$ are just for convenience since multiplying by a scalar still gives a conserved quantity. Substituting these into $L$ and letting $L=-1$ we find
$
  - kappa^2/(1-R/r) + dot(r)^2/(1-R/r) + l^2/(m^2 r^2) = - 1
$
which can in principle be solved. To see what it means we multiply by $1/2 m (1 - R/r)$ and define
$
  E/m equiv (kappa^2 -1)/2
$
giving
$
  underbrace(1/2 m r^2, E_"kin") + underbrace((1- R/r), "GR corr.") underbrace(l^2/(2 m r^2), E_l) - underbrace((G_N m M)/r, E_g) = E
$
which expresses the usual conservation of energy with a relativistic correction. The correction makes the orbit open and solving gives
$
  r = alpha/(1 + e cos[(1-epsilon.alt)phi])
$
with $e$ being the eccentricity and
$
        alpha & = l^2/(G M m^2) = (1 + e) r_"min" \
  epsilon.alt & = (3 R)/(2 alpha)
$
orbit returns to $r_"min"$ not at $phi = 2 pi$ as we'd expect, but at
$
  phi = (2 pi)/(1 - epsilon.alt) tilde.eq 2 pi + (3 pi R)/alpha
$
so the orbit precesses.

#pagebreak()
= Black holes
We write the Schwarzchild metric as
$
  dd(tau^2) = f(r) dd(t^2) - dd(r^2)/f(r) - r^2 dd(Omega^2)
$
where we've defined
$
  f(r) equiv 1 - R/r
$
with the Schwarzchild radius $R = 2 M G$. We notice that at $r = R$ we get a singularity, in many cases this is not a problem since $R < R_"mass"$. What about a point mass or some object with $R_"mass" < R$---what we'd call a black hole?

== The geodesic equation
We want to know what happens when crossing $R$. We consider a free-falling observer, so we use the geodesic equation
$
  dv(x^alpha, tau, 2) + tensor(Gamma, alpha, -mu nu) dv(x^mu, tau) dv(x^nu, tau) = 0
$
or we use the principle of least action to find
$
  0 & = dd(integral [g_(mu nu) dot(x)^mu dot(x)^nu] dd(tau), d: delta) \
    & = dd(integral [-f(r) dot(t)^2 + f^(-1) dot(r)^2], d: delta)
$
where we assume the observer is falling radially $dot(theta)=dot(phi)=0$, so we have
$
  L = f(r) dot(t)^2 - f^(-1) (r) dot(r)^2
$
this only depends on $r$, so again by using Killing vectors we have the conserved quantity
$
  C = g_(mu nu) dot(x)^mu K^nu_((t)) = f(r) dot(t)
$
we can fix $C$ by considering the limit $r -> oo$ in this case $dot(t) -> 1$, $f(r) -> 1$ so $C = 1$ meaning $dot(t) = f^(-1) (r)$. So our Lagrangian becomes
$
  L = - 1/f(r) + f^(-1) (r) dot(r)^2 = - 1
$
since $ 1 = dv(tau^2, tau^2) = 1/dd(tau^2) [- g_(mu nu) dd(x^mu) dd(x^nu) dv(tau^2, tau^2)] = - g_(mu nu) dot(x)^mu dot(x)^nu = -L $
where we use $dd(tau^2) = -g_(mu nu) dd(x^mu) dd(x^nu)$. So
$
  1/f(r) - f^(-1) (r) dot(r)^2 =1
$
giving
$
  dot(r)^2 & = 1 - f(r) \
           & = 1 - (1 - R/r) \
           & = R/r = (2 M G)/r
$
which is nice and simple.

=== Consequences
Now we'd like to know what the free-falling observer experiences by looking at $tau$, and what an observer far from the black hole experience by looking at $t$. Taking the root we find
$
  sqrt(r) dot(r) = - sqrt(2 M G) => 2/(3 sqrt(2 M G)) (r^(3\/2)-r_0^(3\/2)) = tau_0 - tau
$
here we see there is nothing special at $r = R$. We even find for $r = 0$
$
  tau = tau_0 + 2/(3 sqrt(2 M G)) r_0^(3\/2) < oo
$
so the free-falling observer reaches the center of the black hole in finite time. For $t$ we do
$
  dv(r, t) = dot(r)/dot(t) = - sqrt((2 M G)/r) (1 - (2 M G)/r)
$
which has an exact solution, which is rough, so we just look at limits. In the limit $r >> 2 M G$ this becomes
$
  dv(r, t) =^(r << R) - sqrt((2 M G)/r)
$
which is the same as before so $t tilde tau$ if the observer is far from the black hole. In the limit $r = 2 M G$ we find
$
  dv(r, t) = 0
$
so the free-falling observer appears to stop moving. To see more we write
$
  dv(r, t) = - 1/r sqrt((2 M G)/r) (r - 2 M G)
$
and consider $r tilde.eq 2 M G$ then to first order in $r - 2 M G = dd(r, d: delta)$ we can just write
$
  dv(r, t) =^(r tilde 2 M G) - 1/(2 M G) (r- 2 M G)
$
which can be solved
$
  r = 2 M G + C - exp(- t/(2 M G))
$
so as $t -> oo$ we find $r -> 2 M G$, and the free-falling observer again appears to stop, but doing so in infinite time.

What happens to light? For light we have $dd(tau) = 0$ so
$
  (1- (2 M G)/r) dd(t^2) - dd(r^2)/(1-(2 M G)/r) = 0
$
the velocity of the light as seen by the far-away observer is then
$
  abs(dv(r, t)) = 1 - (2 M G)/r
$
for $r >= 2 M G$. So at $r = 2 M G$
$
  abs(dv(r, t)) = 0
$
so light can't escape the black hole---which is why $R$ is known as the horizon of the black hole. Also notice that $g_(00) -> 0$ as $r -> 2 M G$ meaning that a free-falling observer will appear redshifted as it crosses the horizon, and will be infinitely redshifted as it crosses---since the connection between $dd(tau)$ and $dd(t)$ is $g_(00)$.

=== A detour
We'll briefly discuss Hawking radiation since everyone should be aware of it, and it's just really cool.

Everywhere in vacuum we know from QFT that particle-antiparticle pairs can spontaneously appear, if they annihilate sufficiently fast---we can borrow energy from vacuum if we give it back by fast enough---this leads to a vacuum energy. We now consider what happens to these virtual pairs if they form at or near the horizon of a black hole. If one of the virtual particles is within the horizon while the other is outside then they'll never annihilate since the one within the horizon can't escape as we've just seen above. This means the virtual pair becomes real, and the escaping particle becomes radiation. This is only possible if the black hole loses energy, and since the black hole radiates it has a temperature. By Hawking this is given by
$
  k_B T = 1/(4pi) (hbar c)/R
$
with the temperature being called the Hawking temperature. Since black holes are classical objects they have energy $U_"BH" = M c^2$. From thermodynamics we have $T dd(S) = dd(U)$, which can be used to easily find (due to Bekenstein)
$
  dd(S) = (8 pi G k_B M)/(hbar c) dd(M) => S/k_B = (4 pi G M^2)/(hbar c) = 1/4 A/cal(l)_"pl"^2
$
with $A = 4 pi R^2$ being the area of the black hole. So black holes have entropy $prop$ area! We find that as the black hole loses energy and gets smaller then it gets hotter since $T prop R^(-1)$---which acts to increase the radiation---this is absurd! We also see that $S prop A$ even though usually $S prop V$---meaning that if we have some spherical mass which collapses then its entropy goes from being $prop V$ to being $prop A$ which is also absurd. This leads to the principle of holography, which essentially says that there are one-one correspondences between theories on the interior of some space-time to the boundary of said space-time. So everything is determined by the boundary---similar to Stokes theorem. Another thing we might consider is that we start with some pure state, nothing is hidden, then if we were to add a black hole to this pure state, then we could trace out the black hole (and the information) to get a mixed state. What then happens when the black hole eventually radiates away? Then we are back to the initial state but now we've just lost the information we initially traced out---this is the black hole information paradox.

== Kruskal-Szekeres coordinates
In the Schwarzschild metric we have two singularities the first is the singularity at $r = 0$, and the second is the one discussed above at the horizon $r = R$. As we've seen the second is weird in the sense that nothing really happens to the free-falling observer when it crosses the horizon. For this reason we might suspect that this is merely a coordinate artifact, which turns out to be correct. At the horizon there are really no problems with curvature etc. so we should be able to change coordinates to get rid of the singularity---unlike at $r = 0$ where curvature explodes.

The new coordinates are the Kruskal-Szekeres coordinates defined by---note in the actual lecture we used $X -> r'$ and $T -> t'$
$
  X^2 - T^2 = (r/R - 1) exp(r/R) " and " T/X = tanh (t/(2R))
$
where the exponential and hyperbolic tangent act to make our infinities finite. In these coordinates the Schwarzschild metric takes the form
$
  dd(tau^2) = (32 R^3)/(r^2) exp(- r/R) (dd(T^2) - dd(X^2)) - r^2 dd(Omega^2)
$
where $r equiv r(X,T)$. Now there are no problems at the horizon, but the singularity at $r = 0$ persists. Many things about these coordinates are nicely summarized by a Kruskal-Szekeres chart (in this case $R equiv 1$) as in kruskal:
/*
#figure(
  image("Kruskal_diagram_of_Schwarzschild_chart.svg", width: 100%),
  caption: [A Kruskal-Szekeres chart.],
) <kruskal>
*/
We see that for $t -> plus.minus oo => T = plus.minus X$ which corresponds to the two dashed lines. We also see that for $r = R => X^2 = T^2 => X = plus.minus T$ so the dashed line seperating I (exterior of black hole) and II (interior of black hole) is the horizon. The various lines drawn for constant $r$ are all hyperbolic due to $X^2 - T^2 = "const"$, with their direction depending on if $r < R$ or $r > R$. The limit is at $r = 0$ where $X^2 - T^2 = -1$ so for $X = 0 => T = 1$, meaning the singularity is a hyperbola with the minimum at $T = 1$. The two respective mirror regions are denoted III (mirror exterior) and IV (mirror interior or white hole). When $t=0 => T = 0 => X"-axis"$. To see how light moves we know that for light $ dd(tau^2)=0 =>^"no angular part" dd(X^2)=dd(T^2) $
so light moves along $45 degree$ lines. This is important because it means one can't travel from I $->$ III (they are causally disconnected), one can also never escape II (trap), but one is forced to escape IV (anti-trap). However, there seems to be a point connecting I and III, this is what we'll analyze next. But at $T=0$ region I and III appear connected?

=== Wormholes
We consider the slice $T=0$,
$
  X = [r/R-1]^(1\/2) exp(r/(2R))
$
then
$
  dd(T)=0";  " dd(X) = r^(1\/2)/(2 R^(3\/2)) (1-R/r)^(-1\/2) exp(r/(2 R)) dd(r)
$
We fix $theta = 0$, yielding the two-dimensional metric
$
  dd(s^2) = 1/(1-R/r) dd(r^2) + r^2 dd(phi^2)
$
this describes a two-dimensional surface (of revolution)! We want to embed this surface in three-dimensional Euclidean space. We do this by adding a coordinate $(r,phi) -> (r,phi,z)$, which are just cylindrical coordinates with the Euclidean metric
$
  dd(s^2) & = dd(z^2)+dd(r^2)+r^2 dd(phi^2) \
          & = ((dd(z)/dd(r))^2+1) dd(r^2) + r^2 dd(phi^2)
$
matching coefficients give a condition for $z$
$
  (dd(z)/dd(r))^2 +1 = 1/(1-R/r) => dd(z) = plus.minus sqrt(R/(r-R)) dd(r)
$
which can be solved
$
  z^2 (r) = 4 R(r-R)
$
plotting this gives an embedding diagram as in wormhole:
/*
#figure(
  image("Lorentzian_Wormhole.svg", width: 50%),
  caption: [An _embedding diagram_ of a wormhole.],
) <wormhole>
*/
At the throat: $z = 0 => r = R$, so the throat has radius $R$. So the surface at $T=0$ connects regions I and III, meaning in principle there is a path between the two black hole exteriors. We'd like to determine the lifetime of this path. Nothing is special to $T = 0$ about the above procedure, all slices with $T = "const"$ gives a surface with $r_"throat" = r_"min" = r(X=0, T)$. We saw that $T = 0 => r_"min" = R$. For $T = 1$ we find $r_"min" = 0$ meaning the wormhole closes, and for $T > 1$ it pinches off ($r_"min" in CC$). This happens on a timescale
$
  dd(t, d: Delta) tilde R/c
$
which is less then the time it takes for light to travel from I $->$ III.

#pagebreak()
=== A detour
We quickly derive the Kruskal-Szekeres coordinates in terms of the Schwarzschild coordinates.

Following the book we define
$
  p' equiv exp(p/(2 R))";  " q' equiv -exp(- q/(2 R))
$
where (Eddington-Finkelstein)
$
  p equiv c t + R ln abs(r-R) + r";  " q equiv c t - R ln abs(r-R) - r
$
then the Kruskal-Szekeres coordinates are given by
$
  p' = T+X
$
$
  q' = plus.minus (T - X)
$
with $+$ (outside $r > R$):
$
  T_+ = 1/2 (p' + q')";  " X_+ = 1/2 (p' - q')
$
with $-$ (inside $r < R$):
$
  q' = X - T
$
$
  X_- = 1/2 (p' + q') ";  " T_- = 1/2 (p' - q')
$
so the roles switch.


We compute
$
  T_+ &= 1/2 (exp(p/(2R))-exp(-q/(2 R))) \
  &= 1/2 [exp((c t + R ln abs(r-R)+r)/(2R))-exp((-c t + R ln abs(r-R)+r)/(2 R))] \
  &= 1/2 exp(r/(2 R)) exp((ln abs(r-R))/2) [exp((c t)/(2 R))-exp(-(c t)/(2 R))] \
  &= exp(r/(2 R)) sqrt(abs(r-R)) sinh((c t)/(2 R))
$
and
$
  X_+ & = exp(r/(2 R)) sqrt(abs(r-R)) cosh((c t)/(2 R))
$
for $r > R$:
$
  T_+ & = sqrt(r-R) exp(r/(2 R)) sinh ((c t)/(2 R)) \
  X_+ & = sqrt(r-R) exp(r/(2 R)) cosh((c t)/(2 R))
$
for $r < R$:
$
  T_- & = sqrt(R-r) exp(r/(2 R)) cosh ((c t)/(2 R)) \
  X_- & = sqrt(R-r) exp(r/(2 R)) sinh((c t)/(2 R))
$
the $cosh$ and $sinh$ switch since by definition
$
  T^2 - X^2 & = p' q' \
            & = - exp((p-q)/(2 R)) \
            & = - exp((r+ R ln abs(r-R))/(R)) \
            & = - exp(r/R) abs(r-R)
$
so for $r > R$
$
  T^2 - X^2 = (R-r) exp(r/R)
$
and for $r < R$
$
  T^2 - X^2 = -(R-r) exp(r/R)
$
which is satisfied by the above.
