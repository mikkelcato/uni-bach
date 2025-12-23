//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= The Einstein field equations
== The energy-momentum tensor
The Newtonian limit gave
$
  Phi = -(1+g_00)/2
$
Then
$
  nabla^2 Phi & = - 1/2 nabla^2 g_00 \
              & =^! 4 pi G_N rho
$
or
$
  nabla^2 g_00 (bold(x)) = - 8 pi G_N rho (bold(x))
$
with $rho(bold(x))$ being the matter density. We want to generalize this equation. We see by the LHS that $rho(bold(x))$ should generalize to some tensor. The natural choice is the energy-momentum tensor of some relativistic flow
$
  T^(mu nu) = rho (x) U^mu (x) U^nu (x)
$
where $rho(x)$ is the energy density seen by an observer moving with the flow.

To see that $T^(mu nu)$ does describe flow consider $T^(00)$
$
  T^(00) & = rho(x) (dv(x^0, tau))^2 \
         & = rho/(1-bold(v)^2)
$
where we use for $g_(mu nu) = eta_(mu nu)$
$
  dd(tau)^2 = (dd(x^0))^2 (1- bold(v)^2)
$
with $bold(v) = dd(bold(x))\/dd(x^0)$. Taking the Newtonian limit we have
$ T^(00) tilde.eq rho $
Then
$
  nabla^2 g_(00) = - 8 G_N T_(00)
$
The other components of $T^(mu nu)$ are
$
  T^(i j) & = underbracket((rho v_i v_j)/(1-bold(v)^2), "momentum" #linebreak() "current") \
  T^(i 0) & = underbracket((rho v_i)/(1 - bold(v)^2), "momentum" #linebreak() "density") = T^(0 i)
$
$T^(mu nu)$ is the energy-momentum tensor of a system (with $g_(mu nu) = eta_(mu nu)$) with no pressure in the rest-frame. With the above identifications.

Assuming the system is closed then $T^(mu nu)$ is conserved
$
  partial_nu T^(mu nu) = 0
$
Taking $mu = 0$ we find
$
  partial_nu T^(0 nu) &= pdv(, t) (rho/(1-bold(v)^2)) + nabla ((rho bold(v))/(1-bold(v)^2)) \ &=^! 0
$
which is the relativistic Euler equation. Taking $mu = i$ we find
$
  partial_nu T^(i nu) &= pdv(, t) ((rho v_i)/(1-bold(v)^2)) + pdv(, x) ((rho v_i v_x)/(1-bold(v)^2)) + pdv(, y) ((rho v_i v_y)/(1-bold(v)^2)) + pdv(, z) ((rho v_i v_z)/(1-bold(v)^2)) \
  &= rho/(1-bold(v)^2) (pdv(v^i, t) + bold(v) dot nabla v^i) \
  &=^!_("four-velocity constant") 0
$
which is the relativistic Navier-Stokes equation. Then $partial_nu T^(mu nu) = 0$ is a direct consequence of relativistic fluid mechanics.

== The perfect fluid
A perfect fluid includes pressure. We define a perfect fluid by a field $bold(v) (x)$ with the property that any observer moving with the fluid sees it as isotropic everywhere. The rest frame requires
$
     T^00 & = rho \
  T^(i 0) & = 0 = T^(0 i) \
  T^(i j) & = p delta_(i j) = T^(j i)
$
Assuming the fluid has some velocity $bold(v) (x)$ then
$
  T^(mu nu) = p eta^(mu nu) + (p + rho) U^mu U^nu
$
This reduces to the rest-frame $T^(mu nu)$. The generally covariant form is obtained by $eta^(mu nu) -> g^(mu nu)$
$
  T^(mu nu) = p g^(mu nu) + (p + rho) U^mu U^nu
$
which is very important!

The fluid equations in general relativity are
$
  D_nu T^(mu nu) = 0
$
but the total energy-momentum defined by
$
  P^mu = integral dd(x, 3) sqrt(-g) T^(mu 0)
$
is not conserved!

== The curvature tensor
We want some covariant object that vanishes in flat space but not in curved space. This would let us determine if space is flat independent of coordinates. We already have a candidate since
$
  [partial_mu, partial_nu] =^"flat space" 0
$
meaning if
$
  [D_mu, D_nu] eq.not^"curved space" 0
$
we would be done since $D_mu -> partial_mu$ in flat space.

We now show $[D_mu, D_nu] eq.not 0$. By definition
$
  D_kappa T_(mu nu) = partial_kappa T_(mu nu) - tensor(Gamma, lambda, -nu kappa) T_(mu lambda) - tensor(Gamma, lambda, -mu kappa) T_(lambda nu)
$
Let $T_(mu nu) equiv D_mu V_nu$. Then
$
  D_(kappa) D_mu V_nu = partial_kappa D_mu V_nu - tensor(Gamma, lambda, -nu kappa) D_mu V_lambda - tensor(Gamma, lambda, -mu kappa) D_lambda V_nu
$
and similarly
$
  D_mu D_kappa V_nu = partial_mu D_kappa V_nu - tensor(Gamma, lambda, -nu mu) D_kappa V_lambda - tensor(Gamma, lambda, -kappa mu) D_lambda V_nu
$
Then
$
  [D_kappa, D_mu] V_nu &= D_kappa D_mu V_nu - D_mu D_kappa V_nu \ &= partial_kappa D_mu V_nu - partial_mu D_kappa V_nu - tensor(Gamma, lambda, -nu kappa) D_mu V_lambda + tensor(Gamma, lambda, -nu mu) D_kappa V_lambda \
  &= partial_kappa (partial_mu V_nu - tensor(Gamma, sigma, -mu nu) V_sigma) - partial_mu (partial_kappa V_nu - tensor(Gamma, sigma, -kappa nu) V_sigma) \ &#h(10pt)- tensor(Gamma, lambda, -nu kappa) (partial_mu V_lambda - tensor(Gamma, sigma, -mu lambda) V_sigma) + tensor(Gamma, lambda, -nu mu) (partial_kappa V_lambda - tensor(Gamma, sigma, -kappa lambda) V_sigma) \
  &= - partial_kappa tensor(Gamma, sigma, -mu nu) V_sigma - tensor(Gamma, sigma, -mu nu) partial_kappa V_sigma + partial_mu tensor(Gamma, sigma, -kappa nu)V_sigma + tensor(Gamma, sigma, -kappa nu) partial_mu V_sigma \
  &#h(10pt)- tensor(Gamma, lambda, -nu kappa) partial_mu V_lambda + tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda) V_sigma + tensor(Gamma, lambda, -nu mu) partial_kappa V_lambda - tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda) V_sigma \
  &= [-partial_kappa tensor(Gamma, sigma, -mu nu) + partial_mu tensor(Gamma, sigma, -kappa nu) - tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda)+ tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda)] V_sigma \
  &= - tensor(R, sigma, -nu mu kappa) V_sigma
$
where we have defined the _Riemann-Christoffel curvature tensor_
$
  tensor(R, sigma, -nu mu kappa) = partial_kappa tensor(Gamma, sigma, -mu nu) - partial_mu tensor(Gamma, sigma, -kappa nu) + tensor(Gamma, lambda, -nu mu) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -nu kappa) tensor(Gamma, sigma, -mu lambda) =^"flat space" 0
$

A more useful form is
$
  R_(lambda mu nu kappa) &= g_(lambda sigma) tensor(R, sigma, -nu mu kappa) \
  &= 1/2 [pdv(g_(lambda nu), x^kappa, x^mu) +pdv(g_(mu kappa), x^nu, x^lambda)-pdv(g_(mu nu), x^kappa, x^lambda) - pdv(g_(lambda kappa), x^nu, x^mu)] \ &#h(8.4em)+ g_(eta sigma) [tensor(Gamma, eta, -nu lambda) tensor(Gamma, sigma, -mu kappa) - tensor(Gamma, eta, -kappa lambda) tensor(Gamma, sigma, -mu nu)]
$
so $R tilde partial^2 g$ as we would expect. This form also shows the symmetries
$
  R_(lambda mu nu kappa) &= R_(nu kappa lambda mu) "     " (lambda mu)(nu kappa) <-> (nu kappa) (lambda mu) \
  R_(lambda mu nu kappa) &= - R_(mu lambda nu kappa) "   " lambda <-> mu \
  R_(mu lambda nu kappa) &= R_(lambda mu kappa nu) "     " lambda <-> mu, kappa <-> nu \
  R_(lambda mu nu kappa) &= - R_(lambda mu kappa nu) "   " kappa <-> nu \
  0 &= underbracket(R_(lambda mu nu kappa) + R_(lambda kappa mu nu) + R_(lambda nu kappa mu), "cyclic permutation of" (mu nu kappa))
$
with the last being the _first Bianchi identity_. To see these it is simplest to consider a freely falling elevator.

We define the _Ricci tensor_
$
  R_(mu kappa) = g^(lambda nu) R_(lambda mu nu kappa) = tensor(R, lambda, -mu lambda kappa)
$
which is symmetric. We also define the _Ricci scalar_
$
  R = g^(mu kappa) R_(mu kappa) = tensor(R, mu, -mu)
$
We need the _second Bianchi identity_
$
  D_eta R_(lambda mu nu kappa) + D_kappa R_(lambda mu eta nu) + D_nu R_(lambda mu kappa eta) = 0
$
and recall
$
  D_sigma g_(mu nu) = 0
$
Then we can rewrite the second Bianchi identity
$
  0&= g^(lambda nu) (D_eta R_(lambda mu nu kappa) + D_kappa R_(lambda mu eta nu) + D_nu R_(lambda mu kappa eta) ) \ 0&= D_eta R_(mu kappa) - D_kappa R_(mu eta) + D_nu tensor(R, nu, -mu kappa eta) \
  0 &= g^(mu kappa) (D_eta R_(mu kappa) - D_kappa R_(mu eta) + D_nu tensor(R, nu, -mu kappa eta) ) \
  0 &= D_eta R - D_mu tensor(R, mu, -eta) - D_nu tensor(R, nu, -eta) \
  0 &= D_mu (tensor(R, mu, -eta) - 1/2 tensor(delta, mu, -eta) R) \
  0 &= D_mu underbracket((R^(mu nu) - 1/2 g^(mu nu) R), "trace-reversed Ricci tensor")
$
So the _trace-reversed Ricci tensor_ is conserved!

== The field equations
We had
$
  nabla^2 g_(00) (x) = - 8 pi G_N T_(0 0) (x)
$
We write this as a tensor equation
$
  G_(mu nu) (x) = - 8 pi G_N T_(mu nu) (x)
$
where $G_(mu nu)$ depends on $g$ and $partial g$. By dimensional analysis the LHS should have units $[L^(-2)]$ which matches $partial^2 g$. Then an obvious guess would be
$
  G_(mu nu) eq^? R_(mu nu)
$
However we require
$
  D_nu G^(mu nu) =^! 0 " since " D_nu T^(mu nu) = 0
$
implying
$
  G_(mu nu) = A(R_(mu nu) - 1/2 g_(mu nu) R)
$

To determine $A$ consider
$
  G_(00) & = A(R_(00) + 1/2 R) = - 8 pi G_N T_(00)
$
We take the Newtonian limit. Then $T_(i j) -> 0$ giving
$
        0 & tilde.eq R_(i j) - 1/2 g_(i j) R \
  R_(i j) & tilde.eq 1/2 g_(i j) R
$
And $g_(i j) tilde.eq n_(i j)$ giving
$
  R & tilde.eq R_(i i) - R_(00) tilde.eq 3/2 R - R_(0 0) \
  R & tilde.eq 2 R_(00)
$
Then
$
  G_(00) = 2 A R_(00)
$
Similarly $Gamma Gamma tilde 0$ giving
$
  R_(lambda mu nu kappa) tilde.eq 1/2 (pdv(g_(lambda nu), x^kappa, x^mu) - pdv(g_(mu nu), x^nu, x^lambda) - pdv(g_(lambda kappa), x^nu, x^mu) + pdv(g_(mu kappa), x^nu, x^lambda))
$
$g$ is $tilde$ static so
$
     R_(0000) & tilde.eq 0 \
  R_(i 0 j 0) & tilde.eq 1/2 pdv(g_(00), x^i, x^j)
$
Then
$
  G_(00) & tilde.eq 2 A R_(0 0) \
         & tilde.eq 2 A (R_(i 0 i 0) - R_(0 0 0 0)) \
         & tilde.eq A nabla^2 g_(0 0)
$
implying $A =^! 1$.

The _Einstein field equations_ are then
$
  underbracket(R_(mu nu) - 1/2 g_(mu nu) R, G_(mu nu)) = - 8 pi G_N T_(mu nu)
$
where $G_(mu nu)$ is typically called the _Einstein tensor_. A more useful form is
$
  R_(mu nu) = - 8 pi G_N (T_(mu nu) - 1/2 g_(mu nu) tensor(T, alpha, -alpha))
$
We see $T_(mu nu) = 0$ gives $R_(mu nu) = 0$.

The EFE can be extended by including a _cosmological constant_ $Lambda$
$
  R_(mu nu) - 1/2 g_(mu nu) R - Lambda g_(mu nu) = - 8 pi G_N T_(mu nu)
$
since $D_mu g_(mu nu) = 0$! We can interpret $Lambda g_(mu nu)$ as an energy density since we can move it to the RHS and absorb it in $T_(mu nu)$.

= The Einstein-Hilbert action
== Without matter
The simplest action we can write is
$
  S = integral dd(x, 4) sqrt(-g) R
$
this is the _Einstein-Hilbert action_. We write $R = g^(mu nu) R_(mu nu)$ to find
$
  dd(S, d: delta) = integral dd(x, 4) [(dd(sqrt(-g), d: delta))g^(mu nu) R_(mu nu) + sqrt(-g) (dd(g^(mu nu), d: delta))R_(mu nu) + sqrt(-g) g^(mu nu) dd(R_(mu nu), d: delta)]
$
We claim
$
  dd(g^(mu nu), d: delta) &= - g^(mu rho) g^(nu sigma) dd(g_(rho sigma), d: delta) \
  dd(sqrt(-g), d: delta) &= - 1/2 sqrt(-g) g_(mu nu) dd(g^(mu nu), d: delta) \
  dd(R_(mu nu), d: delta) &= D_rho dd(tensor(Gamma, rho, -mu nu), d: delta) - D_nu dd(tensor(Gamma, rho, -mu rho), d: delta)
$
with
$
  dd(tensor(Gamma, rho, -mu nu), d: delta) = 1/2 g^(rho sigma) (D_mu dd(g_(sigma nu), d: delta) + D_nu dd(g_(sigma mu), d: delta)- D_sigma dd(g_(mu nu), d: delta))
$
and $delta Gamma$ being a tensor since it is the difference of two Christoffel symbols. Consider
$
  dd(tensor(Gamma, rho, -mu nu), d: delta) &=^"trivially" 1/2 g^(rho sigma) (partial_mu dd(g_(sigma nu), d: delta) + partial_nu dd(g_(sigma mu), d: delta)-partial_sigma dd(g_(mu nu), d: delta))
$
and since $delta Gamma$ is a tensor it is valid for $partial_mu -> D_mu$. Consider
$
  tensor(R, sigma, -rho mu nu) =^"free fall" partial_mu tensor(Gamma, sigma, -nu rho) - partial_nu tensor(Gamma, sigma, -mu rho)
$
Then
$
  dd(tensor(R, sigma, -rho mu nu), d: delta) &= partial_mu dd(tensor(Gamma, sigma, -nu rho), d: delta) - partial_nu dd(tensor(Gamma, sigma, -mu rho), d: delta) \
  &=^"everywhere" D_mu dd(tensor(Gamma, sigma, -nu rho), d: delta) - D_nu dd(tensor(Gamma, sigma, -mu rho), d: delta)
$
for $mu = sigma$
$
  dd(R_(rho nu), d: delta) = D_mu dd(tensor(Gamma, mu, -nu rho), d: delta)-D_nu dd(tensor(Gamma, mu, -rho mu), d: delta)
$
implying
$
  g^(mu nu) dd(R_(mu nu), d: delta) = D_mu X^mu
$
with
$
  X^mu = g^(rho nu) dd(tensor(Gamma, mu, -rho nu), d: delta) - g^(mu nu) dd(tensor(Gamma, rho, -nu rho), d: delta)
$
Then
$
  dd(S, d: delta) = integral dd(x, 4) sqrt(-g) [(R_(mu nu)-1/2 R g_(mu nu))dd(g^(mu nu), d: delta) + underbracket(D_mu X^mu, "total derivative")]
$
Taking $dd(S, d: delta) = 0$ we find
$
  G_(mu nu) equiv R_(mu nu) - 1/2 R g_(mu nu) = 0
$
which are the vacuum EFE.

A slightly less trivial action is found by multiplying the volume form with a constant
$
  S = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (R - 2 Lambda)
$
Taking $dd(S, d: delta) = 0$ we find
$
  R_(mu nu) - 1/2 R g_(mu nu) = - Lambda g_(mu nu)
$
which are the vacuum EFE with a cosmological constant.

== With matter
The action becomes
$
  S = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (R-2 Lambda) + underbrace(S_M, "matter")
$
and we define
$
  T_(mu nu) = - 2/sqrt(-g) dv(S_M, g^(mu nu), d: delta)
$
Then
$
  dd(S, d: delta) = 1/(16 pi G) integral dd(x, 4) sqrt(-g) (G_(mu nu)+Lambda g_(mu nu)) dd(g^(mu nu), d: delta) - 1/2 integral dd(x, 4) sqrt(-g) T_(mu nu) dd(g^(mu nu), d: delta)
$
Taking $dd(S, d: delta)= 0$ we find
$
  G_(mu nu) + Lambda g_(mu nu) = 8 pi G T_(mu nu)
$
which are the EFE.

= The geodesic as the minimal curve
== The geodesic equation
Consider the line element $dd(s^2) = g_(mu nu) dd(x^mu) dd(x^nu)$. Then
$
  S & = integral dd(s) = integral L dd(tau)
$
where
$
  L equiv sqrt(g_(mu nu) dv(x^mu, tau) dv(x^nu, tau))
$
We parametrize a curve by $x^mu (tau)$. Then the shortest curve between $x_i^mu$ and $x_f^mu$ is found by $dd(S, d: delta) = 0$. We use the _Euler-Lagrange equations_
$
  pdv(L, x^mu) - dv(, tau) pdv(L, dot(x)^mu) = 0
$
The equations of motion from $tilde(L) = L^2$ are the same as those obtained by using $L$. We let $L equiv g_(mu nu) dot(x)^mu dot(x)^nu$. Then
$
  pdv(L, dot(x)^mu) & = 2 g_(mu nu) dot(x)^nu \
       pdv(L, x^mu) & = pdv(g_(nu lambda), x^mu) dot(x)^nu dot(x)^lambda
$
and we find
$
  0 &= pdv(g_(mu nu), x^lambda) dv(x^lambda, tau) dv(x^nu, tau) + g_(mu nu) dv(x^nu, tau, 2) - 1/2 pdv(g_(nu lambda), x^mu) dv(x^nu, tau) dv(x^lambda, tau) \
  0 &=^(g^(sigma mu) (dots)) dv(x^sigma, tau, 2) + 1/2 g^(sigma mu) (2 pdv(g_(mu nu), x^lambda) - pdv(g_(nu lambda), x^mu)) dv(x^lambda, tau) dv(x^nu, tau) \
  0 &= dv(x^sigma, tau, 2) + tensor(Gamma, sigma, -lambda nu) dv(x^lambda, tau) dv(x^nu, tau)
$
which is the geodesic equation!

== A trick
We can easily find $Gamma$ by comparing the integrand of
$
  dd(integral g_(mu nu) dv(x^mu, tau) dv(x^nu, tau) dd(tau), d: delta) = 0
$
with the geodesic equation.

= The time-dependent spherically symmetric metric
== The metric
$dd(tau^2)$ can only depend on rotationally invariant quantities. These are
$
  {t, dd(t), r, r dd(r) = bold(x) dd(bold(x)), dd(r^2) + r^2 (dd(theta^2) + sin^2 theta dd(phi^2)) = dd(bold(x)^2)}
$
We define $dd(Omega^2) equiv dd(theta^2) + sin^2 theta dd(phi^2)$. Then the metric has the form
$
  dd(tau^2) = A dd(t^2) - B dd(r^2) - C dd(r, t) - D r^2 dd(Omega^2)
$
with $A, B, dots$ being functions of $t$ and $r$. We absorb $D$ into $r$ by redefining $r -> r' = r sqrt(D)$
$
  dd(tau^2) = A dd(t^2) - B dd(r^2) - C dd(r, t) - r^2 dd(Omega^2)
$
with new $A, B, dots$ and $C$. Similarly to absorb $C$ we redefine $t -> t'$ by
$
  dd(t') = underbracket(eta(r, t), "independent") [A dd(t)-1/2 C dd(r)]
$
Then
$
  1/(A eta^2) dd(t'^2) = A dd(t^2)- C dd(t, r) + C^2/(4 A) dd(r^2)
$
implying
$
  dd(tau^2) & = 1/(eta^2 A) dd(t^2) - (B+ C^2/(4 A)) dd(r^2) - r^2 dd(Omega^2) \
  dd(tau^2) & = E(r, t) dd(t^2) - F(r,t) dd(r^2) - r^2 dd(Omega^2)
$
We find
$
  g_(r r) &= F";  " g_(theta theta) = r^2";  " g_(phi phi) = r^2 sin^2 theta";  " g_(t t) = - E \
  g^(r r) &= 1/F";  " g^(theta theta) = 1/r^2";  " g^(phi phi) = 1/(r^2 sin^2 theta)";  " g^(t t) = -1/E
$


== The Christoffel symbols
We compute the Christoffel symbols using
$
  dd(integral dd(tau) [E dot(t)^2 - F dot(r)^2 - r^2 dot(theta)^2 - r^2 sin^2 theta dot(phi)^2], d: delta) = 0
$
Consider the $mu = 0$ component of
$
  pdv(L, x^mu) = dv(, tau) pdv(L, dot(x)^mu)
$
We find
$
  dot(t)^2 pdv(E, t) - pdv(F, t) dot(r)^2 &= dv(, tau) (2 E dot(t)) \
  &= 2 E dot.double(t) + 2 dot(t) dv(E, tau) \
  &= 2 E dot.double(t) + 2 dot(t) (pdv(t, tau) pdv(E, t) + pdv(r, tau) pdv(E, r)) \
  &= 2 E dot.double(t) + 2 dot(t)^2 pdv(E, t) + 2 dot(t) dot(r) pdv(E, r)
$
simplifying
$
  0 &= pdv(F, t) dot(r)^2 + 2 E dot.double(t) + dot(t)^2 pdv(E, t) + 2 dot(t) dot(r) pdv(E, r) \
  0 &= dot.double(t) + 1/(2 E) pdv(E, t) dot(t)^2 + 1/E pdv(E, r) dot(t) dot(r) + 1/(2 E) pdv(F, t) dot(r)^2
$
We compare with
$
  0 = dot.double(t) + tensor(Gamma, 0, -mu nu) dot(x)^mu dot(x)^nu
$
giving
$
  tensor(Gamma, t, -t t) & = 1/(2 E) pdv(E, t) \
  tensor(Gamma, t, -r r) & = 1/(2 E) pdv(F, t) \
  tensor(Gamma, t, -r t) & = tensor(Gamma, t, -t r) = underbracket(1/2 times, "symmetric") 1/(E) pdv(E, r)
$
Similarly the $r$-component gives
$
          tensor(Gamma, r, -t r) & = tensor(Gamma, r, -r t) = 1/(2 F) pdv(F, t) \
          tensor(Gamma, r, -r r) & = 1/(2 F) pdv(F, r) \
          tensor(Gamma, r, -t t) & = 1/(2 F) pdv(E, r) \
  tensor(Gamma, r, -theta theta) & = - r/F \
      tensor(Gamma, r, -phi phi) & = - (r sin^2 theta)/F
$
And the $theta$-component gives
$
  tensor(Gamma, theta, -r theta) & = tensor(Gamma, theta, -theta r) = 1/r \
  tensor(Gamma, theta, -phi phi) & = - sin theta cos theta
$
And the $phi$-component gives
$
  tensor(Gamma, phi, -r phi) &= tensor(Gamma, phi, -phi r) = 1/r \
  tensor(Gamma, phi, -theta phi) &= tensor(Gamma, phi, -phi theta) = (cos theta)/(sin theta)
$
== The Ricci tensor
We have
$
  tensor(Gamma, mu, -mu nu) = 1/2 g^(mu sigma) pdv(g_(mu sigma), x^nu) = 1/2 pdv(ln g, x^nu)
$
so
$
  R_(mu kappa) &= partial_kappa tensor(Gamma, sigma, -sigma mu) - partial_sigma tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu sigma) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -mu kappa) tensor(Gamma, sigma, -sigma lambda) \
  &= 1/2 partial_kappa partial_mu ln g - partial_sigma tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu sigma) tensor(Gamma, sigma, -kappa lambda) - 1/2 tensor(Gamma, lambda, -mu kappa) partial_lambda ln g
$
Then
$
  R_(r r) &= 1/(2 E) pdv(E, r, 2)- 1/(4 E^2) (pdv(E, r))^2 - 1/(4 E F) pdv(E, r) pdv(F, r) \
  &#h(12pt)- 1/(r F) pdv(F, r) - 1/(2 E) pdv(F, t, 2) + 1/(4 E^2) pdv(E, t) pdv(F, t) + 1/(4 E F) (pdv(F, t))^2 \
  R_(theta theta) &= -1 + 1/F - r/(2 F^2) pdv(F, r) + r/(2 E F) pdv(E, r) \
  R_(t t) &= -1/(2 F) pdv(E, r, 2) + 1/(4 F^2) pdv(E, r) pdv(F, r) - 1/(r F) pdv(E, r) \
  &#h(12pt)+ 1/(4 E F) (pdv(E, r))^2 + 1/(2 F) pdv(F, t, 2) - 1/(4 F^2) (pdv(F, t))^2 - 1/(4 E F) pdv(E, t) pdv(F, t) \
  R_(t r) &=^"simplified" -1/(r F) pdv(F, t) \
  R_(phi phi) &= sin^2 theta R_(theta theta)
$
The trace-inverted EFE then give equations for $E$ and $F$ in terms of $T_(mu nu)$. As an example since $g_(t r) = 0$ we have
$ 1/(r F) pdv(F, t) = 8 pi G T_(t r) $

= The Schwarzschild solution
== The solution
Consider a point mass with mass $M$ at the origin. We assume $m_"obs" << M$ so we can ignore any backreaction. We use the time-dependent spherically symmetric metric
$ dd(tau^2)=E dd(t^2) - F dd(r^2) - r^2 dd(Omega^2) $
Since $T_(mu nu) = 0$ for all $r eq.not 0$ we have $R_(mu nu) = 0$. Using $R_(r r) = R_(t t) = 0$ we find
$
  0 &= R_(r r)/F + R_(t t)/E \
  0 &=^"time-independent" 1/(r F) underbracket((1/F pdv(F, r) + 1/E pdv(E, r)), =^!0)
$
implying
$
  pdv(ln F, r) = - pdv(ln E, r)
$
so $E(r) F(r) = "const"$. As $r->oo$ the metric should reduce to $eta_(mu nu)$ where $E = F = 1$. Then
$
  E = 1/F
$
for all $r$. Consider
$
       R_(theta theta) & = 0 \
  -1 + E + r pdv(E, r) & = 0
$
implying
$
  pdv(, r) (r E) = 1
$
Then
$
  E = 1 + C/r
$
We take the Newtonian limit
$
  g_(0 0) & = - (1 + 2 Phi) \
          & =^"point mass" - 1 + (2 G_N M)/r
$
Using $E = - g_(0 0)$ we find $C = -2 G_N M$ so
$
  E = 1 - (2 M G_N)/r";   "F = 1/E
$
We have found _the Schwarzschild metric_
$
  dd(tau^2) = (1 - (2 M G_N)/r) dd(t^2) - dd(r^2)/(1-(2 M G_N)/r) - r^2 dd(Omega^2)
$
which is very important! Due to _Birkhoff's theorem_ this describes any spherically symmetric mass distribution.

== Birkhoff's theorem
#thm[Birkhoff's theorem][
  A spherically symmetric gravitational field in empty space must be static with a metric given by the Schwarzschild metric.
]
#proof[
  We assumed $E$ and $F$ were time-independent. This can be shown to always be the case. From $R_(t r) = 0$ we have
  $
    pdv(F, t) = 0
  $
  which is nice. We can write
  $
    E = f(t) (1-(2 M G_N)/r) eq.not E(r)
  $
  and redefine the time-coordinate to absorb $f(t)$ by
  $
    t' = integral sqrt(f(t)) dd(t) => dd(t') = sqrt(f(t)) dd(t)
  $
  So $E$ and $F$ are always time-independent! Then everything done before is valid when $T_(mu nu) = 0$. We can even have weird oscillatory mass distributions if they are spherically symmetric.
]

This explains why gravitational waves are rare since we need to break spherical symmetry for the metric to be changing in time. This also makes it trivial to analyze cavities since these are spherically symmetric and $T_(mu nu) = 0$. But $M = 0$ so in cavities we just have $eta_(mu nu)$!

== Killing vectors
Consider the Lagrangian (with $theta = pi\/2$)
$
  L = g_(mu nu) dot(x)^mu dot(x)^nu = - (1-R/r) dot(t)^2 + 1/(1-R/r) dot(r)^2 + r^2 dot(phi)^2
$
where we define the _Schwarzchild radius_ $R equiv 2 M G_N$.

The Schwarzschild metric is time-independent. Generally if $g_(mu nu)$ is independent of $x^0$ we have
$
  pdv(L, x^0) = 0
$
Then
$
  dv(, tau) underbracket((pdv(L, dot(x)^0)), "constant") =^! 0
$
so
$
  pdv(L, dot(x)^0) = pdv(L, dot(t)) = 2 g_(mu 0) dot(x)^mu
$
is conserved. To ensure covariance we define the _time-like Killing vector_
$
  K^mu_((t)) equiv vecrow(1, 0, 0, 0)
$
Then
$
  gamma_((t)) = g_(mu nu) dot(x)^mu K^nu_((t))
$
is covariant and conserved. For other Killing vectors $K_((alpha))^nu$ the object
$
  gamma_((alpha)) = g_(mu nu) dot(x)^mu K^nu_((alpha))
$
is covariant and conserved if $g_(mu nu)$ is independent of $x^alpha$.

The Schwarzchild metric only depends on $r$ so we have conserved quantities related to both $t$ and $phi$ defined by
$
  kappa &equiv - g_(mu nu) dot(x)^mu K^nu_((t)) = - g_(00) dot(t) = (1- R/r) dot(t) \
  l &equiv m g_(mu nu) dot(x)^mu K^nu_((phi)) = m g_(phi phi) dot(phi) = m r^2 phi
$
We recognize $kappa$ as an _energy_ and $l$ as the angular momentum.

== Precession of orbits
We use the expressions for $kappa$ and $l$ to rewrite $L$. With $L = -1$ we find
$
  - kappa^2/(1-R/r) + dot(r)^2/(1-R/r) + l^2/(m^2 r^2) = - 1
$
We multiply by $1/2 m (1 - R/r)$ and define
$
  E/m equiv (kappa^2 -1)/2
$
Then
$
  underbracket(1/2 m r^2, E_"kin") + underbracket((1- R/r), "GR corr.") underbracket(l^2/(2 m r^2), E_l) - underbracket((G_N m M)/r, E_g) = E
$
This expresses conservation of energy with a relativistic correction. We can solve this to find the equation for an orbit
$
  r = alpha/(1 + e cos[(1-epsilon.alt)phi])
$
with $e$ being the eccentricity and
$
        alpha & = l^2/(G M m^2) = (1 + e) r_"min" \
  epsilon.alt & = (3 R)/(2 alpha)
$
The orbit returns to $r_"min"$ at
$
  phi = (2 pi)/(1 - epsilon.alt) tilde.eq 2 pi + (3 pi R)/alpha
$
so the orbit precesses! Einstein used this to compute the precession of Mercury's orbit which helped show the _correctness_ of general relativity.

= Black holes
We write the Schwarzschild metric as
$
  dd(tau^2) = f(r) dd(t^2) - dd(r^2)/f(r) - r^2 dd(Omega^2)
$
where we define
$
  f(r) equiv 1 - R/r
$
At $r = R$ we find a singularity. Usually this is not a problem since $R < R_"mass"$. But what about a point mass or some object with $R_"mass" < R$? We call such objects black holes for reasons we discuss now.

== The geodesic equation
We consider a free-falling observer. Then
$
  dv(x^alpha, tau, 2) + tensor(Gamma, alpha, -mu nu) dv(x^mu, tau) dv(x^nu, tau) = 0
$
or
$
  0 & = dd(integral [g_(mu nu) dot(x)^mu dot(x)^nu] dd(tau), d: delta) \
  & = dd(integral underbracket([-f(r) dot(t)^2 + f^(-1) dot(r)^2], L), d: delta) dd(tau)
$
where we assume $dot(theta)=dot(phi)=0$. This only depends on $r$ so we have the conserved quantity
$
  C = g_(mu nu) dot(x)^mu K^nu_((t)) = f(r) dot(t)
$
We fix $C$ by considering $r -> oo$ where $dot(t) -> 1$ and $f(r) -> 1$ so $C = 1$ meaning $dot(t) = f^(-1) (r)$. Then
$
  L = - 1/f(r) + f^(-1) (r) dot(r)^2 = - 1
$
implying
$
  1/f(r) - f^(-1) (r) dot(r)^2 =1
$
Then
$
  dot(r)^2 & = 1 - f(r) \
           & = R/r
$
which is nice and simple!

== Crossing the horizon
We want to know what happens when crossing $r = R$. By the above
$
  sqrt(r) dot(r) = - sqrt(2 M G)
$
which has solution
$
  2/(3 sqrt(2 M G)) (r^(3\/2)-r_0^(3\/2)) = tau_0 - tau
$
We see that nothing special at $r = R$. We even find for $r = 0$
$
  tau = tau_0 + 2/(3 sqrt(2 M G)) r_0^(3\/2) < oo
$
So the free-falling observer reaches the center in finite time! To find $t$ we do
$
  dv(r, t) & = dot(r)/dot(t) \
           & = - sqrt((2 M G)/r) (1 - (2 M G)/r)
$
Assuming $r >> 2 M G$ we find
$
  dv(r, t) = - sqrt((2 M G)/r)
$
so $t tilde tau$. Assuming $r = 2 M G$ we find
$
  dv(r, t) = 0
$
so the free-falling observer appears to stop moving. We write
$
  dv(r, t) = - 1/r sqrt((2 M G)/r) (r - 2 M G)
$
and assume $r tilde.eq 2 M G$. Then to first order in $r - 2 M G = dd(r, d: delta)$ we have
$
  dv(r, t) = - 1/(2 M G) (r- 2 M G)
$
Then
$
  r = 2 M G + C - exp(- t/(2 M G))
$
so as $t -> oo$ we find $r -> 2 M G$. The free-falling observer appears to stop at $r = R$ but taking infinite time to do so!

For light $dd(tau) = 0$
$
  (1- (2 M G)/r)^2 dd(t^2) - dd(r^2) = 0
$
The velocity of light as seen by some far-away observer is then
$
  abs(dv(r, t)) = 1 - (2 M G)/r
$
for $r >= 2 M G$. At $r = 2 M G$
$
  abs(dv(r, t)) = 0
$
so light can not escape the black hole! This is why we call $R$ the _horizon_. Also $g_(00) -> 0$ as $r -> 2 M G$ meaning a free-falling observer will appear more and more redshifted.

== Kruskal-Szekeres coordinates
The singularity at $r = R$ is a coordinate artifact. The only _true_ singularity is at $r = 0$ where curvature explodes to $oo$. We can use _Kruskal-Szekeres coordinates_ defined by
$
  X^2 - T^2 = (r/R - 1) exp(r/R) " and " T/X = tanh (t/(2R))
$
to _remove_ the singularity at $r = R$. With these the Schwarzschild metric takes the form
$
  dd(tau^2) = (32 R^3)/(r^2) exp(- r/R) underbracket((dd(T^2) - dd(X^2)), "light travels at" 45 degree) - r^2 dd(Omega^2)
$
where $r equiv r(X,T)$. Things about these coordinates are nicely summarized by a _Kruskal-Szekeres chart_ as in @kruskal:

#figure(
  image("Kruskal_diagram_of_Schwarzschild_chart.svg", width: 100%),
  caption: [A Kruskal-Szekeres chart with $R = 1$.],
) <kruskal>

We see as $t -> plus.minus oo$ we have $T = plus.minus X$ corresponding to the dashed lines. We see that for $r = R$ we have $X^2 = T^2$ so again $X = plus.minus T$. The dashed line seperating I and II is then the horizon. The various lines drawn for constant $r$ are all hyperbolic since $X^2 - T^2 = "constant"$ with their direction depending on if $r < R$ or $r > R$. The limit is $r = 0$ where $X^2 - T^2 = -1$. Then for $X = 0$ we have $T = 1$ meaning the singularity is a hyperbola with the minimum at $T = 1$. The two respective mirror regions are denoted by III and IV. When $t=0$ then $T = 0$ so this corresponds to the $X$-axis. Since light travels at $45 degree$ it is not possible to travel from I to III. Similarly light can never escape II or enter IV.

But when $T=0$ region I and III appear connected? This is what we call a _wormhole_.

== Wormholes
Consider the slice $T=0$
$
  X = (r/R-1)^(1\/2) exp(r/(2R))
$
Then
$
  dd(T)=0";  " dd(X) = r^(1\/2)/(2 R^(3\/2)) (1-R/r)^(-1\/2) exp(r/(2 R)) dd(r)
$
Let $theta = 0$ then
$
  dd(s^2) = (1-R/r)^(-1) dd(r^2) + r^2 dd(phi^2)
$
this describes a two-dimensional surface!

We want to embed this surface in $RR^3$. We do this by adding a coordinate $(r,phi) -> (r,phi,z)$. These are just cylindrical coordinates with the metric
$
  dd(s^2) & = dd(z^2)+dd(r^2)+r^2 dd(phi^2) \
          & = [(dd(z)/dd(r))^2+1] dd(r^2) + r^2 dd(phi^2)
$
By comparison we find
$
  (dd(z)/dd(r))^2 +1 = (1-R/r)^(-1)
$
implying
$
  dd(z) = plus.minus sqrt(R/(r-R)) dd(r)
$
Then
$
  z^2 (r) = 4 R(r-R)
$
plotting this gives an _embedding diagram_ as in @wormhole:

#figure(
  image("Lorentzian_Wormhole.svg", width: 50%),
  caption: [The embedding diagram of a wormhole.],
) <wormhole>

At the throat $z = 0$ where $r = R$. So the surface at $T=0$ connects regions I and III!

We would like to determine the lifetime of the wormhole. The embedding can be done for $T eq.not 0$ in the same way.  All slices $T = "constant"$ give a surface with $ r_"throat" = r_"min" = r(X=0,T) $
For $T = 1$ we find $r_"min" = 0$ meaning the wormhole closes and for $T > 1$ it pinches off with $r_"min" in CC$. This happens on a timescale
$
  dd(t, d: delta) tilde R/c
$
which is smaller than the time it takes for light to travel from I $->$ III.
