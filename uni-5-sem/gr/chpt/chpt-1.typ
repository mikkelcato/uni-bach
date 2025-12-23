//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Special Relativity
== The Lorentz transformations
Consider a system of pointlike particles interacting through gravity. We can write the force on the $N$th particle as
$
  m_N dv(bold(r)_N, t, 2) = sum_M (G m_N m_M )/abs(bold(r)_M-bold(r)_N)^2 (bold(r)_M - bold(r)_N)/abs(bold(r)_M - bold(r)_N)
$
which is just Newtons law of gravitation. This equation is easily seen to be invariant under Galilean transformations
$
  bold(r)' = R bold(r) + bold(v) t + bold(d)";  " t' = t + t_0
$
where $R$ is a member of $"O"(3)$. These form a $10$ parameter group called the Galilean group. Newtonian physics is only valid in frames related by the Galilean group. We call these _inertial frames_. But why are they special? Newton proposed there exists an _absolute frame_ with all inertial frames having some velocity $v_0$ (including $v_0 = 0$) with respect to the absolute frame.

The above becomes a larger problem when we consider electrodynamics. Maxwell's equations are famously not invariant under the Galilean group! An easy way to see this is the constancy of the speed of light $c$. However, they are invariant under Lorentz transformations
$
  X'^alpha = tensor(Lambda, alpha, -beta) X^beta + a^alpha
$
with $tensor(Lambda, alpha, -beta)$ being defined by the condition
$
  tensor(Lambda, alpha, -gamma) tensor(Lambda, beta, -delta) eta_(alpha beta) = eta_(gamma delta)
$
where $eta_(alpha beta)$ is the Minkowski metric. We take this as a symmetry of nature.

Since $a^alpha$ is constant we have
$
  dd(x'^alpha) = tensor(Lambda, alpha, -gamma) dd(x^gamma)
$
We define the proper time $tau$ by
$
  dd(tau^2) equiv - eta_(alpha beta) dd(x^alpha) dd(x^beta) = dd(t^2) - dd(bold(x)^2)
$
since $dd(s^2) = - dd(tau^2)$. We use _natural_ units where $c = 1$. This is invariant
$
  dd(tau'^2) &= - eta_(alpha beta) dd(x'^alpha) dd(x'^beta) \
  &= - eta_(alpha beta) tensor(Lambda, +alpha, -gamma) tensor(Lambda, +beta, -delta) dd(x^gamma) dd(x^delta) \
  &= - eta_(gamma delta) dd(x^gamma) dd(x^delta) = dd(tau^2)
$
Consider
$
  abs(dd(bold(x))/dd(t)) =^"for light" 1
$
implying $dd(tau) = 0$ for light. We just showed $dd(tau') = dd(tau)$ so this speed is a constant!



The Lorentz transformations form the Lorentz group. We are interested in the _proper_ Lorentz group satisfying
$
  tensor(Lambda, +0, -0) >= 1";  " det Lambda =+1
$
this excludes non-physical transformations. Taking $a^alpha =^! 0$ gives the _homogeneous_ proper Lorentz group satisfying
$
  tensor(Lambda, i, -j) = R_(i j)";  " tensor(Lambda, i, -0) = tensor(Lambda, 0, -i) = 0";  " tensor(Lambda, 0, -0) = 0
$
This is very similar to the Galilean group except we include boosts. We will now show this by determining $tensor(Lambda, mu, -nu)$.

Assume an observer $O$ sees a particle at rest while another observer $O'$ sees it having velocity $bold(v)$. Then
$
  dd(x'^alpha) = tensor(Lambda, alpha, -beta) dd(x^beta)
$
but $dd(bold(x)) = 0$ so
$
  dd(x'^i) = tensor(Lambda, i, -0) dd(t)";  " dd(t') = tensor(Lambda, 0, -0) dd(t)
$
and we obtain
$
  tensor(Lambda, i, -0) = v^i tensor(Lambda, 0, -0)
$
where $v^i = dd(x'^i)\/dd(t')$. Then using $eta_(gamma delta) = tensor(Lambda, alpha, -gamma) tensor(Lambda, beta, -delta) eta_(alpha beta)$ we find
$
  -1 & = tensor(Lambda, alpha, -0) tensor(Lambda, beta, -0) eta_(alpha beta) \
  & = sum_i (tensor(Lambda, i, -0))^2 - (tensor(Lambda, 0, -0))^2 = [sum_i (v^i)^2 - 1 ] (tensor(Lambda, 0, -0))^2 \
  tensor(Lambda, 0, -0) &= [1- underbracket(sum_i (v^i)^2, bold(v)^2)]^(-1\/2) equiv gamma
$
so $tensor(Lambda, i, -0) = gamma v^i$.

The rest of these notes use the Einstein summation convention where repeated indices are summed
$
  sum_i (tensor(Lambda, i, -0))^2 &equiv (tensor(Lambda, i, -0))^2 \
  sum_i tensor(Lambda, i, -0) tensor(Lambda, i, -j) &equiv tensor(Lambda, i, -0) tensor(Lambda, i, -j)
$

Consider

$
  0 &= tensor(Lambda, i, -0) tensor(Lambda, i, -j) - tensor(Lambda, 0, -0) tensor(Lambda, 0, -j) \
  &= gamma v^i tensor(Lambda, i, -j) - gamma tensor(Lambda, 0, -j) \
  tensor(Lambda, 0, -j) &= v^i tensor(Lambda, i, -j)
$
so we need to determine $tensor(Lambda, i, -j)$. We impose rotational symmetry about $hat(bold(v))$ and $tensor(Lambda, i, -j) -> bb(1)$ when $bold(v) = 0$. Under these $tensor(Lambda, i, -j)$ takes the form
$
  tensor(Lambda, i, -j) = tensor(delta, i, -j) + A v^i v_j
$
with $A$ being a function of $bold(v)^2$. Consider
$
  delta_(i j) &= tensor(Lambda, k, -i) tensor(Lambda, k, -j) - tensor(Lambda, 0, -i) tensor(Lambda, 0, -j) \
  &= tensor(Lambda, k, -i) tensor(Lambda, k, -j) - v^m tensor(Lambda, m, -i) v^n tensor(Lambda, n, -j) \
$
We compute
$
  tensor(Lambda, k, -i) tensor(Lambda, k, -j) &= delta_(i j) + (2 A + A^2 bold(v)^2) v_i v_j \
  v^k tensor(Lambda, k, -i) &= v_i (1 + A bold(v)^2)
$
substituting these give
$
  0 & = A^2 + 2/bold(v)^2 A - gamma^2/(bold(v)^2)
$
this is a quadratic in $A$ giving
$
  A = (gamma-1)/bold(v)^2
$
We find
$
  tensor(Lambda, i, -j) = tensor(delta, i, -j) + (gamma-1)/bold(v)^2 v^i v_j";  " tensor(Lambda, 0, -j) & = gamma v_j
$

== Tensors in special relativity
By definition we have
$
  dd(x'^alpha) = tensor(Lambda, alpha, -beta) dd(x^beta)
$
we call anything that transforms like $dd(x^beta)$ a contravariant four-vector. Then
$
  dd(x^beta) = (tensor(Lambda, alpha, -beta))^(-1) dd(x'^alpha)
$
implying
$
  pdv(x^beta, x'^alpha) = (tensor(Lambda, alpha, -beta))^(-1)
$
Then we have
$
  partial_alpha' equiv pdv(, x'^alpha) = (tensor(Lambda, alpha, -beta))^(-1) pdv(, x^beta)
$
we call anything that transforms like $partial_a'$ a covariant four-vector. These therefore transform inversely by definition. We can write $ (tensor(Lambda, alpha, -beta))^(-1) = tensor(Lambda, -alpha, beta) = eta_(alpha delta) eta^(gamma beta) tensor(Lambda, delta, -gamma) $
since
$
  tensor(Lambda, -alpha, gamma) tensor(Lambda, alpha, -beta) = eta_(alpha delta) eta^(gamma epsilon) tensor(Lambda, delta, -epsilon) tensor(Lambda, alpha, -beta) = eta_(epsilon beta) eta^(gamma epsilon) = tensor(delta, gamma, -beta)
$
where we use $eta^(beta delta) eta_(alpha delta) = tensor(delta, beta, -alpha)$. This makes contractions invariant
$
  U'_alpha V'^alpha = tensor(Lambda, -alpha, gamma) tensor(Lambda, alpha, -beta) U_gamma V_beta = tensor(delta, beta, -gamma) U_gamma V^beta = U_gamma V^gamma
$
which is very useful! All contravariant $V^beta$ have a covariant friend given by
$
  V_alpha equiv eta_(alpha beta) V^beta
$
and similarly all covariant $U_beta$ have a contravariant friend $U^alpha equiv eta^(alpha beta) U_beta$. Then we can raise and lower indices with $eta^(alpha beta)$ and $eta_(alpha beta)$
$
  eta^(alpha beta) V_beta = eta^(alpha beta) eta_(beta gamma) V^gamma = V^alpha
$
To see $V_alpha$ is covariant consider
$
  V'_alpha = eta_(alpha beta) V'^beta = eta_(alpha beta) tensor(Lambda, beta, -gamma)V^gamma = eta_(alpha beta) eta^(gamma delta) tensor(Lambda, beta, -gamma) V_delta = tensor(Lambda, -alpha, delta)V_delta
$
similarly one can show $U^alpha$ is contravariant.

Above we see four-vectors have one index. Tensors are more general objects with multiple indices. These transform in the obvious way
$
  tensor(T, gamma, -alpha beta) arrow tensor(T', gamma, -alpha beta) = tensor(Lambda, gamma, -delta) tensor(Lambda, -alpha, epsilon) tensor(Lambda, -beta, rho) tensor(T, delta, -epsilon rho)
$
We can contract tensors by $tensor(T, alpha gamma) equiv tensor(T, alpha, -beta, gamma beta)$

$
  tensor(T', alpha gamma) = tensor(T', alpha, -beta, gamma beta) &= tensor(Lambda, alpha, -delta) tensor(Lambda, -beta, epsilon) tensor(Lambda, gamma, -rho) tensor(Lambda, beta, -kappa) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(delta, -kappa, epsilon) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, -epsilon, rho epsilon) = tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, rho)
$
We can also take the direct product of two tensors giving a new tensor
$
  tensor(T, alpha, -beta, gamma) equiv tensor(A, alpha, -beta) tensor(B, gamma)
$
this is how we take derivatives
$
  tensor(T, -alpha, beta gamma) equiv partial_alpha tensor(T, beta gamma)
$
The linear combination of tensors is also a tensor
$
  tensor(T, alpha, -beta) equiv a tensor(R, alpha, -beta) + b tensor(S, alpha, -beta)
$
These imply that $tensor(delta, alpha, -beta)$ is a tensor since $eta_(alpha beta)$ is a tensor. Similarly raising or lowering an index preserves _tensorness_

= The equivalence principle
== Newtonian field theory
Newtonian gravity is summarized by
$
  bold(F) = - G_N (m M)/r^2 hat(bold(r)) = m bold(g)
$
where we define the gravitational field
$
  bold(g) = - G_N M/r^2 hat(bold(r))
$
We can integrate this to find
$
  underbracket(integral.cont_S bold(g) dot dd(bold(A)), "gravitational flux" #linebreak() "through" S) = - 4 pi G_N M_"enclosed"
$
using Gauss' law we can rewrite this as
$
  integral div bold(g) dd(V) = - 4 pi G_N integral rho dd(V)
$
this is true for any volume $V$ meaning
$
  div bold(g) = - 4 pi G_N rho
$
We define $bold(g) equiv - grad Phi$ with $Phi$ being the gravitational potential giving
$
  laplacian Phi = 4 pi G_N rho
$
this is the Newtonian field equation. This equation describes how some matter distribution $rho$ shapes $Phi$. Newton's second law can also be written in terms of $Phi$ giving
$
  dv(bold(r), t, 2) = - grad Phi
$
this equation describes how a particle moves given $Phi$.

These two equations make up the Newtonian field theory of gravity. We immediately see these are incompatible with special relativity since time and space are not treated equally. Actually it is a static theory due to the _action at a distance_ description which underlies Newtonian gravity.

== The principle
We have two kinds of mass in Newtonian physics. The inertial mass $m_i$ and the gravitational mass $m_g$. These are defined by
$
  bold(F) = m_i bold(a)";  " bold(F) = m_g bold(g)
$
The basic statement of the equivalence principle is $m_i =^! m_g$ meaning we can identify $ dot.double(bold(x)) = bold(g) $

Assuming $dot.double(bold(x)) = bold(g)$ we can always do a coordinate transformation to a local frame with no acceleration $dot.double(bold(y)) = 0$ by
$
  bold(y) = bold(x) - 1/2 bold(g) t^2
$
as we will see this is very important! We call these locally inertial frames for freely falling elevators.

Einstein generalized this to all of physics by the strong equivalence principle.

_In any gravitational field it is possible to select a locally inertial system such that the laws of physics are the same as in special relativity._

So we can always do a coordinate transformation to a freely falling elevator with no gravity!

= The geodesic equation
== The metric
Above we had
$
  dd(tau)^2 = - eta_(alpha beta) dd(y)^alpha dd(y)^beta
$
under the assumption of no gravity. With gravity this is only true locally $y^alpha -> y^alpha (x)$
$
  dd(tau)^2 &= - eta_(alpha beta) dd(y)^alpha (x) dd(y)^beta (x) \
  &=^"chain rule" - eta_(alpha beta) pdv(y^alpha (x), x^mu) pdv(y^beta (x), x^nu) dd(x)^mu dd(x)^nu \
  &= - g_(mu nu) (x) dd(x)^mu dd(x)^nu
$
this is a global relation. We have defined the metric $g$
$
  g_(mu nu) (x) = eta_(alpha beta) pdv(y^alpha (x), x^mu) pdv(y^beta (x), x^nu)
$
which is now space-dependent meaning space becomes non-Euclidean. Therefore $g$ contains the effect gravity!

== The equation
Locally we have
$
  dv(y^alpha (x), tau, 2) = 0
$
We can rewrite this as
$
  0 & = dv(, tau) (dv(y^alpha (x), tau)) \ 0 &= dv(, tau) (pdv(y^alpha, x^mu) pdv(x^mu (tau), tau)) \
  0 & = pdv(y^alpha, x^mu) dv(x^mu, tau, 2) + pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau) \
  0 & = pdv(x^lambda, y^alpha) (pdv(y^alpha, x^mu) dv(x^mu, tau, 2) + pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau)) \
  0 & = tensor(delta, lambda, -mu) dv(x^mu, tau, 2) + pdv(x^lambda, y^alpha) pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau) \
  0 & = dv(x^lambda, tau, 2) + tensor(Gamma, lambda, -mu nu) dv(x^mu, tau) dv(x^nu, tau) equiv dot.double(x)^lambda + tensor(Gamma, lambda, -mu nu) dot(x)^mu dot(x)^nu
$
This is called the _geodesic equation_. We have defined
$
  tensor(Gamma, lambda, -mu nu) equiv pdv(x^lambda, y^alpha) pdv(y^alpha, x^mu, x^nu) =^"flat space" 0
$
which is the _Christoffel symbol_. With our definition it is _torsion-free_ meaning symmetric in the lower indices $tensor(Gamma, lambda, -mu nu) = tensor(Gamma, lambda, -nu mu)$.

Since the Christoffel symbol vanishes in flat space there must be some dependence on $g$. We have
$
  tensor(Gamma, lambda, -mu nu) pdv(y^beta, x^lambda) & = pdv(y^beta, x^mu, x^nu)
$
Then
$
  pdv(g_(mu nu), x^lambda) &= pdv(, x^lambda) (eta_(alpha beta) pdv(y^alpha, x^mu) pdv(y^beta, x^nu)) \
  &= eta_(alpha beta) pdv(y^alpha, x^lambda, x^mu) pdv(y^beta, x^nu) + eta_(alpha beta) pdv(y^beta, x^lambda, x^nu) pdv(y^alpha, x^mu) \
  &=^"by above" eta_(alpha beta) pdv(y^beta, x^nu) tensor(Gamma, rho, -lambda mu) pdv(y^alpha, x^rho) + eta_(alpha beta) pdv(y^alpha, x^mu) tensor(Gamma, rho, -lambda nu) pdv(y^beta, x^rho) \
  &= g_(rho nu) tensor(Gamma, rho, -mu lambda) + g_(rho mu) tensor(Gamma, rho, -nu lambda)
$
Using this and symmetry we have
$
  pdv(g_(mu nu), x^lambda) + pdv(g_(lambda nu), x^mu) - pdv(g_(mu lambda), x^nu) = 2 g_(sigma nu) tensor(Gamma, rho, -lambda mu)
$
To isolate $Gamma$ we need the inverse metric $g^(mu nu)$ satisfying
$
  g_(mu nu) (x) g^(nu sigma) (x) = tensor(delta, sigma, -mu)
$
We claim
$
  g^(nu sigma) = eta^(alpha beta) pdv(x^nu, y^alpha) pdv(x^sigma, y^beta)
$
This is simple to show
$
  g_(mu nu) g^(nu sigma) &= eta_(gamma delta) eta^(alpha beta) pdv(y^gamma, x^mu) underbracket(pdv(y^delta, x^nu) pdv(x^nu, y^alpha), tensor(delta, delta, -alpha)) pdv(x^sigma, y^beta) = underbracket(eta_(gamma delta) eta^(delta beta), tensor(delta, beta, -gamma)) pdv(y^gamma, x^mu) pdv(x^sigma, y^beta) = tensor(delta, sigma, -mu)
$

Then we find
$
  tensor(Gamma, lambda, -mu nu) = 1/2 g^(lambda sigma) [pdv(g_(nu sigma), x^mu) + pdv(g_(mu sigma), x^nu) - pdv(g_(mu nu), x^sigma)] = 1/2 g^(lambda sigma) [partial_mu g_(nu sigma) + partial_nu g_(mu sigma) underbracket(- partial_sigma g_(mu nu), "symmetric")]
$
where the dependence on $g$ is manifest.

== The Newtonian limit
For the above to be consistent we should be able to recover Newtonian gravity in what we call _the Newtonian limit_. In this limit all velocities are small
$
  abs(dv(bold(x), tau)) << 1
$
and the metric $g$ is static. We also assume gravity is weak meaning we can write $ g_(mu nu) (x) = eta_(mu nu) + h_(mu nu) (x) $
with $abs(h) << 1$.

Assuming small velocities the geodesic equation reduces to
$
  0 & = dot.double(x)^mu + tensor(Gamma, mu, -nu lambda) dot(x)^nu dot(x)^lambda \
    & tilde.eq dot.double(x)^mu + tensor(Gamma, mu, -00) (dv(t, tau))^2
$
Assuming a static $g$ the Christoffel symbol becomes
$
  tensor(Gamma, mu, -00) &= 1/2 g^(mu sigma) [partial_0 g_(0 sigma) + partial_0 g_(0 sigma) - partial_sigma g_(00)] \
  &tilde.eq^("static" g) - 1/2 g^(mu sigma) partial_sigma g_00
$

Assuming gravity is weak we have to order $cal(O)(h)$
$
  tensor(Gamma, mu, -00) tilde.eq - 1/2 eta^(mu sigma) partial_sigma h_00
$
Since $h_00$ does not depend on time we find $tensor(Gamma, 0, -00) = 0$ meaning
$
  dv(t, tau) = "constant"
$
by the geodesic equation. Then for $mu = i$
$
  0 & tilde.eq dv(bold(x), tau, 2) - 1/2 (dv(t, tau))^2 grad h_00 (bold(x)) \
    & eq dv(bold(x), t, 2) - 1/2 grad h_00 (bold(x))
$
By comparison we identify
$
  h_00 = - 2 Phi + underbracket("constant", "can absorb")
$
meaning
$
  g_00 = - (1 + 2 Phi)
$
which is quite nice!

This already leads to non-trivial results. Consider two clocks at rest. Both satisfy
$
  dd(tau)^2 = - g_(00) dd(t)^2
$
Then since $omega prop dd(tau)^(-1)$ we have to order $cal(O)(Phi)$
$
  (omega_2-omega_1)/omega_1 & = sqrt((-g_00 (x_1))/(-g_00 (x_2)))- 1 \
                            & eq sqrt((1+2 Phi(x_1))/(1+2 Phi(x_2))) - 1 \
                            & tilde.eq^"Taylor" Phi(x_1) - Phi(x_2)
$
so
$
  dd(omega, d: Delta)/omega_1 = - dd(Phi, d: Delta)
$
This is gravitational redshift!

= The principle of general covariance
== The principle
We want to express the laws of physics such that they hold in all frames. This is called _general covariance_. Specifically we want to write the laws of special relativity in covariant form. Then by definition they would be valid in all frames including those with gravity. This is the principle of general covariance.

Before we see how this is done consider the geodesic equation
$
  dot.double(x)^lambda + underbracket(tensor(Gamma, lambda, -mu nu) dot(x)^mu dot(x)^nu, "gravity") = 0
$
At all points $x = tilde(x)$ we can by the equivalence principle define a freely falling elevator where
$
  g_(mu nu) (tilde(x)) = eta_(mu nu)
$
This requires the Christoffel symbol vanishes at $x = tilde(x)$ meaning
$
  evaluated(partial_sigma g_(mu nu) (x))_(x=tilde(x)) = 0
$
But we do not require
$
  evaluated(partial_sigma partial_rho g_(mu nu) (x))_(x=tilde(x)) = 0
$
so if we move away from $x=tilde(x)$ the Christoffel symbol is not guaranteed to vanish. And in fact if this were true then $g_(mu nu) = eta_(mu nu)$ everywhere! When we quantify curvature later this becomes important.


== Tensors in general relativity
We saw before how tensors in special relativity are objects that transform in a specific way under Lorentz transformations. This notion generalizes to general tensors being objects that transform in a specific way under general coordinate transformations.

Scalars are objects like $dd(tau)$ and $phi(x)$. These are invariant under $x -> x'$ meaning they transform trivially.

$dd(x^mu)$ transforms as
$
  dd(x'^mu) =^"chain rule" pdv(x'^mu, x^nu) dd(x^nu)
$
anything transforming like $dd(x^mu)$ is a contravariant vector. Anything transforming inversely of this is a covariant vector
$
  A'_mu = pdv(x^nu, x'^mu) A_nu
$
This is the inverse since contracting gives a scalar
$
  A'_mu U'^mu = underbracket(pdv(x^nu, x'^mu) pdv(x'^mu, x^sigma), tensor(delta, nu, -sigma)) A_nu U^sigma = A_nu U^nu
$
We can also form a covariant vector by differentiating a scalar
$
  partial_mu' phi' (x') eq^"chain rule" pdv(x^nu, x'^mu) underbracket(partial_nu phi(x), A_mu)
$
But we will find that $partial_mu U^nu$ is not generally a tensor.

General tensors transform in the obvious way. As an example take the metric $g_(mu nu)$
$
  g'_(mu nu) (x') = eta_(alpha beta) pdv(y^alpha, x'^mu) pdv(y^beta, x'^nu) = eta_(alpha beta) pdv(y^alpha, x^rho) pdv(y^beta, x^sigma) pdv(x^rho, x'^mu) pdv(x^sigma, x'^nu) = pdv(x^rho, x'^mu) pdv(x^sigma, x'^nu) g_(rho sigma) (x)
$
so it is covariant. Similarly the inverse metric $g^(mu nu)$ is contravariant. Using the metric we define
$
  tilde(T)_(mu nu) equiv overbracket(underbracket(g_(mu sigma) g_(nu rho), 2 times "covariant") underbracket(T^(sigma rho), "contravariant"), "covariant") equiv T_(mu nu)
$
By comparison $tilde(T)_(mu nu)$ is covariant since the indices $sigma rho$ _cancel_ leaving the covariant $mu nu$. We write $tilde(T)_(mu nu) equiv T_(mu nu)$ with $T_(mu nu)$ and $T^(mu nu)$ representing the same object. Similarly $tensor(delta, mu, -nu)$ is a mixed tensor $tensor(delta, mu, -nu) = g^(mu sigma) g_(sigma nu)$ since the $sigma$ _cancels_.

We define $g equiv det g_(mu nu)$. This is used to define the invariant measure
$
  underbracket(sqrt(-g), "inverse Jacobian") dd(x, 4)
$

== The covariant derivative
As mentioned $partial_mu A^nu$ is not a tensor. This is a problem since we would like to take derivatives. Consider
$
  "invariant" & = dv(phi(x), tau) \
              & = pdv(phi, x^mu) (dv(x^mu, tau))
$
This implies
$
  "invariant" &= dv(phi, tau, 2) \ &= pdv(phi, x^mu, x^nu) dv(x^nu, tau) dv(x^mu, tau) + pdv(phi, x^mu) dv(x^mu, tau, 2)
$
Assuming $x^mu (tau)$ is described by the geodesic equation we find
$
  underbracket(dv(phi, tau, 2), "scalar") = overbracket([pdv(phi, x^mu, x^nu)-tensor(Gamma, sigma, -mu nu) pdv(phi, x^sigma)], "must be covariant") times underbracket(dv(x^mu, tau) dv(x^nu, tau), "contravariant")
$
We call the object in $[dots]$ the covariant derivative. Let $V_mu = partial_mu phi$ then $[dots]$ becomes
$
  D_mu V_nu = partial_mu V_nu - tensor(Gamma, sigma, -mu nu) V_sigma
$
and similarly
$
  D_mu V^nu = partial_mu V^nu + tensor(Gamma, nu, -mu sigma) V^sigma
$
which generalizes in the obvious way.

We see $D_mu -> partial_mu$ in the freely falling elevator since $Gamma$ vanishes. This also implies that both $tensor(Gamma, lambda, -mu nu)$ and $partial_mu$ are not tensors.

Consider $D_sigma g_(mu nu)$. For a freely falling elevator we have
$
  D_sigma g_(mu nu) & =^"free fall"_("possible by" #linebreak() "equivalence principle") partial_sigma g_(mu nu) =^"flat locally" 0
$
Since $D_sigma g_(mu nu)$ is written in covariant form then the above is true in all frames when it is true in one frame! This relies on us being able to pick a freely falling elevator, this is possible by the equivalence principle. Generally any equation that is true in flat space becomes true in all frames by the replacement $partial_mu -> D_mu$.

Consider
$
  dv(V^nu, tau, d: D) &equiv dv(x^mu, tau) D_mu V^nu \ &= dv(V^nu, tau) + tensor(Gamma, nu, -mu sigma) V^sigma dv(x^mu, tau)
$
We have $D_tau -> d_tau$ in the freely falling elevator. We can now immediately derive the geodesic equation by
$
  dv(V^mu, tau, d: D) = 0
$
which reduces to
$
  dv(x^mu, tau, 2) = 0
$
in the freely falling elevator.


== Electrodynamics
As an example of the procedure mentioned above we consider classical electrodynamics.



We define the antisymmetric Faraday tensor
$ F_(mu nu) equiv partial_mu A_nu - partial_nu A_mu $
and the four-current $J^mu equiv vecrow(rho, bold(J))$. Then Maxwell's equations can be written as
$
  partial_mu F^(mu nu) &= - J^nu \
  0 &= partial_mu F_(nu gamma) + partial_(gamma) F_(mu nu) + partial_nu F_(gamma mu)
$
To include gravity we replace $partial_mu -> D_mu$
$
  D_mu F^(mu nu) & = - J^nu \
  0 & = D_mu F_(nu lambda) + D_lambda F_(mu nu) + D_nu F_(lambda mu) =^(Gamma "cancel") partial_mu F_(nu gamma) + partial_gamma F_(mu nu) + partial_nu F_(gamma mu)
$
We can show
$
  D_mu F^(mu nu) = 1/sqrt(-g) partial_mu (sqrt(-g) F^(mu nu))
$
which holds for any antisymmetric $F^(mu nu)$. Then by $D_mu F^(mu nu) = - J^nu$ we have
$
  partial_mu (sqrt(-g) F^(mu nu)) = - sqrt(-g) J^nu
$
Then
$
  underbracket(partial_nu partial_mu, "symmetric") overbracket((sqrt(-g) F^(mu nu)), "antisymmetric") = 0
$
implying
$
  partial_nu (sqrt(-g) J^nu) = 0
$
This is the covariant form of $partial_mu J^mu = 0$.

Consider the Lorentz force
$
  f^mu = e tensor(F, mu, -nu) U^nu
$
this is already covariant!
