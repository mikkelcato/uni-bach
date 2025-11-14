//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Introduction
Relativity is in short the idea that only relative motion is measurable, this can be stated as a symmetry---namely that our equations should be unchanged under coordinate transformations. Special relativity is the symmetry with respect to inertial frames, and general relativity extends this to general frames.

#pagebreak()
= Special Relativity
== Galileo and Lorentz
If we consider a system of point particles interacting gravitationally we can write
$
  m_N dv(arrow(r)_N, t, 2) = G sum_M (m_N m_M (arrow(r)_M-arrow(r)_N))/abs(arrow(r)_M-arrow(r)_N)^3
$
as the force on the particle $N$. This equation is invariant under the Galilean transformations
$
  arrow(r)' = R arrow(r) + arrow(v) t + arrow(d)",  " t' = t + t_0
$
where $R$ is a member of $O(3)$---a rotation group---the other symbols have the obvious meanings. The Galilean transformations evidently form a $10$-parameter group, called the Galilean group. Newton was confused as to why there are many other transformations which do not leave Newtonian physics invariant. Newton's laws only hold in frames related by the Galilean transformations---the inertial frames. Newton wanted to know what made these frames special---he believed there must be some absolute space and that inertial frames are at rest or move with uniform velocity within it.

The problems start for real with Maxwell's electromagnetism, which is not invariant under Galilean transformations---e.g. because the speed of light is a universal constant. Initially it was believed that light propagated in the ether with electromagnetism only being valid in the frame at rest with respect tothe ether---this was disproven by Michelson and Morley. Einstein solved this problem by replacing Galilean invariance with Lorentz invariance, which correspond to the Lorentz transformations
$
  X'^alpha = tensor(Lambda, +alpha, -beta) X^beta + a^alpha
$
$tensor(Lambda, +alpha, -beta)$ is defined such that $tensor(Lambda, +alpha, -gamma) tensor(Lambda, +beta, -delta) tensor(eta, -alpha, -beta) = tensor(eta, -gamma, -delta)$ with $eta_(alpha beta)$ being the usual Minkowski metric defined as $eta_(alpha beta) = g_(alpha beta) equiv e_alpha dot e_beta$---proceeding in this manner we derive everything from this symmetry, which is essentially the cleanest way to do it. Since $a^alpha$ is constant we find
$
  dd(x'^alpha) = tensor(Lambda, +alpha, -gamma) dd(x^gamma)
$
we define the proper time (with $c equiv 1$)
$
  dd(tau^2) equiv dd(t^2)-dd(arrow(x)^2) = - eta_(alpha beta) dd(x^alpha) dd(x^beta)
$
this is from $dd(s^2) = -c^2 dd(t^2)+dd(arrow(x)^2)$ with $dd(arrow(x))=0 => dd(t)=dd(tau)$ giving $dd(s^2) = -c^2 dd(tau^2)$. This is invariant
$
  dd(tau'^2) &= - eta_(alpha beta) dd(x'^alpha) dd(x'^beta) \
  &= - eta_(alpha beta) tensor(Lambda, +alpha, -gamma) tensor(Lambda, +beta, -delta) dd(x^gamma) dd(x^delta) \
  &= - eta_(gamma delta) dd(x^gamma) dd(x^delta) \
  &= dd(tau^2)
$
This gives us a universal speed of light, since we have for light
$
  abs(dv(arrow(x), t)) = c = 1 => dd(tau) = 0
$
so light travels along null-geodesics---as will be explained later. Since $dd(tau)$ is invariant we have $dd(tau')=0$ so $ abs(dv(arrow(x)', t')) = 1 $ i.e. the speed of light is the same in all coordinates.

These transformations form the Lorentz group, typically we are interested in the proper Lorentz groups, which is a subgroup satisfying
$
  tensor(Lambda, +0, -0) >= 1",  " det Lambda =+1
$
this excludes non-fundamental symmetries corresponding to non-physical transformations. Further taking $a^alpha = 0$ gives us the homogeneous proper Lorentz groups, this has a subgroup with
$
  tensor(Lambda, +i, -j) = R_(i j)",  " tensor(Lambda, +i, -0)=tensor(Lambda, +0, -i) = 0",  " tensor(Lambda, 0, -0) = 1
$
so the difference between this and the Galilean group are the boosts. Let's assume an observer $cal(O)$ sees a particle at rest, while another observer $cal(O)'$ sees it having velocity $arrow(v)$, then
$
  dd(x'^alpha)=tensor(Lambda, alpha, -beta) dd(x^beta)
$
but we know $dd(arrow(x))=0$ so
$
  dd(x'^i) = tensor(Lambda, i, -0) dd(t)",  " dd(t') = tensor(Lambda, 0, -0) dd(t)
$
combining these give
$
  v^i = dv(x'^i, t') => tensor(Lambda, i, -0) = v^i tensor(Lambda, 0, -0)
$
and from $eta_(gamma delta) = tensor(Lambda, alpha, -gamma) tensor(Lambda, beta, -delta) eta_(alpha beta)$ we get
$
  -1 = tensor(Lambda, alpha, -0) tensor(Lambda, beta, -0) eta_(alpha beta) = sum_i (tensor(Lambda, i, -0))^2 - (tensor(Lambda, 0, -0))^2
$
which upon plugging in the previous gives
$
  tensor(Lambda, 0, -0) = gamma",  " tensor(Lambda, i, -0) = gamma v^i
$
the other elements are not uniquely determined, but we can pick
$
  tensor(Lambda, i, -j) = delta_(i j) + v_i v_j (gamma - 1)/v^2",  " tensor(Lambda, 0, -j) = gamma v_j
$

== Force, energy and momentum
We can generalize Newton's second law with
$
  f^alpha = m dv(x^alpha, tau, 2)
$
if a particle is a rest then $dd(tau)=dd(t)$ and this reduces to the usual non-relativistic force, and since $dd(x'^alpha)=tensor(Lambda, alpha, -beta)dd(x^beta)$ with $dd(tau)=dd(tau')$ this implies that $f'^alpha = tensor(Lambda, alpha, -beta) f^beta$, so it is a four-vector---this is nice because we can construct the four-force from the non-relativistic force.

Similarly we can define four-momentum
$
  p^alpha = m dv(x^alpha, tau)
$
we can then write Newton's law with a net force as
$
  dv(p^alpha, tau) = f^alpha
$
this is nice because it is covariant. Using $dd(tau) = sqrt(dd(t^2)-dd(arrow(x)^2)) = sqrt(1-v^2) dd(t) = dd(t)\/gamma$ and $arrow(v) = dd(arrow(x))\/dd(t)$ we can find
$
  p^0 = m gamma equiv E",  " arrow(p)=m gamma arrow(v)
$
in the low-velocity limit these become
$
  arrow(p) = m arrow(v) + cal(O)(v^3)",  " E = m + 1/2 m v^2 + cal(O)(v^4)
$
By linearity of the Lorentz transformation, if $p^alpha$ is conserved then $p'^alpha$ is conserved.

== Four-vectors and Tensors
Any four-vector transforming as
$
  V^alpha arrow V'^alpha = tensor(Lambda, alpha, -beta) V^beta", for" x^alpha arrow x'^alpha = tensor(Lambda, alpha, -beta)x^beta
$
is contravariant, and any four-vector transforming inversely
$
  U_alpha arrow U'_alpha = tensor(Lambda, -alpha, beta) U_beta
$
with $tensor(Lambda, -alpha, beta) equiv eta_(alpha gamma) eta^(beta delta) tensor(Lambda, gamma, -delta)$ and $eta^(beta delta)=eta_(beta delta)$ with $eta^(beta delta) eta_(alpha delta) = tensor(delta, beta, -alpha)$ is covariant. These are inverses, meaning
$
  tensor(Lambda, -alpha, gamma) tensor(Lambda, alpha, -beta) = eta_(alpha delta) eta^(gamma epsilon) tensor(Lambda, delta, -epsilon) tensor(Lambda, alpha, -beta) = eta_(epsilon beta) eta^(gamma epsilon) = tensor(delta, -beta, gamma)
$
this means that contractions become invariant
$
  U'_alpha V'^alpha = tensor(Lambda, -alpha, gamma) tensor(Lambda, +alpha, -beta) U_gamma V^beta = tensor(delta, -beta, gamma) U_gamma V^beta = U_beta V^beta
$
which is very nice. Every contravariant four-vector has a covariant buddy and vice versa given by
$
  V_alpha equiv eta_(alpha beta) V^beta",  " U^alpha equiv eta^(alpha beta) U_beta
$
so we can raise and lower indices,
$
  eta^(alpha beta) V_beta = eta^(alpha beta) eta_(beta gamma) V^gamma = V^alpha
$
we can check that $V_alpha$ is actually covariant,
$
  V'_alpha = eta_(alpha beta) V'^beta = eta_(alpha beta) tensor(Lambda, beta, -gamma)V^gamma = eta_(alpha beta) eta^(gamma delta) tensor(Lambda, beta, -gamma) V_delta = tensor(Lambda, -alpha, delta)V_delta
$
similarly we can check that $U^alpha$ is contravariant,
$
  U'^alpha = eta^(alpha beta) U'_beta = eta^(alpha beta) tensor(Lambda, -beta, gamma) U_gamma = eta^(alpha beta) eta_(gamma delta) tensor(Lambda, -beta, gamma) U^delta = tensor(Lambda, alpha, -delta) U^delta
$

Tensors are essentially generalized four-vectors, so they have more indices, where contra- and covariant indices transform as one would expect,
$
  tensor(T, gamma, -alpha beta) arrow tensor(T', gamma, -alpha beta) = tensor(Lambda, gamma, -delta) tensor(Lambda, -alpha, epsilon) tensor(Lambda, -beta, rho) tensor(T, delta, -epsilon rho)
$
indices within a tensor can also be contracted $tensor(T, alpha gamma) equiv tensor(T, alpha, -beta, gamma beta)$,
$
  tensor(T', alpha gamma) = tensor(T', alpha, -beta, gamma beta) &= tensor(Lambda, alpha, -delta) tensor(Lambda, -beta, epsilon) tensor(Lambda, gamma, -rho) tensor(Lambda, beta, -kappa) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(delta, -kappa, epsilon) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, -epsilon, rho epsilon) = tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, rho)
$
we can also make many other construction which are also tensors, e.g. by taking derivatives
$
  tensor(T, -alpha, beta gamma) equiv pdv(, x^alpha) tensor(T, beta gamma)
$
is a tensor when $tensor(T, beta gamma)$ is a tensor---this is more obvious if we write,
$
  partial_alpha equiv pdv(, x^alpha)
$
which is covariant. The direct product of two tensors is a tensor,
$
  tensor(T, alpha, -beta, gamma) equiv tensor(A, alpha, -beta) tensor(B, gamma)
$
the linear combination of tensors is a tensor,
$
  tensor(T, alpha, -beta) equiv a tensor(R, alpha, -beta) + b tensor(S, alpha, -beta)
$
the Minkowski metric is by definition a tensor, implying that $tensor(delta, alpha, -beta)$ is a tensor---this also implies that raising or lowering indices $tensor(T, -alpha, delta, -gamma) equiv eta^(delta beta) tensor(T, -alpha beta gamma)$ conserves tensor-ness. The Levi-Civita tensor $epsilon^(alpha beta gamma delta)=-epsilon_(alpha beta gamma delta)$ is a tensor. The zero tensor $tensor(T, alpha, -beta)=0$ is always zero, $tensor(T', alpha, -beta)=0$, since $tensor(T, alpha, -beta)=tensor(S, alpha, -beta) => tensor(T', alpha, -beta)=tensor(S', alpha, -beta)$.

#pagebreak()
= The equivalence principle
== Newtonian field theory
Newton's theory of gravity can be summarized by
$
  arrow(F) = - G_N (m M)/r^2 hat(r) = m arrow(g)
$
where we introduce the gravitational field
$
  arrow(g) = - G_N M/r^2 hat(r)
$
this can be expressed as
$
  integral.cont_S arrow(g) dot dd(arrow(A)) = - 4 pi G_N M
$
where the LHS represents the gravitational flux through any closed surface $S$ and $M$ is the total mass enclosed. As is familiar this can be rewritten as follows
$
  integral div arrow(g) dd(V) = - 4 pi G_N integral rho dd(V)
$
this must hold for any volume, and if we define $arrow(g) equiv - grad Phi$, with $Phi$ the gravitational potential, we get the field equation
$
  laplacian Phi = 4 pi G_N rho
$
from Newton's second law $arrow(F) = m arrow(a)$ we also get
$
  dv(arrow(r), t, 2) = - grad Phi
$
Importantly the Newtonian field theory of gravitation is clearly incompatible with special relativity as time and space are not treated equally---in fact it is a static theory, this reflects the underlying physics admitting an action at a distance description, implying an infinite signal speed.

#pagebreak()
== The principle
In Newtonian mechanics we have two kinds of mass the inertial mass $m_i$ and the gravitational mass $m_g$ defined by
$
  arrow(F) = m_i arrow(a)"  and  "arrow(F) = m_g arrow(g)
$
the most basic statement of the equivalence principle (EP) is that $m_i = m_g$, which lets us identify $dot.double(arrow(x)) = arrow(g)$---notably this is an empirical law, which is why we call it a principle.

The principle implies that we can always do a coordinate transformation to some local frame with no acceleration by
$
  arrow(y) = arrow(x) - 1/2 arrow(g) t^2 => dot.double(arrow(y)) = 0
$
if $dot.double(arrow(x)) = arrow(g)$---so locally $arrow(g) tilde arrow(a)$.

What Einstein did was generalize this to all of physics by the strong equivalence principle (strong EP)---which can be stated as:

_"In any gravitational field it is possible to select a locally inertial system (i.e. we can pick $arrow(y)$ with no gravity)---a freely falling elevator---such that the laws of physics are the same as in special relativity---i.e. no gravity."_

From this principle everything in general relativity essentially follows, as will be seen the next section.

== The metric and the geodesic equation
In special relativity we had
$
  dd(tau)^2 = - eta_(alpha beta) dd(y)^alpha dd(y)^beta
$
here space is flat and there is no gravity ($x_i = x^i$). With gravity the previous is still true locally $y^alpha -> y^alpha (x)$ so
$
  dd(tau)^2 = - eta_(alpha beta) dd(y)^alpha (x) dd(y)^beta (x)
$
to get a global relation we have by the chain rule
$
  dd(tau)^2 &= - eta_(alpha beta) pdv(y^alpha (x), x^mu) pdv(y^beta (x), x^nu) dd(x)^mu dd(x)^nu \
  &= - g_(mu nu) (x) dd(x)^mu dd(x)^nu
$
so now the metric $g_(mu nu)$ has become space dependent, with
$
  g_(mu nu) (x) = eta_(alpha beta) pdv(y^alpha (x), x^mu) pdv(y^beta (x), x^nu)
$
now the metric obviously becomes space-dependent, and hence space becomes non-euclidean---so $g_(mu nu)$ contains gravity.

Locally we have
$
  dv(y^alpha (x), tau, 2) = 0
$
which can be written as
$
  0 & = dv(, tau) (dv(y^alpha (x), tau)) \
  & = dv(, tau) (pdv(y^alpha, x^mu) pdv(x^mu (tau), tau)) \
  & = pdv(y^alpha, x^mu) dv(x^mu, tau, 2) + pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau) \
  & = pdv(x^lambda, y^alpha) (pdv(y^alpha, x^mu) dv(x^mu, tau, 2) + pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau)) \
  & = tensor(delta, lambda, -mu) dv(x^mu, tau, 2) + pdv(x^lambda, y^alpha) pdv(y^alpha, x^mu, x^nu) dv(x^mu, tau) dv(x^nu, tau) \
  & = dv(x^lambda, tau, 2) + tensor(Gamma, lambda, -mu nu) dv(x^mu, tau) dv(x^nu, tau)
$
where we define the Christoffel symbol (or affine connection) as
$
  tensor(Gamma, lambda, -mu nu) = pdv(x^lambda, y^alpha) pdv(y^alpha, x^mu, x^nu)
$
which is zero for flat space. This is our equation of motion in general relativity, and it is just the geodesic equation---as we see later.

Since the Christoffel symbol is zero for flat space it in some way encodes gravity, meaning there must be some relation between $Gamma tilde g$. To this end note that
$
  tensor(Gamma, lambda, -mu nu) pdv(y^beta, x^lambda) &= pdv(y^alpha, x^mu, x^nu) pdv(x^lambda, y^alpha) pdv(y^beta, x^lambda) \
  &= pdv(y^alpha, x^mu, x^nu) tensor(delta, beta, -alpha) \
  &= pdv(y^beta, x^mu, x^nu)
$
then
$
  pdv(g_(mu nu), x^lambda) &= pdv(, x^lambda) (eta_(alpha beta) pdv(y^alpha, x^mu) pdv(y^beta, x^nu)) \
  &= eta_(alpha beta) pdv(y^alpha, x^lambda, x^mu) pdv(y^beta, x^nu) + eta_(alpha beta) pdv(y^beta, x^lambda, x^nu) pdv(y^alpha, x^mu) \
  &= eta_(alpha beta) pdv(y^beta, x^nu) tensor(Gamma, rho, -lambda mu) pdv(y^alpha, x^rho) + eta_(alpha beta) pdv(y^alpha, x^mu) tensor(Gamma, rho, -lambda nu) pdv(y^beta, x^rho) \
  &= g_(rho nu) tensor(Gamma, rho, -mu lambda) + g_(rho mu) tensor(Gamma, rho, -nu lambda)
$
using the symmetry of the Christoffel symbol we can then find
$
  pdv(g_(mu nu), x^lambda) + pdv(g_(lambda nu), x^mu) - pdv(g_(mu lambda), x^nu) = 2 g_(sigma nu) tensor(Gamma, rho, -lambda mu)
$
now we define the inverse metric $g^(mu nu)$ by
$
  g_(mu nu) (x) g^(nu sigma) (x) = tensor(delta, sigma, -mu)
$
we claim that
$
  g^(nu sigma) = eta^(alpha beta) pdv(x^nu, y^alpha) pdv(x^sigma, y^beta)
$
is the inverse. This is easily checked
$
  g_(mu nu) g^(nu sigma) &= eta_(gamma delta) pdv(y^gamma, x^mu) pdv(y^delta, x^nu) eta^(alpha beta) pdv(x^nu, y^alpha) pdv(x^sigma, y^beta) \
  &= eta_(gamma delta) eta^(alpha beta) tensor(delta, delta, -alpha) pdv(y^gamma, x^mu) pdv(x^sigma, y^beta) \
  &= eta_(gamma delta) eta^(delta beta) pdv(y^gamma, x^mu) pdv(x^sigma, y^beta) \
  &= tensor(delta, beta, -gamma) pdv(y^gamma, x^mu) pdv(x^sigma, y^beta) \
  &= pdv(y^beta, x^mu) pdv(x^sigma, y^beta) \
  &= tensor(delta, sigma, -mu)
$
so the inverse exists and we can write the Christoffel symbol as
$
  tensor(Gamma, lambda, -mu nu) = 1/2 g^(lambda sigma) [pdv(g_(nu sigma), x^mu) + pdv(g_(mu sigma), x^nu) - pdv(g_(mu nu), x^sigma)]
$
and this clearly shows that the metric encodes gravity. Which is why as we'll see solving the Einstein field equations give us the metric---which then gives us the path of motion by the geodesic equation.

== Newtonian limit
We want to check if we can recover Newton's law of gravity in the Newtonian limit. In this limit all velocities are small $abs(dd(arrow(x))\/dd(tau)) << 1$, and $g_(mu nu)$ is time-independent---so the field is static.

In this limit we have to lowest order (by the first assumption)
$
  dv(x^lambda, tau, 2) + tensor(Gamma, lambda, -mu nu) dv(x^mu, tau) dv(x^nu, tau) = 0 => dv(x^mu, tau, 2) + tensor(Gamma, mu, -00) (dv(t, tau))^2 tilde.eq 0
$
by the second assumption
$
  tensor(Gamma, mu, -00) = 1/2 g^(mu sigma) [pdv(g_(0 sigma), x^0) + pdv(g_(mu 0), x^0) - pdv(g_(00), x^sigma)] tilde.eq - 1/2 g^(mu sigma) pdv(g_(00), x^sigma)
$
we assume gravity is small, so we write
$
  g_(mu nu) (x) = eta_(mu nu) + h_(mu nu) (x)",  " abs(h_(mu nu)) << 1
$
so to leading order in $h_(mu nu)$ we have
$
  tensor(Gamma, mu, -00) tilde.eq - 1/2 eta^(mu sigma) pdv(h_(00), x^sigma)
$
given $h_(00)$ is time-independent then $tensor(Gamma, 0, -00) = 0$, so we find
$
  dv(t, tau) = "const"
$
for $mu = i$ then
$
  dv(arrow(x), tau, 2) - 1/2 (dv(t, tau))^2 grad h_(00) (arrow(x)) tilde.eq 0 => dv(arrow(x), t, 2) = 1/2 grad h_(00) (arrow(x))
$
by comparison with Newton's law of gravity we identify
$
  h_(00) = - 2 Phi + "const"
$
or
$
  g_(00) = - (1 + 2 Phi(x))
$
this immediately leads to non-trivial results by considering a clock in free fall
$
  dd(tau)^2 = - eta_(00) dd(t)^2_"falling" tilde.eq - g_(00) dd(t)^2
$
so $dd(tau) = sqrt(-g_(00)) dd(t)$. This leads to a redshift
$
  omega_2/omega_1 = sqrt((g_00 (x_2))/(g_00 (x_1))) => (Delta omega)/omega_1 = Phi(x_2) - Phi(x_1)
$

#pagebreak()
= The principle of general covariance
We'd like to express the laws of physics such that they hold in all coordinate frames---we want them to be generally covariant. The idea is that if we can express the laws of special relativity in covariant form, then they will be valid in all coordinate frames---in this way we get laws of physics that by definition are valid in every frame and thus satisfy the equivalence principle.

Using the equivalence principle, we can to each point define a freely falling elevator, wherein the laws of physics are those of special relativity. Recall
$
  dv(x^lambda, tau, 2) + tensor(Gamma, lambda, -mu nu) dv(x^mu, tau) dv(x^nu, tau) = 0
$
at every point $x = tilde(x)$ we can pick a freely falling elevator such that
$
  g_(mu nu) (tilde(x)) = eta_(mu nu)
$
for this to hold---i.e. no deviation from free fall--- we need the Christoffel symbols to vanish at $x = tilde(x)$ meaning that
$
  pdv(g_(mu nu) (x), x^sigma)_(x=tilde(x)) = 0
$
but we don't require that
$
  pdv(g_(mu nu) (x), x^sigma, x^rho)_(x=tilde(x)) = 0
$
this means that if we move away from $x=tilde(x)$ the Christoffel symbols don't necessarily vanish anymore.

== Tensors in relativity
To write our laws of physics in a systematic way to ensure covariance we need the language of tensors---since these objects transform nicely under general coordinate transformations.

Scalars like $dd(tau)$ or $phi(x)$ are invariant under $x -> x'$.

Note that $dd(x^mu)$ transforms like
$
  dd(x'^mu) = pdv(x'^mu, x^nu) dd(x^nu)
$
we say that stuff that transforms like $dd(x^mu)$ is a contravariant vector---e.g. so $U^mu$ is a contravariant vector if under $x -> x'$ we have
$
  U'^mu = pdv(x'^mu, x^nu) U^nu
$

Covariant vectors are defined as those that transform inversely of this, so a covariant vector transforms as
$
  A'_mu = pdv(x^nu, x'^mu) A_nu
$
this is the inverse since a contraction of a contravariant and covariant vector gives a scalar
$
  A'_mu U'^mu = pdv(x^nu, x'^mu) pdv(x'^mu, x^sigma) A_nu U^sigma = tensor(delta, nu, -sigma) A_nu U^sigma = A_nu U^nu
$
we can also form a covariant vector by differentiation of a scalar
$
  pdv(phi' (x'), x'^mu) = pdv(phi(x), x^nu) pdv(x^nu, x'^mu)
$
so $A_mu = pdv(phi(x), x^mu) = partial_mu phi(x)$ is covariant. A general tensor transforms in the obvious way. As an example take the metric tensor $g_(mu nu)$ which is covariant
$
  g'_(mu nu) (x') = eta_(alpha beta) pdv(y^alpha, x'^mu) pdv(y^beta, x'^nu) = eta_(alpha beta) pdv(y^alpha, x^rho) pdv(y^beta, x^sigma) pdv(x^rho, x'^mu) pdv(x^sigma, x'^nu) = pdv(x^rho, x'^mu) pdv(x^sigma, x'^nu) g_(rho sigma) (x)
$
similarly $g^(mu nu)$ is contravariant, and we can define
$
  tilde(T)_(mu nu) = g_(mu sigma) g_(nu rho) T^(sigma rho)
$
by comparison we see that $tilde(T)_(mu nu)$ must be covariant since the indices $sigma, rho$ "cancel" leaving just the covariant $mu, nu$, and we let $tilde(T)_(mu nu) = T_(mu nu)$. This also shows that the Kronecker-delta is a mixed tensor $tensor(delta, mu, -nu) = g^(mu sigma) g_(sigma nu)$.

We define $g equiv det g_(mu nu)$ giving $g' = abs(pdv(x, x'))^2 g$ where $abs(pdv(x, x')) = abs(det(pdv(x^mu, x'^nu)))$ is the Jacobian determinant. Then
$
  sqrt(-g') dd(x', 4) = sqrt(-g') abs(pdv(x', x)) dd(x, 4) = sqrt(-g) dd(x, 4)
$
so the measure $sqrt(-g) dd(x, 4)$ is invariant.

== Derivatives
In general the derivative of a tensor is a non-tensor, this is a problem. Instead consider some $x^nu (tau)$ then
$
  "invariant" = dv(phi(x), tau) = pdv(phi, x^mu) (dv(x^mu, tau))
$
implying
$
  "invariant" = dv(phi, tau, 2) = pdv(phi, x^mu, x^nu) dv(x^nu, tau) dv(x^mu, tau) + pdv(phi, x^mu) dv(x^mu, tau, 2)
$
if $x^mu (tau)$ is the trajectory of a free-fall then using the geodesic equation gives
$
  dv(phi, tau, 2) = (pdv(phi, x^mu, x^nu)-tensor(Gamma, sigma, -mu nu) pdv(phi, x^sigma)) times dv(x^mu, tau) dv(x^nu, tau)
$
the right factor is contravariant, so the left factor must be covariant given their contraction is invariant. With $V_mu = pdv(phi, x^mu)$ we define the covariant derivative
$
  D_mu V_nu = partial_mu V_nu - tensor(Gamma, sigma, -mu nu) V_sigma
$
similarly
$
  D_mu V^nu = partial_mu V^nu + tensor(Gamma, nu, -mu sigma) V^sigma
$
which generalizes in the obvious way. In the freely falling local elevator this reduces to the ordinary derivative since the Christoffel symbol will vanish---so $D_mu -> partial_mu$. This also implies that the Christoffel symbol is not a tensor.

In a local elevator we have
$
  partial_sigma g_(mu nu) = 0 => D_sigma g_(mu nu) = 0
$
since the covariant form reduces to the special relativity equation---$D_sigma -> partial_sigma$ so the covariant form holds in the local elevator, and thus it must hold in any frame---similarly any equation can be made generally covariant by $partial_sigma -> D_sigma$ since $D_sigma -> partial_sigma$ in the local elevator.

By projecting the covariant derivative to the tangent of the curve $x^mu (tau)$ by multiplying it by $dd(x^mu)\/dd(tau)$ we obtain
$
  dv(V^nu, tau, d: D) equiv dv(V^nu, tau) + tensor(Gamma, nu, -mu sigma) V^sigma dv(x^mu, tau)
$
notice that $D_tau -> d_tau$ in the local elevator, this lets us derive the geodesic equation. In special relativity the equation of motion can be written as
$
  dv(V^mu, tau) = 0 => dv(V^mu, tau, d: D) = 0
$
which using $V^mu = dv(x^mu, tau)$ becomes the geodesic equation.

== Electromagnetism
From electrodynamics we recall the antisymmetric Faraday tensor $ F_(mu nu) = partial_mu A_nu - partial_nu A_mu $ and the four-current $J^mu = vecrow(rho, bold(J))$. Then we know that the two inhomogeneous Maxwell equations can be written in compact form as
$
  partial_mu F^(mu nu) = - J^nu
$
likewise the two homogeneous Maxwell equations become
$
  partial_mu F_(nu gamma) + partial_(gamma) F_(mu nu) + partial_nu F_(gamma mu) = 0
$
which we call the Bianchi identity---this generally holds for tensors of the form $T_(mu nu) = partial_mu B_nu - partial_nu B_mu$ since terms cancel due to antisymmetry and the Schwarz identity. Two include gravity we simply replace $partial_mu -> D_mu$ giving
$
  D_mu F^(mu nu) & = - J^nu \
               0 & = D_mu F_(nu lambda) + D_lambda F_(mu nu) + D_nu F_(lambda mu)
$
But, due to the cyclical nature of the Bianchi identity and the antisymmetricy of $F_(mu nu)$, then all Christoffel symbols due to the covariant derivative cancel giving
$
  partial_mu F_(nu lambda) + partial_lambda F_(mu nu) + partial_nu F_(lambda mu) = 0
$
so it is unchanged. Also using the identity
$
  D_mu F^(mu nu) = 1/sqrt(-g) partial_mu (sqrt(-g) F^(mu nu))
$
which holds for any antisymmetric tensor---we obtain
$
  partial_mu (sqrt(-g) F^(mu nu)) = - sqrt(-g) J^nu
$
differentiating both sides the LHS gives
$
  partial_nu partial_mu (sqrt(-g) F^(mu nu)) = 0
$
since $partial_nu partial_mu$ is symmetric while $F^(mu nu)$ is antisymmetric. This implies
$
  partial_nu (sqrt(-g) J^nu) = 0
$
which is the covariant version of four-current conversation, as opposed to in special relativity where $partial_nu J^nu = 0$.

Finally recall the Lorentz force
$
  f^mu = e tensor(F, mu, -nu) U^nu
$
this is already covariant! Given we write $F_(mu nu) = D_mu A_nu - D_nu A_mu$.


