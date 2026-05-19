#import "chpt-temp.typ": *
#show: chpt-note.with()

= Fields in Curved Space
== The scalar field
We consider a scalar field in curved space. The Lagrangian for such a field is
$
  cal(L) = - 1/2 sqrt(-g) [g^(mu nu) partial_mu phi.alt partial_nu phi.alt + (m^2 + xi R) phi.alt^2]
$
with $xi$ being the _non-minimal coupling_. We have the metric
$
  dd(s^2) = g_(mu nu) dd(x^mu, x^nu)";  " mu, nu = 0,1 dots, d-1
$
and we define
$
  g equiv det g_(mu nu)
$
The term $xi R$ is the only possible local scalar coupling in $cal(L)$ by dimensional analysis. The inclusion of other, higher dimensional couplings would require the introduction of some mass scale.

By varying the action
$
  S = integral dd(x, [n]) cal(L)
$
we trivially find the equation of motion#footnote[This is simply the Klein-Gordon equation with the mass $m$ shifted.]
$
  (square + m^2 + xi R) phi.alt = 0
$
and we can write#footnote[Using $(-,+,+,+)$.]
$
  square phi.alt & equiv g^(mu nu) nabla_mu nabla_nu phi.alt \
  & eq 1/sqrt(-g) partial_mu (sqrt(-g) g^(mu nu) partial_nu phi.alt)
$
We say the field is _minimally coupled_ since when $xi -> 0$ we recover flat space. When
$
  xi = 1/4 (n-2)/(n-1) tilde^(n=4) 1/6
$
The field has _conformal coupling_ since the equations of motion with $m = 0$ become invariant under _conformal transformations_ $g_(mu nu) -> Omega^2 g_(mu nu)$.

We seek an orthonormal basis for the solutions to the above equation. This means we need a covariant definition of the inner product. We define
$
  (phi.alt_1, phi.alt_2) = - i integral_Sigma phi.alt_1 (bold(x)) accent(nabla, \u{2194})_mu phi.alt_2^* (bold(x)) sqrt(-g_Sigma) dd(Sigma^mu)
$
with $dd(Sigma^mu) = n^mu dd(Sigma)$ where $n^mu$ is the _future-directed_ unit vector orthogonal to the _spacelike_ hypersurface $Sigma$ and
$
  phi.alt_1 accent(nabla, \u{2194})_mu phi.alt_2^* equiv phi.alt_1 nabla_mu phi.alt_2^* - phi.alt_2^* nabla_mu phi.alt_1 tilde j_mu
$
Consider
$
  D_mu j^mu &= nabla_mu (phi.alt_1 nabla_mu phi.alt_2^* - phi.alt_2^* nabla_mu phi.alt_1) \
  &= phi.alt_1 square phi.alt_2^* - phi.alt_2^* square phi.alt_1 \
  &= phi.alt_1 (-m^2 - xi R) phi.alt_2^* - phi.alt_2^* (-m^2 - xi R) phi.alt_1 \
  &= 0
$
meaning $j_mu$ is conserved. We now take two spacelike hypersurfaces $Sigma_1$ and $Sigma_2$ that bound a spacetime region $cal(V)$. Then by Gauss' theorem
$
  integral_cal(V) dd(x, [d]) sqrt(-g) D_mu j^mu = integral_(dd(cal(V), d: partial)) dd(Sigma_mu) j^mu
$
The LHS vanishes so
$
  0 &= integral_Sigma_2 dd(Sigma_mu) j^mu - integral_Sigma_1 dd(Sigma_mu) j^mu + underbracket("boundary terms", -> 0)
$
implying
$
  integral_(Sigma_1) dd(Sigma_mu) j^mu = integral_(Sigma_2) dd(Sigma_mu) j^mu tilde Q(Sigma)
$
meaning $(phi.alt_1,phi.alt_2)$ is independent of the choice of $Sigma$.
With Minkowski spacetime we recover the usual _Klein-Gordon inner product_
$
  (phi.alt_1, phi.alt_2) =^"flat" -i integral_t phi.alt_1 accent(partial_t, \u{2194}) phi.alt_2^* dd(bold(x), 3)
$
The orthonormal modes satisfy
$
  (u_i, u_j) = delta_(i j)";  " (u_i^*, u_j^*) = - delta_(i j)";  " (u_i, u_j^*) = 0
$
which are called the _Wronskian conditions_. Consider the Minkowski modes
$
  u_bold(k) tilde e^(plus.minus i bold(k) dot bold(x) minus.plus i omega t)
$
These satisfy the Wronskian conditions meaning the above definition is a valid covariant generalisation of the inner product.#footnote[By the _principle of general covariance_.]

Then some general solution can be written as
$
  phi.alt (bold(x)) = sum_i (a_i u_i (bold(x)) + a_i^dagger u_i^* (bold(x)))
$
Quantizing the field is done as before by imposing the commutation relations
$
  [a_i, a_j^dagger] = delta_(i j)";   others" = 0
$

== Bogolubov transformations
The Minkowski metric is time-independent implying the existence of a $pdv(, t)$ _Killing vector_ orthonormal to the spacelike hypersurface where $t$ is constant. The Minkowski modes are eigenfunctions of $pdv(, t)$ with eigenvalues $-i omega$. We have something similar for spatial translations and boosts and the Minkowski vacuum is the state invariant under the Poincaré group. This defines a _unique_ vacuum. However in curved space Poincaré invariance breaks meaning there is no unique vacuum in general. The existence of Killing vectors can amend this.#footnote[We can view this as general covariance making coordinates physically meaningless. This implies no unique mode decomposition exists and as a consequence there is no well-defined frequency or conserved energy.]

We want to compare two different orthonormal sets $u_i$ and $overline(u)_j$ and the resulting equivalent expansions
$
  phi.alt = sum_i (a_i u_i + a_i^dagger u_i^*)";  " phi.alt = sum_j (overline(a)_j overline(u)_j + overline(a)_j^dagger overline(u)_j^*)
$
We are lead to two definitions of the vacuum
$
  a_i ket(0) = 0";  " overline(a)_j ket(overline(0)) = 0
$
The bases are assumed to be complete implying we can write
$
  overline(u)_j = sum_i (alpha_(j i) u_i + beta_(j i) u_i^*)";  " u_i = sum_j (alpha_(j i)^* overline(u)_j - beta_(j i) overline(u)_j^*)
$
These are called the _Bogolubov transformations_ with $alpha_(i j)$ and $beta_(i j)$ being the _Bogolubov coefficients_. We immediately find
$
  alpha_(i j) = (overline(u)_i, u_j)";  " beta_(i j) = - (overline(u)_i, u_j^*)
$
Then we have
$
  sum_i (a_i u_i + a_i^dagger u_i^*) = sum_i (overline(a)_i overline(u)_i + overline(a)_i^dagger overline(u)_i^*)
$
Taking $(dots, u_i)$ we obtain
$
  a_i &= sum_j [overline(a)_j (overline(u)_j,u_i) + overline(a)_j^dagger (overline(u)_j^*, u_i)] = sum_j (alpha_(j i) overline(a)_j + beta_(j i)^* overline(a)_j^dagger)
$
similarly
$
  overline(a)_j = sum_i (alpha_(j i) a_i + beta_(j i)^* a_i^dagger)
$
We see that for $beta_(i j) eq.not 0$
$
  a_i ket(overline(0)) = sum_j beta_(j i)^* ket(overline(1)_j) eq.not 0
$
meaning the vacua are _truly different_. We also have for $N_i = a_i^+ a_i$
$
  braket(overline(0), N_i, overline(0)) = sum_j abs(beta_(j i))^2
$
implying the vacuum $ket(overline(0))$ contains particles of the $u_i$ mode!#footnote[This has many consequences. One such example is structure formation due to particle production when the vacuum changes with time in an expanding Universe. Another example is the _Unruh effect_ where an accelerated observer observes particles in the Minkowski vacuum. We will discuss these in due time.]

== Killing vectors
We consider an integral curve $x^mu (t)$ and define the vector field
$
  V^mu (x) equiv evaluated(dv(x^mu (t), t))_x tilde "tangents"
$
The _Lie derivative_ $scr(L)_V$ measures how any field changes as we drag it along the flow generated by $V^mu$. We simply define#footnote[See e.g. _Nakahara_ for details.]
$
     scr(L)_V f & = V^mu partial_mu f \
  scr(L)_V U^mu & = V^nu partial_nu U^mu - U^nu partial_nu V^mu
$
The Lie derivative of $g_(mu nu)$ is
$
  scr(L)_V g_(mu nu) & = nabla_mu V_nu + nabla_nu V_mu
$
To see this consider a displacement along the direction given by some vector $X^mu$
$
  x'^mu = x^mu + epsilon X^mu
$
Then
$
  g_(mu nu) &= (delta_mu^rho + epsilon partial_mu X^rho)(delta_nu^sigma + epsilon partial_nu X^sigma) underbracket(g_(rho sigma) (x^lambda + epsilon X^lambda), g_(rho sigma) (x) + epsilon X^lambda partial_lambda g_(rho sigma) (x)+ dots) \
  &= g_(mu nu) (x) + epsilon (g_(mu sigma) partial_nu X^sigma + g_(nu sigma) partial_mu X^sigma + X^lambda partial_lambda g_(mu nu)) + cal(O) (epsilon^2) \
  &= g_(mu nu) + epsilon scr(L)_X g_(mu nu)
$
We say $xi^mu equiv X^mu$ is a _Killing vector_ if the above operation is an _isometry_ meaning
$
  scr(L)_xi g_(mu nu) = nabla_mu xi_nu + nabla_nu xi_mu = 0
$
which is the _Killing equation_. This implies conserved quantities since $g_(mu nu)$ is invariant along the _flow_ generated by $xi_mu$.

We can explicitly pick coordinates $vecrow(y, x^1, x^2, x^3)$ where
$
  xi^mu = (partial_y)^mu tilde "basis vector along" y
$
or $xi^mu = delta^mu_y$. This collapses the definition of the Lie derivative as
$
  scr(L)_xi -> partial_y
$
implying $partial_y g_(mu nu) = 0$. Consider the case mentioned before where $y = t$ meaning $xi = partial_t$ and
$
  partial_t g_(mu nu) = 0 tilde "stationary spacetime"
$
which is nice. Also since $nabla_mu T^(mu nu) = 0$ we have
$
  nabla_mu (xi_nu T^(mu nu)) = 1/sqrt(-g) partial_mu (sqrt(-g) xi_nu T^(mu nu)) = 0
$
meaning the current $j^mu = xi_nu T^(mu nu)$ is conserved. Assuming the spacetime is stationary we have $xi^mu = delta^mu_t$ implying
$
  nabla_mu tensor(T, mu, -t) = 0 tilde "local energy conservation"
$
Likewise $xi = partial_phi.alt$ gives angular momentum conservation.
