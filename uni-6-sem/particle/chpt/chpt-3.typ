#import "../../temp.typ": *
#show: chpt-note.with()

= Special Relativity
== The Lorentz transformations
According to special relativity the laws of physics are valid in all inertial frames.#footnote[By inertial we mean non-accelerating.] Now, consider two inertial frames $S$ and $S'$, with $S'$ having a uniform velocity $bold(v)$ relative to $S$. With coordinates where $bold(v)$ is along $hat(x)$ these frames are related by
$
  x' = gamma (x-beta c t)";  " y' = y";  " z' = z";  " c t' = gamma (c t - beta x),
$
with
$
  gamma & eq 1/sqrt(1-beta^2)",  " beta & eq v/c.
$
We call these the Lorentz transformations. #footnote[The inverse transformations are given by $beta -> - beta$]
== Four-vectors
We define the position four-vector $x^mu$
$
  x^mu = vecrow(c t, x, y, z)
$
By comparison with above the Lorentz transformations become
$
  x'^0 = gamma (x^0 - beta x^1)",  etc."
$
We write $x'^mu = tensor(Lambda, mu, -nu) x^nu$#footnote[Using the summation convention.] where $tensor(Lambda, mu, -nu)$ represents the Lorentz transformations for a given $bold(v)$.
We define a general four-vector $A^mu$ as an object transforming like $x^mu$
$
  A^mu = tensor(Lambda, mu, -nu) A^nu
$
We would like invariant quantities. These can be found by contractions#footnote[All indices cancel.] using the metric $g_(mu nu)$
$
  "scalar" = g_(mu nu) A^mu A^nu equiv A_mu A^mu equiv A^2
$
We work in flat Minkowski spacetime where $g_(mu nu) = eta_(mu nu)$.#footnote[We use $eta_(mu nu) = diag(1, -1, -1, -1)$.] The spacetime interval is the most important example
$
  dd(s^2) = g_(mu nu) dd(x^mu, x^nu) = c^2 dd(t^2) - dd(x^2) - dd(y^2) - dd(z^2)
$
We can extend the definition of a four-vector to define general second-rank tensors $T^(mu nu)$ as objects transforming like
$
  T'^(mu nu) = tensor(Lambda, mu, -kappa) tensor(Lambda, nu, -sigma) T^(kappa sigma)
$
The definition of higher rank tensors is obvious. We can raise and lower indices as implied before using $g_(mu nu)$
$
  tensor(T, mu, -nu) = g_(nu kappa) tensor(T, mu, kappa) tilde "mixed"
$
We can also use $g_(mu nu)$ to take the trace
$
  tensor(T, mu, -mu) = g_(mu nu) tensor(T, mu, nu) tilde "scalar"
$

== The four-momentum
We define the proper time $tau$ by
$
  dd(tau) = dd(t)/gamma
$
This is the time measured in the instantaneous rest-frame.#footnote[We use $dd(s^2) = c^2 dd(tau^2) = c^2 dd(t^2) - dd(bold(r)^2)$.]

We define velocity as
$
  bold(v) equiv dv(bold(x), t) tilde "lab frame"
$
Then the proper velocity is defined by
$
  bold(eta) equiv dv(bold(x), tau) = gamma bold(v)
$
This quantity is useful since $dd(tau)$ is a scalar. This motivates the four-velocity $eta^mu$
$
  eta^mu equiv dv(x^mu, tau)
$
The related invariant is $eta_mu eta^mu = c^2$. We can now define the four-momentum $p^mu$
$
  p^mu equiv m eta^mu
$
With $E equiv gamma m c^2$ we can write $p^mu$ as
$
  p^mu = vecrow(E/c, bold(p))
$
with $bold(p) = gamma m bold(v)$. The related invariant is
$
  p_mu p^mu = E^2/c^2 - bold(p)^2 =^"rest frame" m^2 c^2
$
We can expand $E$ to obtain
$
  E = underbracket(m c^2, "rest energy") + underbracket(1/2 m v^2 + dots, "kinetic energy")
$
When we consider particles with $m = 0$ we have $p_mu p^mu = 0$ implying
$
  E = abs(bold(p))c
$

We care about $p^mu$ since it is always conserved in any physical process!

== Examples
Consider two lumps with mass $m$ colliding at $3/5 c$. After the collision the lumps stick together.

We obviously have $bold(p)_1 = - bold(p)_2$ meaning $bold(p)_M = 0$. Then since $E_1 + E_2 = E_M$ we have
$
  M c^2 = 2 E_1 = (2 m c^2)/sqrt(1- (3/5)^2) = 5/2 m c^2
$
implying
$
  M = 5/2 m > 2 m
$

A particle of mass $M$ at rest decays into two pieces with mass $m$.

By $bold(p)_M = bold(p)_1 + bold(p)_2 = 0$ we have $bold(p)_1 = -bold(p)_2$. Then since $E_M = 2 E_1$ we have
$
  M = (2 m)/sqrt(1 - beta^2) => beta^2 = 1 - ((2 m)/M)^2
$
We need $M >= 2 m$ for the process to make sense. We refer to $M = 2m$ as the threshold.#footnote[When $beta = 0$.]

A $pi$ at rest decays into $mu + nu$.

We have
$
  p_pi = p_mu + p_nu => p_nu = p_pi - p_mu
$
implying
$
  p_nu^2 = p_pi^2 + p_mu^2 - 2 p_pi dot p_mu
$
Where $p_nu^2 tilde 0$ and
$
  p_pi^2 = m_pi^2 c^2";  " p_mu^2 = m_mu^2 c^2";  " p_pi dot p_mu = (E_pi E_mu)/c^2
$
Then
$
  0 = m_pi^2 c^2 + m_mu^2 c^2 - 2 m_pi E_mu => E_mu tilde.eq (m_pi^2 + m_mu^2)/(2 m_pi) c^2
$
Similarly with $p_mu = p_pi - p_nu$ we find
$
  m_mu^2 c^2 = m_pi^2 c^2 - 2 m_pi E_nu
$
using $E_nu tilde.eq abs(bold(p)_mu) c$ we have
$
  abs(bold(p)_mu) tilde.eq (m_pi^2 - m_mu^2)/(2 m_pi) c
$

A $p$ collides with a $p$ at rest by $p + p -> p + p +p + overline(p)$. What is the threshold energy?

We define the center-of-momentum frame to be the frame where $bold(p)_"total" = 0$. We consider $p^mu_"total"$ before the collision in the lab-frame
$
  p^mu_"total" = vecrow((E+m c^2)/c, abs(bold(p)), 0, 0)
$
and the $p'^mu_"total"$ after the collision in the CM-frame
$
  p'^mu_"total" = vecrow(4 m c, bold(0)) tilde 4 p "at rest"
$
Then we use $p_mu p^mu = p'_mu p'^mu$ to obtain
$
  E = 7 m c^2
$
Consider some complicated frame $S$ with
$
  E_"total" = sum_i gamma_i m_i c^2";  " bold(p)_"total" = sum_i gamma_i m_i bold(v)_i
$
Then
$
  abs(bold(p)'_"total") = gamma (abs(bold(p)_"total") - beta E_"total"/c)
$
The CM-frame $S'$ is then given by
$
  beta = (abs(bold(p)_"total") c)/E_"total"
$
We always have $beta < 1$ meaning this frame always exists.
