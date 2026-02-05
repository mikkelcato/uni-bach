#import "chpt-temp.typ": *
#show: chpt-note.with()

= Note on General Relativity
The object of interest in general relativity is the _metric_ field $g_(mu nu)$. To determine it we solve the _Einstein field equations_
$
  underbracket(R_(mu nu) - 1/2 g_(mu nu) R, "geometry of spacetime") = 8 pi G_N underbracket(T_(mu nu), "matter")
$
where $T_(mu nu)$ is the _energy-momentum tensor_ and $R_(mu nu)$ is the _Ricci tensor_ with $R equiv g^(mu nu) R_(mu nu)$ being the _Ricci scalar_. The Ricci tensor is defined by
$
  R_(mu kappa) equiv g^(lambda nu) R_(lambda mu nu kappa)
$
with
$
  tensor(R, sigma, -mu nu kappa) equiv partial_kappa tensor(Gamma, sigma, -mu nu) - partial_nu tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu nu) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -mu kappa) tensor(Gamma, sigma, -nu lambda)
$
being the _Riemann-Christoffel curvature tensor_. Where
$
  tensor(Gamma, lambda, -mu nu) equiv 1/2 g^(lambda sigma) (partial_mu g_(sigma nu) + partial_nu g_(sigma mu) - partial_sigma g_(mu nu))
$
is the _Christoffel symbol_ or _affine connection_.

== The Einstein-Hilbert action
We would like some action
$
  S_H = integral dd(x, 4) sqrt(-g) cal(L)_H
$
with $cal(L)_H$ being a scalar and $g equiv det g_(mu nu)$. The simplest choice is $cal(L)_H = R$. We minimize this with respect to $g_(mu nu)$ by $dd(S_H, d: delta)=0$ giving
$
  1/sqrt(-g) dv(S_H, g^(mu nu), d: delta) = R_(mu nu) - 1/2 R g_(mu nu) = 0
$
To include matter we add $S_M$
$
  S = 1/(8 pi G_N) S_H + S_M
$
We obtain
$
  1/(sqrt(-g)) dv(S, g^(mu nu), d: delta) = 1/(8 pi G_N) [R_(mu nu) - 1/2 g_(mu nu) R] + 1/sqrt(-g) dv(S_M, g^(mu nu), d: delta) =^! 0
$
Leading us to the definition of the energy-momentum tensor
$
  T_(mu nu) = - 1/sqrt(-g) dv(S_M, g^(mu nu), d: delta)
$
As a final remark we can always add another constant to $cal(L)_H$ by $R -> R-2 Lambda$ giving an additional term
$
  R_(mu nu) - 1/2 g_(mu nu) R + Lambda g_(mu nu) = 0
$
We call $Lambda$ the cosmological constant. We could also treat this as the energy density of vacuum by
$
  R_(mu nu) - 1/2 g_(mu nu) R = underbracket(- Lambda g_(mu nu), T_(mu nu))
$

#pagebreak()
= Canonical Quantization
== Classical field theory
Before attempting to understand quantum field theory on a curved background we should understand quantum field theory on a flat background.


We consider the theory of a real scalar field $phi.alt(bold(x), t)$. This is described by
$
  cal(L) = 1/2 partial_mu phi.alt partial^mu phi.alt - 1/2 m^2 phi.alt^2
$
Which gives the _Klein-Gordon equation_
$
  (partial_mu partial^mu + m^2) phi.alt(bold(x), t) = 0
$
To quantize our field we need the Hamiltonian formulation. We define the _canonical conjugate field_ $pi$ by
$
  pi (bold(x),t) equiv pdv(cal(L), dot(phi.alt)(bold(x),t)) = dot(phi.alt)(bold(x),t)
$
Then we can write
$
  H & = integral dd(x, 3) [pi (bold(x),t) dot(phi.alt)(bold(x),t)-cal(L)] \
    & = 1/2 integral dd(x, 3) [ pi^2 + (nabla phi.alt)^2 + m^2 phi.alt^2]
$
Which is just the Legendre transform.

== Quantization
We promote $phi.alt$ and $pi$ to operators, and impose the _continuous_ generalization of the usual commutation relations
$
  [phi.alt(bold(x), t), pi(bold(y), t)] & = i delta^((3)) (bold(x)-bold(y)) \
  [phi.alt(bold(x), t),phi.alt(bold(y), t)] & = [pi(bold(x), t),pi(bold(y), t)] = 0
$
We would like to determine the spectrum of our Hamiltonian. How we proceed is in no way obvious. We try expanding $phi.alt$ as
$
  phi.alt (bold(x),t) = integral dd(bold(k), 3)/(2pi)^(3\/2) e^(i bold(k) dot bold(x)) phi.alt_bold(k) (t)
$
with $phi.alt_bold(k)^* = phi.alt_(-bold(k))$. Then the Klein-Gordon equation becomes
$
  [pdv(, t, 2) + (bold(k)^2+m^2)] phi.alt_bold(k) (t) = 0
$
implying each mode becomes a harmonic oscillator with frequency $omega_bold(k) equiv sqrt(bold(k)^2 + m^2)$! We can now write our Hamiltonian as
$
  H = 1/2 integral dd(k, 3) [pi_bold(k) pi_(-bold(k)) + omega_bold(k)^2 phi.alt_bold(k) phi.alt_(-bold(k))]
$
Since every mode is a harmonic oscillator it is convenient to introduce ladder operators
$
  a_bold(k) (t) &equiv sqrt(omega_bold(k)/2) (phi.alt_bold(k) + i pi_bold(k)/omega_bold(k)) \
  a_bold(k)^dagger (t) &equiv sqrt(omega_bold(k)/2) (phi.alt_(-bold(k))- i pi_(-bold(k))/omega_bold(k))
$
The equations of motion for $phi.alt$ and $pi$ become (since $dot(phi.alt)_bold(k) = pi_bold(k)$)
$
  dv(, t) a_bold(k)^dagger = i omega_bold(k) a_bold(k)^dagger";  " dv(, t) a_bold(k) = - i omega_bold(k) a_bold(k)
$
Which has the solution
$
  a_bold(k)^dagger (t) = a_bold(k)^dagger (0) e^(i omega_bold(k) t)";  " a_bold(k) (t) = a_bold(k) (0) e^(-i omega_bold(k) t)
$
Where the ladder operators must satisfy
$
          [a_bold(k), a^dagger_bold(k')] & = delta^((3)) (bold(k)-bold(k)') \
  [a^dagger_bold(k),a^dagger_(bold(k)')] & = [a_bold(k),a_bold(k)'] = 0
$
We can now write
$
  phi.alt_bold(k) &= 1/sqrt(2 omega_bold(k)) (a_bold(k) e^(-i omega_bold(k) t) + a_(-bold(k))^dagger e^(i omega_bold(k) t)) \
  pi_bold(k) &= i sqrt(omega_bold(k)/2) (a_(-bold(k))^dagger e^(i omega_bold(k) t) - a_bold(k) e^(-i omega_bold(k) t))
$
Then our Hamiltonian becomes
$
  H &= integral dd(k, 3) omega_bold(k)/2 [a_bold(k) a^dagger_bold(k) + a_bold(k)^dagger a_bold(k)] \
  &= integral dd(k, 3) [omega_bold(k) a_bold(k)^dagger a_bold(k) + underbracket(omega_bold(k)/2 delta^((3))(0), "vacuum energy" #linebreak() -> oo)]
$
We can easily compute
$
  [H,a_bold(k)^dagger] = omega_bold(k) a_bold(k)^dagger";  " [H,a_bold(k)] = - omega_bold(k) a_bold(k)
$
Vacuum is defined by $a_bold(k) ket(0) = 0$ for all $bold(k)$. Any excited state can then be built by acting with $a^dagger$ on $ket(0)$. To see this consider
$
  H (a_bold(k)^dagger ket(0)) = omega_bold(k) a_bold(k)^dagger ket(0) + underbracket(a_bold(k)^dagger H ket(0), "dropping infinity")
$
This is consistent with interpreting $ket(bold(k)) = a_bold(k)^dagger ket(0)$ as a one-particle state with momentum $bold(k)$ and energy $E = omega_bold(k)$.

We can now express $phi.alt (bold(x),t)$ as
$
  phi.alt (bold(x),t) &= integral dd(k, 3)/(2 pi)^(3\/2) 1/sqrt(2 omega_bold(k)) [a_bold(k) e^(-i omega_bold(k) t+i bold(k) dot bold(x))+ a_bold(k)^dagger e^(i omega_bold(k) t - i bold(k) dot bold(x))]
$
Typically this form is written immediately (with commutation relations for $a$ and $a^dagger$) to _shortcut_ quantization.
