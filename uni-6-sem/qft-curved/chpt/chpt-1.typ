#import "chpt-temp.typ": *
#show: chpt-note.with()

= Note on General Relativity
The object of interest in general relativity is the _metric_ field $g_(mu nu)$. To determine it we solve the _Einstein field equations_
$
  underbracket(cal(R)_(mu nu) - 1/2 g_(mu nu) cal(R), "geometry of spacetime") = 8 pi G_N underbracket(T_(mu nu), "matter")
$
where $T_(mu nu)$ is the _energy-momentum tensor_ and $cal(R)_(mu nu)$ is the _Ricci tensor_ with $cal(R) equiv g^(mu nu) cal(R)_(mu nu)$ being the _Ricci scalar_. The Ricci tensor is defined by
$
  cal(R)_(mu kappa) equiv g^(lambda nu) cal(R)_(lambda mu nu kappa)
$
with
$
  tensor(cal(R), sigma, -mu nu kappa) equiv partial_kappa tensor(Gamma, sigma, -mu nu) - partial_nu tensor(Gamma, sigma, -mu kappa) + tensor(Gamma, lambda, -mu nu) tensor(Gamma, sigma, -kappa lambda) - tensor(Gamma, lambda, -mu kappa) tensor(Gamma, sigma, -nu lambda)
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
with $cal(L)_H$ being a scalar and $g equiv det g_(mu nu)$. The simplest choice is $cal(L)_H = cal(R)$. We minimize this with respect to $g_(mu nu)$ by $dd(S_H, d: delta)=0$ giving
$
  1/sqrt(-g) dv(S_H, g^(mu nu), d: delta) = cal(R)_(mu nu) - 1/2 cal(R) g_(mu nu) = 0
$
To include matter we add $S_M$
$
  S = 1/(8 pi G_N) S_H + S_M
$
We obtain
$
  1/(sqrt(-g)) dv(S, g^(mu nu), d: delta) = 1/(8 pi G_N) [cal(R)_(mu nu) - 1/2 g_(mu nu) cal(R)] + 1/sqrt(-g) dv(S_M, g^(mu nu), d: delta) =^! 0
$
Leading us to the definition of the energy-momentum tensor
$
  T_(mu nu) = - 1/sqrt(-g) dv(S_M, g^(mu nu), d: delta)
$
As a final remark we can always add another constant to $cal(L)_H$ by $cal(R) -> cal(R)-2 Lambda$ giving an additional term
$
  cal(R)_(mu nu) - 1/2 g_(mu nu) cal(R) + Lambda g_(mu nu) = 0
$
We call $Lambda$ the cosmological constant. We could also treat this as the energy density of vacuum by
$
  cal(R)_(mu nu) - 1/2 g_(mu nu) cal(R) = underbracket(- Lambda g_(mu nu), T_(mu nu))
$

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

= The Casimir Effect
Assuming vacuum energy acts like any other energy then according to relativity it will gravitate. The vacuum energy density is infinite so this implies an infinite gravitational field! However since we exist something about this is wrong.

Above we found
$
  E_0 & = braket(0, H, 0) \
      & = integral dd(k, 3) omega_bold(k)/2 delta^((3)) (0) tilde "diverges"
$
which is problematic. We might argue the integral has an upper limit since particles with $lambda_"Compton"$ smaller than the Planck length $l_p$ would become black holes since their Schwarzschild radius $r_s$ would be larger than $lambda_"Compton"$. However introducing a UV cutoff $Lambda tilde M_p$ still yields $epsilon_0 tilde M_p^4$ which is enormous. This is called _regularising_ the integral and we see the divergence reappears if we remove the cutoff by $Lambda -> oo$. Typically regularisation is followed by _renomalisation_ which is simply a prescription for how to deal with the large arbitrary cutoff dependent contribution to some physical quantity. Above this quantity is the cosmological constant.

== The effect
The _Casimir effect_ is an experimentally verified prediction of the above procedure which we now derive. Consider a massless scalar field $phi.alt$ in Minkowski $1+1$ space-time. We now consider two plates separated by a distance $L$. We require the field vanishes on the plates
$
  evaluated(phi.alt)_(x=0) = evaluated(phi.alt)_(x=L) = 0
$
and consider solutions to the equation of motion
$
  partial_t^2 phi.alt - partial_x^2 phi.alt = 0
$
These are given by
$
  phi.alt = sum_(n=1)^oo (A_n e^(-i omega_n t) + B_n e^(i omega_n t)) sin omega_n x
$
where $omega_n equiv n pi\/L$.

Without the plates we would have the following mode expansion
$
  phi.alt = integral dd(k)/(sqrt(2 pi)) 1/sqrt(2 omega_k) [a_k e^(-i omega_k t + i k x) + a_k^dagger e^(i omega_k t - i k x)]
$
However the plane wave modes are no longer viable after introducing the plates. The plates lead to a discrete subset of orthonormal functions being the basis for the mode expansion. We introduce the basis
$
  g_n (x) = sqrt(2/L) sin omega_n x";  " integral_0^L g_m (x) g_n (x) dd(x) = delta_(m n)
$
The continuous $k$ has been replaced by the discrete $omega_n$ and now the basis satisfies the boundary conditions. We can now expand the field in terms of $e^(-i omega_n t) g_n (x)$ instead of plane waves to obtain
$
  phi.alt = sqrt(2/L) sum_(n=1)^oo (sin omega_n x)/sqrt(2 omega_n) [a_n e^(-i omega_n t) + a_n^dagger e^(i omega_n t)]
$
We see that fewer modes are allowed between the plates as compared to outside. The only modes allowed are those that form standing waves. This leads to the vacuum energy between the plates being smaller as compared to outside!

== The energy
We want to compute the vacuum energy between the plates. Using
$
  H = 1/2 integral_0^L dd(x) [(partial_t phi.alt)^2 + (partial_x phi.alt)^2]
$
and
$
  braket(0, a_m a_n^dagger, 0) = delta_(m n)";  " "others" = 0
$
We obtain
$
  epsilon_0 & equiv 1/L braket(0, H, 0) = pi/(2 L^2) sum_(n=1)^oo n tilde "diverges"
$
As $L -> oo$ this diverges and
$
  epsilon_0^"free" = lim_(L-> oo) epsilon_0 (L)
$
This energy should not be infinite and should be $tilde 0$! #footnote[Assuming we ignore dark energy which is negligible.] We are led to
$
  dd(epsilon(L), d: Delta) = epsilon_0 (L) - epsilon_0^"free" tilde oo - oo
$
Which is absurd.

We regulate $epsilon_0$ and $epsilon_0^"free"$ by a cutoff hoping that the physical observables we compute become finite and independent of the cutoff. We introduce an exponential cutoff
$
  epsilon_0 (L,alpha) = pi/(2 L^2) sum_(n=1)^oo n exp(- (n alpha)/L)
$
with $alpha$ being the cutoff parameter. We can write$ epsilon_0 (L,alpha) & = -pi/(2 L) pdv(, alpha) sum_(n=1)^oo exp(-(n alpha)/L) \
                    & = pi/(8 L^2) sinh^(-2) (alpha/(2 L)) $
Taking the limit $alpha -> 0$ where the cutoff vanishes we then have
$
  epsilon_0 (L, alpha) = underbracket(pi/(2 alpha^2), tilde "divergent") - pi/(24 L^2) + underbracket(cal(O) (alpha^2), "vanishes")
$
Consider now $dd(epsilon (L), d: Delta) = epsilon_0 (L) - epsilon_0^"free"$ then
$
  dd(epsilon_"renorm", d: Delta) (L) &= lim_(alpha -> 0) [epsilon_0 (L,alpha)- lim_(L-> oo) epsilon_0 (L,alpha)] \
  &= lim_(alpha -> 0) (pi/(2 alpha^2) - pi/(24 L^2) + dots - pi/(2 alpha^2)) \
  &= - pi/(24 L^2)
$
Then the zero-point energy in the presence of plates is lower by a finite amount as compared to the free case!

This vacuum energy is a form of potential energy depending on the distance $L$ between the plates. Then by definition we have a force
$
  F_"Casimir" & = - dv(, L) dd(E, d: Delta) \
              & = - dv(, L) (L dd(epsilon_"renorm", d: Delta)) \
              & = - pi/(24 L^2)
$
This force is negative meaning the plates get pulled together. The increased number of modes outside the plates gives rise to a pressure!

== A note on the cosmological constant
To determine the physical vacuum energy in curved space we will follow a procedure similar to the above by subtracting the flat space vacuum energy. We will find that in the _Friedmann-Robertson-Walker_ spacetime the vacuum energy of a scalar field is the same as in flat space. Then we expect the renormalised cosmological constant to vanish. However we do observe a tiny cosmological constant
$
  Lambda/G_N tilde M_p^2 H_0^2 tilde 10^(-120) M_p^4
$
This is much smaller than expected with a cutoff at the Planck scale $tilde M_p^4$ and non-zero. We might expect this vacuum energy has some dynamical origin which we refer to as _dark energy_.

Removing the UV divergence as above can be interpreted physically as not knowing the theory above the UV cutoff. We remove some negative contribution to the vacuum energy given by the unkown physics above the cutoff. This is a problem since the cosmological constant becomes very sensitive to these unknown UV physics. We can write
$
  Lambda_"renorm" = Lambda_"bare" + Lambda_"counter-term"
$
where $Lambda_"bare" tilde M_p^4$ and $Lambda_"counter-term"$ is some counter term which depends on the UV physics. This counter term is supposed to cancel $Lambda_"bare"$ with a precision on the order of $10^120$. Then any small change in the UV physics such as a phase transition, symmetry breaking etc. will change the IR cosmological constant enormously! This _fine-tuning_ seems absurd and is typically referred to as the _cosmological constant problem_.

As a concrete example consider $ Lambda_"renorm"^"TeV" = Lambda_"bare"^"TeV" + Lambda_"CT"^"TeV" $

When $tilde 100 "GeV"$ we have the electroweak phase transition which changes the vacuum energy by some significant $dd(V, d: Delta)$ meaning $ Lambda_"obs" -> Lambda_"obs"^"TeV" + dd(V, d: Delta) $
This breaks our previous fine-tuning. We could instead consider $ Lambda_"renorm"^"meV" = Lambda_"bare"^"meV" + Lambda_"CT"^"meV" $ When $tilde "meV"$ all _fundamental_ phase transitions have occured so nothing in principal should break our fine-tuning. However now $Lambda_"bare"^"meV"$ is implied to know everything about what has happened at larger scales while also being very sensitive to what happened. This is quite non-sensical and runs in to many problems. We strongly believe that scales do decouple. The minute details of e.g. quantum gravity should not change the dynamics of relativity at very large scales.
