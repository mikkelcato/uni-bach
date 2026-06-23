#import "../../temp.typ": *
#show: chpt-note.with()

= The Unruh Effect
We have seen how the _natural_ vacuum state can change with time. We will now see how the _natural_ vacuum can also be completely observer dependent, even in Minkowski spacetime.

== Rindler spacetime
We consider two observers in Minkowski spacetime. One is at rest while the other is accelerated by some constant acceleration. We would like to relate the vacuum states of these two observers.

We need a _comoving_ coordinate system for the accelerated observer. The metric describing the rest frame of the accelerated observer is called Rindler spacetime. Assuming the acceleration is along $x$ we have
$
  dd(s^2) =^"Minkowski" - dd(t^2) + dd(x^2)
$
We transform our coordinates as#footnote[Here $(tau, xi)$ are the rest frame coordinates.]
$
  t & = 1/a e^(a xi) sinh a tau \
  x & = 1/a e^(a xi) cosh a tau
$
implying
$
  dd(s^2) = e^(2 a xi) (dd(tau^2) - dd(xi^2))
$
Consider an observer along $x^(mu) (tau)$. Then
$
  u^mu u_mu = -1";  " u^mu = dv(x^mu, tau)
$
and
$
  a^mu a_mu equiv alpha^2";  " a^mu = dv(u^mu, tau)
$
Within the observers rest frame we have $u^mu = (1, 0)$ implying $(u^0)^2 = 1$ so
$
  dv(, tau) (u^mu u_mu) = 0
$
implying $a^0 = 0$ so $alpha$ is the proper acceleration as seen by the observer and $a^mu (tau) = vecrow(0, alpha)$. Then
$
  u^0 = sqrt(1 + (u^1)^2)";  " dv(u^1, tau) = alpha sqrt(1 + (u^1)^2)
$
implying
$
  dv(x, tau) & = sinh alpha tau \
  dv(t, tau) & = cosh alpha tau
$
or
$
  x(tau) & = 1/alpha cosh alpha tau \
  t(tau) & = 1/alpha sinh alpha tau
$
Which we recognise above. Then the metric for a uniformly accelerated observer in Minkowski spacetime with
$
  1/alpha equiv 1/a e^(a xi)
$
is
$
  dd(s^2) = e^(2 a xi) (- dd(tau^2) + dd(xi^2))
$
Assuming $xi = "constant"$ we have
$
  x^2 - t^2 = 1/alpha^2
$
which is a hyperbole. Also, there are now horizons at $x = plus.minus t$. We call these the Rindler horizons.

== Quantum fields in Rindler space
Consider the action for a scalar field with $m = 0$
$
  S[phi.alt] = - 1/2 integral dd(x, 2) sqrt(-g) g^(mu nu) partial_mu phi.alt partial_nu phi.alt
$
which is conformally invariant under
$
  g_(mu nu) -> tilde(g)_(mu nu) = Omega^2 g_(mu nu)
$
We have
$
  S[phi.alt] &=^"rest Minkowski" 1/2 integral dd(t, x) [(partial_t phi.alt)^2 - (partial_x phi.alt)^2] \
  S[phi.alt] &=^"comoving Rindler" 1/2 integral dd(tau, xi) [(partial_tau phi.alt)^2 - (partial_xi phi.alt)^2]
$
so
$
  pdv(phi.alt, t, 2) - pdv(phi.alt, x, 2) = 0";  " tilde(U)_k = 1/sqrt(2 k) e^(i k x - i omega t) \
  pdv(phi.alt, tau, 2) - pdv(phi.alt, xi, 2) = 0";  " U_k = 1/sqrt(2 k) e^(i k xi - i omega tau)
$
We obtain two mode expansions
$
  phi.alt(t, x) &=^"rest Minkowski" integral dd(k)/sqrt(2 pi) 1/sqrt(2 abs(k)) [a_k e^(-i abs(k)t+i k x) + a_k^dagger e^(i abs(k)t-i k x)] \
  phi.alt(tau, xi) &=^"comoving Rindler" integral dd(k)/sqrt(2 pi) 1/sqrt(2 abs(k)) [b_k e^(-i abs(k) tau+i k xi) + b_k^dagger e^(i abs(k)tau-i k xi)]
$
implying
$
  a_k ket(0_"M") = 0";  " b_k ket(0_"R") = 0
$
with $ket(0_"M") eq.not ket(0_"R")$.

We now define _light-cone coordinates_
$
  overline(u) = t-x";  " overline(v) = t + x \
  u = tau- xi";  " v = t + xi
$
Then
$
  dd(s^2) = - dd(overline(u), overline(v)) = - e^(a (v-u)) dd(u, v)
$
We also have
$
  pdv(phi.alt(overline(u), overline(v)), overline(u), overline(v)) = 0";  " pdv(phi.alt(u, v), u, v) = 0
$
implying
$
  phi.alt(overline(u), overline(v)) & = A(overline(u)) + B(overline(v)) \
                      phi.alt(u, v) & = P(u) + Q(v)
$
We would like these functions. Using $omega = abs(k)$ for the Minkowski fields we have
$
  phi.alt(overline(u), overline(v)) = integral_0^oo dd(omega)/sqrt(2 pi) 1/sqrt(2 omega) [underbracket(a_omega e^(-i omega overline(u)) + a_omega^dagger e^(i omega overline(u)), A(overline(u))) + overbracket(a_(-omega) e^(-i omega overline(v)) + a_(-omega)^dagger e^(i omega overline(v)), B(overline(v)))]
$
and taking $Omega = abs(k)$ for the Rindler fields we have
$
  phi.alt(u, v) = integral_0^oo dd(Omega)/sqrt(2 pi) 1/sqrt(2 Omega) [underbracket(a_Omega e^(-i Omega u) + a_Omega^dagger e^(i Omega u), P(u)) + overbracket(a_(-Omega) e^(-i Omega v) + a_(-Omega)^dagger e^(i Omega v), Q(v))]
$
Consider
$
  phi.alt(u, v) = A(overline(u)(u)) + B(overline(v)(v))
$
implying
$
  A(overline(u)(u)) & = P(u)
$
or
$
  integral dd(omega)/(sqrt(2 pi)) 1/sqrt(2 omega) [a_omega e^(-i omega overline(u)) + a_omega^dagger e^(i omega overline(u))] = integral dd(Omega)/sqrt(2 pi) 1/sqrt(2 Omega) [b_Omega e^(-i Omega u) + b_Omega^dagger e^(i Omega u)]
$
We Fourier transform in $u$. The RHS becomes
$
  integral dd(u)/sqrt(2 pi) e^(i Omega u) P(u) = 1/sqrt(2 abs(Omega)) cases(b_Omega",  "& Omega > 0, b_abs(Omega)^dagger",  "& Omega < 0)
$
and the LHS becomes
$
  integral dd(u)/sqrt(2 pi) e^(i Omega u) A(overline(u)) &= integral_0^oo dd(omega)/sqrt(2 omega) [a_omega F(omega,Omega) + a_omega^dagger F(-omega, Omega)]
$
with
$
  F(omega,Omega) = integral dd(u)/(2 pi) e^(i Omega u - i omega overline(u)) = integral dd(u)/(2 pi) exp[i Omega u + (i omega)/a e^(-a u)]
$
We find by inspection
$
  b_Omega = integral_0^oo dd(omega) [alpha_(omega Omega) a_omega + beta_(omega Omega) a_omega^dagger]
$
with
$
  alpha_(omega Omega) &= sqrt(Omega/omega) F(omega,Omega)";  " beta_(omega Omega) = sqrt(Omega/omega) F(-omega,Omega)
$
Also, we have
$
  [a_omega, a_(omega')^dagger] = delta(omega-omega')";  " [b_Omega,b_(Omega')^dagger] = delta(Omega-Omega')
$
which implies
$
  integral dd(omega) [alpha_(omega Omega) alpha_(omega Omega')^* - beta_(omega Omega) beta_(omega Omega')^*] = delta(Omega-Omega')
$
Then
$
  expval(N_Omega) &= braket(0_"M", b_Omega^dagger b_Omega, 0_"M") \ &= integral dd(omega, omega') braket(0_"M", [underbracket(alpha_(omega Omega)^* a_omega^dagger, 0) + beta_(omega Omega)^* a_omega] [underbracket(alpha_(omega' Omega) a_(omega'), 0) + beta_(omega' Omega) a_(omega')^dagger], 0_"M") \
  &= integral dd(omega, omega') braket(0_"M", beta_(omega Omega)^* beta_(omega' Omega) a_omega a_(omega')^dagger, 0_"M") \
  &= integral dd(omega) abs(beta_(omega Omega))^2 \
  &= n_Omega delta(0)
$

== The Unruh temperature
We had
$
  beta_(omega Omega) = sqrt(Omega/omega) F(-omega, Omega)
$
Consider
$
  F(omega,Omega) = integral dd(u)/(2 pi) exp[i Omega u + (i omega)/a e^(-a u)]
$
We substitute
$
  t = - (i omega)/a e^(-a u)";  " dd(t) = -a t dd(u)
$

Then
$
  F(omega,Omega) &= - 1/(2 pi a) integral dd(t)/t exp[(i Omega)/a ln (omega/(i a t)) -t ] \
  &= -1/(2 pi a) integral dd(t)/t exp[(i Omega)/a (ln omega/a - (i pi)/2 - ln t)- t] \
  &= -1/(2 pi a) exp[(i Omega)/a ln omega/a + (pi Omega)/(2 a)] integral dd(t) t^(-1-i Omega\/a) e^(-t) \
  &= underbracket(-, "contour-dependent")1/(2 pi a) exp[(i Omega)/a ln omega/a + (pi Omega)/(2 a)] Gamma(-(i Omega)/a)
$
We have
$
  ln(-omega) = underbracket(ln(-1), i pi) + ln omega
$
implying
$
  F(omega,Omega) = F(-omega,Omega) e^(pi Omega\/a)
$
Then
$
  delta(Omega-Omega') &= integral_0^oo dd(omega)/omega sqrt(Omega Omega') [F(omega,Omega) F^* (omega,Omega') - F(-omega,Omega) F^* (-omega,Omega')] \
  &= {exp[(pi Omega+pi Omega')/a]-1} integral_0^oo dd(omega)/omega sqrt(Omega Omega') F (-omega,Omega) F^*(-omega,Omega')
$
implying
$
  integral_0^oo dd(omega)/omega sqrt(Omega Omega') F(-omega,Omega) F^* (-omega, Omega') = {exp[(2 pi Omega)/a]-1}^(-1) delta(Omega-Omega')
$
Then
$
  expval(N_Omega) & = integral_0^oo dd(omega) abs(beta_(omega Omega))^2 \
                  & = integral_0^oo dd(omega)/omega Omega abs(F(-omega,Omega))^2 \
                  & = {exp[(2 pi Omega)/a]-1}^(-1) delta(0)
$
so
$
  n_Omega = {exp[(2 pi Omega)/a]-1}^(-1)
$
which is the Bose-Einstein distribution!

We compare this distribution with the thermal spectrum of a $m=0$ particle with energy $E = Omega$
$
  n(E) = {exp(E/T)-1}^(-1)
$
where $T$ is the _Unruh temperature_
$
  T_"Unruh" eq a/(2 pi)
$
Then a uniformly accelerated observer in Minkowski spacetime experiences a thermal spectrum with the temperature $T_"Unruh"$. This is called the _Unruh effect_.

With all units we have
$
  k_B T_"Unruh" = (a hbar)/(2 pi c) tilde "insanely small"
$
implying the Unruh effect is very inefficient. However, it can becomes strong in certain domains.

#pagebreak()
= Hawking Radiation
We know staying outside a black hole requires an acceleration. Then any observer which is not in free-fall experiences the emission of radiation. This is referred to as _Hawking radiation_.

== Rindler-esque derivation of Hawking radiation
Consider the Schwarzschild metric
$
  dd(s^2) = (1 - (2 M G)/r) dd(t^2) - (1- (2 M G)/r)^(-1) dd(r^2) - r^2 dd(Omega^2)
$
We see the metric diverges at the Schwarzschild radius $r_*$ defined by
$
  r_* equiv 2 M G
$
However, the curvature is finite at $r_*$ so this is not a _true_ singularity. We refer to the surface $r=0$ as the true singularity. Any observer falling into the black hole would experience nothing special at $r_*$. But, to a distant observer the infalling observer would _freeze_ at the horizon, eventually disappearing due to redshift. Also, as $r -> r_*$ we find $c_"apparent" -> 0$. This implies nothing can escape upon crossing the horizon. This is why we call them _black_.

We define $rho$ by
$
  rho = integral_(r_*)^r sqrt(g_(r r) (r')) dd(r')
$
Then
$
  dd(s^2) = (1- (2 M G)/(r(rho))) dd(t^2) - dd(rho^2) - r(rho)^2 dd(Omega^2)
$
Close to the horizon we have
$
  rho tilde.eq sqrt(r(r-2 M G)) + 2 M G sqrt(r/(2 M G)-1)
$
implying
$
  dd(s^2) tilde.eq rho^2 (dd(t)/(4 M G))^2 - dd(rho^2) - r^2 (rho) dd(Omega^2)
$
We define $xi$ by
$
  rho = 1/a e^(a xi)
$
Then
$
  dd(s^2) = e^(2 a xi) (dd(t^2) - dd(xi^2))
$
with
$
  a = 1/(4 M G)
$
We recognise this as the Rindler metric for an accelerated observer with proper acceleration
$
  alpha = a e^(-a xi)
$
With $alpha = a$ we have $xi = 0$.

We then immediately find the _Hawking temperature_ by
$
  T_H & = a/(2 pi) \
      & = 1/(8 pi M G)
$
implying black holes evaporate over time.

Note, we refer to the vacuum for the free-falling observer as the _Kruskal vacuum_ $ket(0_"K")$ and the vacuum for an observer at constant $r$ as the _Boulware vacuum_ $ket(0_"B")$. These are analogous to the Minkowski vacuum $ket(0_"M")$ and the Rindler vacuum $ket(0_"R")$.

== in de Sitter spacetime
We can write the de Sitter metric as
$
  dd(s^2) = (1-r^2/r_dd(S)^2) dd(t^2) - (1- r^2/r_dd(S)^2)^(-1) dd(r^2) - r^2 dd(Omega^2)
$
where
$
  r_dd(S) equiv H^(-1)
$
is the de Sitter radius. We can as before determine
$
  T_dd(S) = H/(2 pi) tilde dd(phi.alt_k, d: delta)
$

== Thermodynamics of black holes
We consider black holes to be spherical bodies with radius $r_*$ and surface temperature $T_H$. We can then write the luminosity of a black hole using the Stefan-Boltzmann law
$
  L = gamma sigma T_H^4 A
$
Then
$
  dv(M, t) =-L
$
implying
$
  M(t) = M_0 (1-t/t_L)^(1\/3)
$
with
$
  t_L =(5120 pi M_0^3)/gamma tilde "lifetime"
$
where $M_0$ is the initial black hole mass. We see $t_L tilde M_0^3$ so black holes have to be quite small for this effect to be noticable.

We can determine the entropy of a black hole using
$
  dd(E) = dd(M) = T_H dd(S_"BH") tilde 1^"st" "law"
$
Assuming $S_"BH" prop A$ as Bekenstein did we write
$
  S_"BH" = alpha 16 pi M^2
$
and find
$
  alpha =^! 1/4
$
Which gives the Bekenstein-Hawking entropy
$
  S_"BH" = A/4 tilde (k_B A)/(4 cal(l)_p^2)
$
which is weird since we normally consider $S = ln N$ with $N$ being the number of indistinguishable microstates.

We also find the black hole equation of state
$
  E = 1/(8 pi T_H)
$
implying
$
  C_"BH" = pdv(E, T) = - 1/(8 pi T^2) < 0
$
So black holes become colder upon absorbing heat!

The _area theorem_ states
$
  dv(A, t) >= 0
$
in any _classical_ interaction. Where we note Hawking radiation is a quantum process. This is related to the _generalised $2^"nd"$ law of thermodynamics_
$
  dd(S_"total") = dd(S_"matter") + dd(S_"BH") >= 0
$
which is just the $2^"nd"$ law of thermodynamics where we include $S_"BH"$.

== Consequences
Consider collapsing ordinary matter to form a black hole. Before collapse we have $S prop V$ while after collapse we have $S_"BH" prop A$. This implies a loss of information and a violation of the GLS given $V$ is large enough. This is referred to as the _information paradox_ and is usually resolved by either sacrifing unitarity or locality.

This also leads to the notion of _holography_ which sacrifices locality and postulates
$
  S_"matter" < A/4
$
i.e. the entropy of matter is bounded by the area which encloses it. Then as the matter is collapsed to form a black hole we have no issues.
