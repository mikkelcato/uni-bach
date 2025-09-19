//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *problems in astrophysics and cosmology*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    problems from class
  ],
)

= Lecture 2
== The fluid equation
Derive the fluid equation
$
  dot(epsilon) + 3 dot(a)/a (epsilon + P) = 0
$
from the two Friedmann equations.

We have
$
  dot(a)^2/a^2 = (8 pi G)/(3 c^2) epsilon - (kappa c^2)/R_0^2 1/a^2
$
and
$
  dot.double(a)/a = - (4 pi G)/(3 c^2) (epsilon + 3 P)
$
we take the time derivative of the first
$
  (2 dot.double(a))/a dot(a)/a -2 dot(a)^2/a^2 dot(a)/a = (8 pi G)/(3c^2) dot(epsilon) + 2 (kappa c^2)/R_0^2 1/a^2 dot(a)/a
$
clean it up
$
  (4 pi G)/(3 c^2) dot(epsilon) + dot(a)/a ((kappa c^2)/R_0^2 1/a^2 + dot(a)^2/a^2 dot(a)/a - dot.double(a)/a dot(a)/a) = 0
$
and subbing in both Friedmann equations
$
  0&=(4 pi G)/(3 c^2) dot(epsilon) + dot(a)/a ((kappa c^2)/R_0^2 1/a^2 + dot(a)^2/a^2 - dot.double(a)/a)\
  &= (4 pi G)/(3 c^2) dot(epsilon) + dot(a)/a ((kappa c^2)/R_0^2 1/a^2 + (8 pi G)/(3 c^2) epsilon - (kappa c^2)/R_0^2 1/a^2 + (4 pi G)/(3 c^2) (epsilon + 3 P) ) \
  &= dot(epsilon) + dot(a)/a ( 2 epsilon + (epsilon + 3 P) ) \
  &= dot(epsilon) + 3 dot(a)/a (epsilon + P)
$

#pagebreak()
== Solution of the non-relativistic Friedmann equation
Assume the universe only contains non-relativistic particles with $epsilon = c^2 rho$ and $P=0$, and is flat $kappa = 0$.

We start with the fluid equation
$
  dot(epsilon) + 3 dot(a)/a ( epsilon + P ) = 0
$
which becomes
$
   dv(rho (a(t)), t) + 3 dot(a)/a rho(a) & = 0 \
  dv(rho, a) dot(a) + 3 dot(a) /a rho(a) & = 0
$
this can be solved for $rho(a)$
$
                  dv(rho, a) & =- 3 rho(a)/a \
              dd(rho)/rho(a) & = -3 dd(a)/a \
             ln rho(a)/rho_0 & = - 3 ln(a/a_0) \
  a_0 = 1 => ln rho(a)/rho_0 & = ln a^(-3) \
                      rho(a) & = rho_0/a^3
$
as we'd expect. Now we use the first Friedmann equation
$
  (dot(a)/a)^2 = (8 pi G)/(3 c^2) c^2 rho - (kappa c^2)/(R_0^2 a^2)
$
which becomes
$
  (dot(a)/a)^2 & = (8 pi G)/(3 a^3) rho_0 \
        dot(a) & = sqrt((8 pi G rho_0)/3) 1/a^(1/2)
$
or
$
  a^(1/2) dd(a) & = sqrt((8 pi G rho_0)/3) dd(t) \
        a^(3/2) & = 3/2 sqrt((8 pi G rho_0)/3) t \
              a & = root(3, 6 pi G rho_0) t^(2/3)
$

== Equation of state parameter for non-relativistic matter and radiation, Ryden P4.5

Wavelength increases as $lambda prop a$. The total energy density of a gas of particles can be written $epsilon = n E$, with $n$ being number density and $E$ energy, assuming all particles have same $p$ and $m$ then
$
  E = (m^2 c^4 + p^2 c^2)^(1\/2) = (m^2 c^4 + h^2 c^2\/lambda^2)^(1\/2)
$
find $w(a)$.

We know $P = w epsilon$. $n = N/V$ but $V prop a^3$ so $n prop a^(-3)$. The fluid equation
$
  dot(epsilon) + 3 dot(a)/a (1 + w) epsilon = 0
$
$
  n &= (3N)/(4 pi R_0^3 a^3) => dot(n) = (-9 N dot(a))/(4 pi R_0^3 a^4) = - 3 dot(a)/a n
$
$
  E = (m^2 c^4 + (h^2 c^2)/lambda_0^2 1/a^2)^(1\/2) => dot(E) &= 1/2 (m^2 c^4 + (h^2 c^2)/lambda_0^2 1/a^2)^(-1/2) (-2 (h^2 c^2)/lambda_0^2 1/a^2 dot(a)/a) \
  &= (m^2 c^4 - E^2)/E dot(a)/a
$
$
  dot(epsilon) & = E dot(n) + dot(E) n \
               & = -3 dot(a)/a E n + (m^2 c^4 - E^2)/E n dot(a)/a \
               & = dot(a)/a E n ((m^2 c^4 - E^2)/E^2 - 3)
$
so from the fluid equation
$
  0 & = dot(a)/a epsilon ((m^2 c^4 - E^2)/E^2 -3) + 3 dot(a)/a (1 + w) epsilon \
  0 & = (m^2 c^4 - E^2)/E^2 + 3 w \
  w & = 1/3 (E^2-m^2 c^4)/E^2 = 1/3 [(1 - (m^2 c^4)\/E^2)/1]
$
in the non-relativistic limit $E^2 tilde.eq m^2 c^4$ so $w = 0$, and in the relativistic limit $E^2 >> m^2 c^4$ so $w = 1\/3$.

== Radius of curvature in Einsteins universe, Ryden P4.3
If $rho = 2.7 times 10^(-27) "kg""m"^(-3)$ what is $R_0$ of Einstein's static universe. How long would it take for a photon to circumnavigate such a universe?

In a static universe $a = a_0 = 1$ and we require $dot(a) = dot.double(a)=0$. We have the second Friedmann equation
$
  dot(a)/a = - (4 pi G)/(3 c^2) (epsilon + 3 P) + Lambda/3
$
which becomes
$
  0 = -(4 pi G)/3 rho + Lambda/3 => Lambda = 4 pi G rho
$
We have the first Friedmann equation
$
  (dot(a)/a)^2 = (8 pi G)/(3 c^2) epsilon - (kappa c^2)/R_0^2 1/a^2 + Lambda/3
$
or with our assumptions ($epsilon = c^2 rho$),
$
  0 & = (8 pi G)/3 rho - (kappa c^2)/R_0^2 + Lambda/3 \
  0 & = 4 pi G rho - (kappa c^2)/R_0^2 \
    & => R_0 = (c)/(2 sqrt(pi G rho))
$
since we require $kappa = + 1$, otherwise the second line makes no sense. With the given $rho$,
$
  R_0 tilde.eq 2 times 10^26 "m"
$
so for a photon with $v = c$ to circumnavigate this universe would take
$
  t_"c" = (2 pi R_0)/c = 4.18 times 10^18 "s" tilde.eq 133 "Gyr"
$

== Pertubation of the static universe, Ryden P4.2
We again consider Einsteins static universe, there the attractive force from $rho$ cancels the repulsive force from $Lambda = 4 pi G rho$. Suppose some $rho -> "radiation"$. Will the universe start to expand or contract?

The conversion would preserve $epsilon$. But radiation has an associated pressure with $w = 1\/3$ as opposed to matter with $w = 0$. Thus $P$ would increase and $P > 0$ (since in Einstein's static universe it is filled with only matter). Then
$
  dot.double(a)/a = - (4 pi G)/(3 c^2) (epsilon + 3 P) + Lambda/3
$
since $epsilon$ is preserved this term just cancels as before but now
$
  dot.double(a)/a = - (4 pi G P)/c^2 < 0
$
so $dot.double(a) < 0$ meaning the universe would contract.

#pagebreak()
= Lecture 3

