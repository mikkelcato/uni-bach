//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures
#import "@preview/mannot:0.3.0": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Semiconductor crystals
Semiconductors are very important materials and lay the groundwork for all modern electronics. Qualitatively we can define semiconductors as materials with electrical resistivity $10^(-2) -> 10^9 "ohm-cm"$, with strong temperature dependence. We define the band gap as the energy difference between the bottom of the conduction band and top of the valence band. At $T = 0"K"$ the conduction band is vacant, while the valence band is filled. As the temperature increases electrons can be thermally excited to the conduction band, leaving holes in the valence band. This is called intrinsic conductivity, since it occurs in pure semiconductors, and is largely proportional to $E_g\/k_B T$---where $E_g$ can be deduced using optical absorption\*.

== Equations of motion
We'd like to know the motion of an electron in an energy band. We consider the motion of a wave packet in an electric field, and we assume this wave packet is concentrated about wavevector $k$. We know $omega = epsilon.alt\/hbar$, so the velocity of the wavepacket is
$
  bold(v) = 1/hbar nabla_bold(k) epsilon.alt (bold(k))
$
the effect of the crystal is contained in the dispersion relation $epsilon.alt(bold(k))$.

We can write the work done on the electron (in one-dimension) in two ways
$
  dd(epsilon.alt, d: delta) &= underbrace(- e E, F) v_g dd(t, d: delta) \
  &= dv(epsilon.alt, k) dd(k, d: delta) = hbar v_g dd(k, d: delta)
$
implying $hbar dot(k) = -e E$. Actually we can write
$
  markrect(hbar dot(bold(k)) = bold(F), outset: #.3em)
$
#proof[
  We consider the Bloch functions $psi_bold(k)$ to $epsilon.alt_bold(k)$ and $bold(k)$:
  $
    psi_bold(k) = sum_bold(G) C(bold(k)+bold(G)) exp[i(bold(k)+bold(G))dot bold(r)]
  $
  then the momentum of an electron in Bloch state $bold(k)$ is
  $
    bold(P) &= braket(bold(k), -i hbar nabla, bold(k)) \
    &= sum_bold(G) hbar (bold(k)+bold(G)) abs(C(bold(k)+bold(G)))^2 \
    &=^(sum abs(dots)^2 = 1) hbar (bold(k) + sum_bold(G) bold(G) abs(C(bold(k)+bold(G)))^2)
  $
  We now consider the momentum transfer between electron and lattice when the electron goes from the state $bold(k)$ to the state $bold(k) + dd(bold(k), d: Delta)$---we imagine we just have a single electron. If the electron were free then
  $
    underbrace(bold(J), "impulse") = dd(bold(p)_"tot", d: Delta) = dd(bold(p)_"el", d: Delta) = hbar dd(bold(k), d: Delta)
  $
  if the electron interact with the periodic potential of the lattice then
  $
    bold(J) = dd(bold(p)_"tot", d: Delta) = dd(bold(p)_"lat", d: Delta)+dd(bold(p)_"el", d: Delta)
  $
  by the above we know
  $
    dd(bold(p)_"el", d: Delta) = hbar dd(bold(k), d: Delta) + sum_bold(G) hbar bold(G) [(nabla_bold(k) abs(C(bold(k)+bold(G)))^2) dot dd(bold(k), d: Delta)]
  $
  for $dd(bold(p)_"lat", d: Delta)$ we consider an electron reflected. If it has initial momentum $hbar bold(k)$ then after reflection it has $hbar (bold(k)+bold(G))$, so the lattice acquires momentum $-bold(G) hbar$. This means that when $psi_bold(k) -> psi_(bold(k)+dd(bold(k), d: Delta))$
  $
    dd(bold(p)_"lat", d: Delta) = - hbar sum_bold(G) bold(G) underbrace([(nabla_bold(k) abs(C(bold(k)+bold(G)))^2) dot dd(bold(k), d: Delta)], "each gets reflected")
  $
  So we obtain
  $
    dd(bold(p)_"tot", d: Delta) = hbar dd(bold(k), d: Delta) => hbar dot(bold(k)) = bold(F)
  $
]


=== Holes
As mentioned in the beginning electrons can be thermally excited, and when this happens they leave behind holes. These vacant orbitals are of great importance, since they act as if they had charge $+e$!

#proof[
  By symmetry it follows that
  $
    sum bold(k) =^! 0 "for filled band"
  $
  since if $bold(k)$ is filled, then by inversion $-bold(k)$ is also filled. This implies
  $
    bold(k)_h = - bold(k)_e
  $
  if the band is symmetric under $bold(k) <-> - bold(k)$ then
  $
    epsilon.alt_e (bold(k)_e) = epsilon.alt_e (- bold(k)_e) = - epsilon.alt_h (-bold(k)_e) = - epsilon.alt_h (bold(k)_h)
  $
  this implies
  $
    bold(v)_h = bold(v)_e
  $
  since $nabla epsilon.alt_h (bold(k)_h) = nabla epsilon.alt_e (bold(k)_e)$. Below we show $m^(-1) prop dv(epsilon.alt, k, 2)$ this implies
  $
    m_h = - m_e
  $
  if $bold(k)_h = - bold(k)_e$ and $bold(v)_h = bold(v)_e$ it immediately follows that holes act like $+e$ by substitution:
  $
    hbar dot(bold(k))_e & = - e (bold(E) + 1/c bold(v)_e times bold(B)) \
    hbar dot(bold(k))_h & = + e (bold(E) + 1/c bold(v)_h times bold(B))
  $
]

=== Effective mass
Recall the free electron dispersion relation
$
  epsilon.alt = hbar^2/(2m) k^2
$
we notice that $m^(-1)$ determines the curvature of $epsilon.alt$. In the previous part we found the dispersion relation for an electron near the band gap:
$
  epsilon.alt (K) &= epsilon.alt_c + hbar^2/(2 m_e) K^2";   " m_e/m = [(2 lambda)/U - 1]^(-1) "upper band"\
  epsilon.alt (K) &= epsilon.alt_v - hbar^2/(2 m_h) K^2";   " m_h/m = [(2 lambda)/U + 1]^(-1) "lower band"
$
we now define the effective mass $m^*$. We can write
$
  dv(v, t) = 1/hbar dv(epsilon.alt, k, 2) underbrace(dv(k, t), F\/hbar)
$
and we obtain
$
  F = underbracket(hbar^2/(dif^2 epsilon.alt\/dif k^2), equiv m^*) dv(v, t)
$
or
$
  markrect(1/m^* = 1/hbar^2 dv(epsilon.alt, k, 2), outset: #.3em)
$
so the electron would respond to forces as if it had mass $m^*$---this happens due to Bragg reflections affecting the momentum.

== Intrinsic carrier concentration
Assuming $epsilon.alt-mu >> k_B T$ we can write the probability that an orbital is occupied as
$
  f_e tilde.eq exp[(mu-epsilon.alt)/(k_B T)]
$
the energy of an electron in the conduction band is
$
  epsilon.alt_k = E_c + (hbar^2)/(2m_e) k^2
$
and
$
  D_e (epsilon.alt) = 1/(2 pi^2) ((2 m_e)/hbar^2)^(3\/2) (epsilon.alt - E_c)^(1\/2)
$
which we use to obtain the concentration of electrons in the conduction band:
$
  n & = integral_(E_c)^oo D_e (epsilon.alt) f_e (epsilon.alt) dd(epsilon.alt) \
    & = 2 ((m_e k_B T)/(2 pi hbar^2))^(3\/2) exp[(mu - E_c)/(k_B T)]
$
We can also compute the concentration of holes $p$. Using $f_h = 1 - f_e$ we find
$
  f_h tilde.eq exp[(epsilon-mu)/(k_B T)]
$
if they behave like particles with $m^* = m_h$ we have
$
  D_h (epsilon.alt) = 1/(2 pi^2) ((2 m_h)/hbar^2)^(3\/2) (E_v - epsilon.alt)^(1\/2)
$
the concentration of holes in the valence band is then
$
  p & = integral_(-oo)^(E_v) D_h (epsilon.alt) f_h (epsilon.alt) dd(epsilon.alt) \
    & = 2 ((m_h k_B T)/(2 pi hbar^2))^(3\/2) exp[(E_v-mu)/(k_B T)]
$
combining the two equations
$
  n p= 4 ((k_B T)/(2 pi hbar^2))^3 (m_c m_h)^(3\/2) exp(- E_g/(k_B T))
$
where $E_g = E_c-E_v$---note this also holds for impure semiconductors, we've only assumed $mu-epsilon.alt >> k_B T$.

If we have an intrinsic semiconductor then $n_i = p_i$:
$
  n_i = p_i = 2 ((k_B T)/(2 pi hbar^2))^(3\/2) (m_e m_h)^(3\/4) exp(- E_g/(2 k_B T))
$
we can also set $n = p$ giving
$
  mu = 1/2 E_g + 3/4 k_B T ln m_h/m_e =^(m_h = m_e) 1/2 E_g
$
as measured from the top of the valence band. So the Fermi level lies exactly in the middle of the forbidden gap. (this can also be derived by requiring $f_e (E_c) = f_h (E_v)$)

We define
$
  sigma equiv (n e mu_e + p e mu_h)
$
where $mu_e$ and $mu_h$ is the mobility $mu equiv abs(v) E^(-1)$. We've previously found the drift velocity as $v = q tau E m^(-1)$ where $tau$ is the collision time, so we obtain:
$
  sigma = e^2 (tau_e/m_e n + tau_h/m_h p)
$

== Impurity conductivity
Impurities change the electrical properties of semiconductors---the process of adding impurities is called doping---these can add electrons (donors or n-type) or take electrons (acceptors or p-type) and leave holes.

For a donor state the extra electron moves in the Coulomb potential $e\/epsilon r$ of the impurity ion. Using the Bohr model one can find ionization energy
$
  E_d = (e^4 m_e)/(2 epsilon^2 hbar^2)
$
when ionized the electron is _donated_ to the conduction band. An acceptor state takes an electron leaving a hole. When an acceptor is ionized a hole is freed---requiring energy $E_a$.


#pagebreak()
= Fermi Surfaces
\*225-239

== Zone schemes

== Fermi surfaces & electron orbits

== Tight-binding model
