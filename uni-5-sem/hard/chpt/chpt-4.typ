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
The Fermi surface is simply the surface of constant $epsilon.alt_F$ in $bold(k)$-space---it acts to seperate unfilled and filled orbitals at $T = 0"K"$. All electrical properties are thus determined by the shape and volume of the Fermi surface, since current is due to changes in the occupancy of states near the Fermi surface.

So we want to know how to construct the Fermi surface, which requires zone schemes.

== Zone schemes
=== Reduced zone scheme
We can always pick $bold(k)$ of any Bloch function to lie within the first Brillouin zone. This procedure is called _mapping_ the band in the _reduced zone scheme_.

If we have some Bloch function $psi_bold(k)' = e^(i bold(k)'dot bold(r)) u_bold(k)'$ with $bold(k)'$ outside the first Brillouin zone, then we can always find $bold(G)$ such that $bold(k) = bold(k)' + bold(G)$ is inside the first Brillouin zone. Then
$
  psi_bold(k)' = e^(i bold(k)' dot bold(r)) u_bold(k)' = e^(i bold(k) dot bold(r)) underbracket((e^(-i bold(G) dot bold(r)) u_bold(k)'), "periodic") = e^(i bold(k) dot bold(r)) u_bold(k) = psi_bold(k)
$
For free electrons we can also use the scheme. In this case any energy $epsilon.alt_bold(k)'$ with $bold(k)'$ outside the first Brillouin zone is equal to $epsilon.alt_bold(k)$ with $bold(k)$ inside the first Brillouin zone with $bold(k) = bold(k)' + bold(G)$. So for each band we only need to solve for the energies within the first Brillouin zone. We may therefore find different energies, corresponding to different bands, for the same $bold(k)$. In fact the wavefunctions with different energies, but same $bold(k)$, will be independent! This is due to the coefficients $C(bold(k)+bold(G))$ differing for different bands, so we should really write the Bloch function as:
$
  psi_(n,bold(k)) = e^(i bold(k) dot bold(r)) u_(n,bold(k)) = sum_bold(G) C_n (bold(k)+bold(G)) e^(i (bold(k)+bold(G)) dot r)
$
for the $n$th band.

=== Periodic zone scheme
Given we can translate any band to the first Brillouin zone from other zones, then we should be able to do the opposite and translate the first Brillouin zone into all other zones. The energy of a band $epsilon.alt_bold(k)$ becomes a periodic function
$
  epsilon.alt_(bold(k)) = epsilon.alt_(bold(k)+bold(G))
$
and these are understood to both refer to the same band---this is the _periodic zone scheme_.

The last common scheme is the _extended zone scheme_, here you essentially just plot $epsilon.alt_bold(k)$.

== Fermi surfaces & electron orbits
=== Construction
We go through an example construction of the Fermi surface for a square lattice:

We start by finding just the first three Brillouin zones. This is done by using $bold(G)_1$, $bold(G)_2$ and $bold(G)_3$ in Figure ...

The free electron Fermi surface is determined by
$
  epsilon.alt_F = hbar^2/(2 m) bold(k)_F^2 => "circle"
$
so the Fermi surface is a circle of radius $k_F$, meaning electrons occupy all states with $k_x^2 + k_y^2 <= k_F^2$. We map this region with the reduced zone scheme as in Figure ...

so parts of the Fermi surface fall in the second, third and fourth Brillouin zones, while the first Brillouin zone is completely occupied.

With the free electron Fermi surface in hand we can consider what happens in the nearly free electron case.\*

=== Orbits
Previously we showed that electrons move according to $hbar dot(bold(k)) = bold(F)$ in the case of a static magnetic field $bold(B)$ this can be written as
$
  dot(bold(k)) = - e/(hbar^2 c) nabla_bold(k) epsilon.alt times bold(B)
$
this tells us that the electron moves in a direction normal to $nabla_bold(k) epsilon.alt$, meaning on a surface of constant energy!---so an electron on the Fermi surface will move in a curve on the Fermi surface. This leads to three types of orbits as in Figure ...

Generally orbits that enclose filled states are electron orbits, and orbits that enclose empty states are hole orbits, while orbits that move from zone to zone without closing are open orbits.


