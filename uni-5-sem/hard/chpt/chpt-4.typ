//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Semiconductor crystals
Semiconductors are very important materials and lay the groundwork for all modern electronics. As we have seen certain materials have band gaps. For semiconductors the band gap between the conduction band and the valence band is $tilde 2 "eV"$. At $T = 0"K"$ the conduction band is vacant while the valence band is filled. As the temperature increases electrons can be thermally excited to the conduction band, leaving holes in the valence band. This is called intrinsic conductivity since it occurs in all semiconductors. All of this follows from the Fermi-Dirac distribution as we will see.

Conductivity can also be induced by absorption of a photon giving an electron energy $+hbar omega$. This can happen directly (direct bandgap) or indirectly (indirect bandgap) with the help of a phonon giving energy $+ hbar K$. Another way to induce conductivity is by adding impurities.

== Equations of motion
We would like to determine the motion of an electron in an energy band. We consider the motion of a wavepacket in an electric field. By assumption this wavepacket is concentrated about some wavevector $k$. We know $omega = hbar^(-1) epsilon.alt$ so the velocity of the wavepacket is
$
  bold(v) = 1/hbar nabla_bold(k) epsilon.alt (bold(k))
$
the effect of the crystal is within $epsilon.alt(bold(k))$. We can then write the work done on the electron (in one-dimension) in two ways
$
  dd(epsilon.alt, d: delta) &= underbrace(- e E, F) v dd(t, d: delta) \
  dd(epsilon.alt, d: delta) &= dv(epsilon.alt, k) dd(k, d: delta) = hbar v dd(k, d: delta)
$
implying $hbar dot(k) = -e E$. This generalizes to
$
  hbar dot(bold(k)) = bold(F)
$
#proof[
  We consider the Bloch functions $psi_bold(k)$ to $epsilon.alt_bold(k)$ and $bold(k)$:
  $
    psi_bold(k) = sum_bold(G) C(bold(k)+bold(G)) e^(i(bold(k)+bold(G))dot bold(r))
  $
  then the momentum of an electron in a Bloch state $bold(k)$ is
  $
    bold(P) &= braket(bold(k), -i hbar nabla, bold(k)) \
    &= sum_bold(G) hbar (bold(k)+bold(G)) abs(C(bold(k)+bold(G)))^2 \
    &=^(sum abs(dots)^2 = 1) hbar (bold(k) + sum_bold(G) bold(G) abs(C(bold(k)+bold(G)))^2)
  $
  We now consider the momentum transfer between an electron and the lattice when the electron goes from the state $bold(k)$ to the state $bold(k) + dd(bold(k), d: Delta)$. Taking the electron to be free then
  $
    underbrace(bold(J), "impulse") = dd(bold(p)_"tot", d: Delta) = dd(bold(p)_"el", d: Delta) = hbar dd(bold(k), d: Delta)
  $
  if the electron interacts with the periodic potential of the lattice then
  $
    bold(J) = dd(bold(p)_"tot", d: Delta) = dd(bold(p)_"lat", d: Delta)+dd(bold(p)_"el", d: Delta)
  $
  by the above we know
  $
    dd(bold(p)_"el", d: Delta) = hbar dd(bold(k), d: Delta) + sum_bold(G) hbar bold(G) [(nabla_bold(k) abs(C(bold(k)+bold(G)))^2) dot dd(bold(k), d: Delta)]
  $
  To find $dd(bold(p)_"lat", d: Delta)$ we consider an electron being reflected. The momentum increase of the electron is $+ bold(G) hbar$. So the lattice acquires momentum $-bold(G) hbar$, meaning when $psi_bold(k) -> psi_(bold(k)+dd(bold(k), d: Delta))$
  $
    dd(bold(p)_"lat", d: Delta) = - hbar sum_bold(G) bold(G) underbrace([(nabla_bold(k) abs(C(bold(k)+bold(G)))^2) dot dd(bold(k), d: Delta)], "each gets reflected")
  $
  So we obtain
  $
    dd(bold(p)_"tot", d: Delta) = hbar dd(bold(k), d: Delta) => hbar dot(bold(k)) = bold(F)
  $
]

=== Holes
As mentioned in the beginning electrons can be thermally excited, when this happens they leave behind vacant orbitals. We call these holes and they act as if they had charge $+e$!

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
  This implies
  $
    bold(v)_h = bold(v)_e
  $
  since $nabla epsilon.alt_h (bold(k)_h) = nabla epsilon.alt_e (bold(k)_e)$. Below we show $m^(-1) prop dv(epsilon.alt, k, 2)$ implying
  $
    m_h = - m_e
  $
  Given $bold(k)_h = - bold(k)_e$ and $bold(v)_h = bold(v)_e$ it immediately follows that holes act like $+e$ by substitution into $hbar dot(k) = bold(F)$
  $
    hbar dot(bold(k))_e & = - e (bold(E) + 1/c bold(v)_e times bold(B)) \
    hbar dot(bold(k))_h & = + e (bold(E) + 1/c bold(v)_h times bold(B))
  $
]

=== Effective mass
Recall the free electron dispersion
$
  epsilon.alt_"free" = hbar^2/(2m) k^2
$
we notice that $m^(-1)$ determines the curvature of $epsilon.alt_"free"$. Previously we found the dispersion for a nearly free electron near the band gap
$
  epsilon.alt (k) &= epsilon.alt_c + hbar^2/(2 m_e) k^2";   " m_e/m = [(2 lambda)/U - 1]^(-1) "upper band"\
  epsilon.alt (k) &= epsilon.alt_v - hbar^2/(2 m_h) k^2";   " m_h/m = [(2 lambda)/U + 1]^(-1) "lower band"
$
We now define the effective mass $m^*$. Consider
$
  dv(v, t) = 1/hbar dv(epsilon.alt, k, 2) underbrace(dv(k, t), F\/hbar)
$
so we obtain
$
  F = underbrace(hbar^2/(dif^2 epsilon.alt\/dif k^2), equiv m^*) dv(v, t)
$
or
$
  1/m^* = 1/hbar^2 dv(epsilon.alt, k, 2)
$
so the electron responds to forces (the crystal potential) as if it had mass $m^*$.

As an example for the tight-binding model we found in the project that
$
  m^* = hbar^2/(2 t a^2)
$
with
$
  dd(epsilon.alt, d: Delta) = 4 t
$

== Carrier conductivity
Assuming $epsilon.alt-mu >> k_B T$ we can write the probability that an orbital is occupied as
$
  f_e tilde.eq exp[(mu-epsilon.alt)/(k_B T)]
$
the energy of an electron in the conduction band is
$
  epsilon.alt_k = E_c + (hbar^2)/(2m_e) k^2
$
with $m_e$ being the effective mass. Recall the density of states in three-dimensions
$
  D_e (epsilon.alt) = 1/(2 pi^2) ((2 m_e)/hbar^2)^(3\/2) (epsilon.alt - E_c)^(1\/2)
$
which we use to obtain the concentration of electrons in the conduction band
$
  n & = integral_(E_c)^oo D_e (epsilon.alt) f_e (epsilon.alt) dd(epsilon.alt) \
    & = 2 ((m_e k_B T)/(2 pi hbar^2))^(3\/2) exp[(mu - E_c)/(k_B T)] \
    & equiv N(T) exp[(mu - E_c)/(k_B T)]
$
We can also compute the concentration of holes $p$. Using $f_h = 1 - f_e$ we find
$
  f_h tilde.eq exp[(epsilon.alt-mu)/(k_B T)]
$
taking them to behave like particles with effective mass $m_h$ we have
$
  D_h (epsilon.alt) = 1/(2 pi^2) ((2 m_h)/hbar^2)^(3\/2) (E_v - epsilon.alt)^(1\/2)
$
The concentration of holes in the valence band is then
$
  p & = integral_(-oo)^(E_v) D_h (epsilon.alt) f_h (epsilon.alt) dd(epsilon.alt) \
    & = 2 ((m_h k_B T)/(2 pi hbar^2))^(3\/2) exp[(E_v-mu)/(k_B T)] \
    & equiv P(T) exp[(E_v - mu)/(k_B T)]
$
combining the two equations gives the law of mass action
$
  n p & = 4 ((k_B T)/(2 pi hbar^2))^3 (m_e m_h)^(3\/2) exp(- E_g/(k_B T)) \
      & = N(T) P(T) exp(- E_g/(k_B T))
$
where $E_g = E_c-E_v$. For an intrinsic semiconductor with $n_i = p_i$ we find
$
  n_i = p_i = 2 ((k_B T)/(2 pi hbar^2))^(3\/2) (m_e m_h)^(3\/4) exp(- E_g/(2 k_B T))
$
We can also set $n = p$ to determine $mu$
$
  mu = 1/2 E_g + 3/4 k_B T ln m_h/m_e =^(m_h = m_e) 1/2 E_g
$
as measured from the top of the valence band.

We can also define the Drude conductivity by
$
  sigma equiv (n e mu_e + p e mu_h)
$
where $mu_e$ and $mu_h$ is the mobility $mu equiv abs(v) E^(-1)$. We previously found the drift velocity as $ v = q tau E m^(-1) $ where $tau$ is the collision time. So we obtain
$
  sigma = e^2 (tau_e/m_e n + tau_h/m_h p)
$

=== Extrinsic conductivity
Adding impurities or doping a semiconductor increases conductivity. This can be accomplished by adding donors or acceptors.

We consider adding donors to our semiconductor. The extra electron from this donor will move under the Coulomb potential due to the donor ion $e (epsilon r)^(-1)$. From the Bohr model one can then determine the ionization energy
$
  E_d = (e^4 m_e)/(2 epsilon^2 hbar^2) tilde 10 "meV"
$
this will be a small number. The Bohr radius can also be determined
$
  a_d =(epsilon hbar^2)/(m_e e^2) tilde^"typical" 50 "Å"
$
This is large and leads to impurities having overlapping orbitals. This overlap is equivalent to the appearance of a new impurity band below the conduction band. And hopping between impurities leads to conductivity. Given $E_d tilde k_B T$ electrons can be thermally excited to the impurity band. This type of doped semiconductor is generally called $n$-type.

Similarly we can consider adding acceptors to our semiconductor. These add holes leading to conductivity and an impurity band given by $E_a$ above the valence band. This type of doped semiconductor is called $p$-type.


== Junctions
We consider combining a $n$- and $p$-type semiconductor, giving a $p n$-junction. For a $n$-type semiconductor $mu_n$ is just below $E_c$ and for a $p$-type semiconductor $mu_p$ is just above $E_v$. When combined and assuming thermal equilibrium we eventually have $mu_n = mu_p = mu$. This equilibrium is established by diffusion of electrons from the $n$-type to the $p$-type. When done the depleted zone is established wherein there are no charge carriers, as a consequence an electric field from the $p$-type to the $n$-type is also established stopping further diffusion. We can treat the electric fields as creating a potential of size $phi$ which an electron from the $n$-type has to overcome if it wants to enter the $p$-type.

We have two types of competing currents. The first is the thermal current $J_(i,t)$ due to thermally excited carriers diffusion across the barrier, this is independent of the external field. The second is the recombination current $J_(i, r)$ due to carriers diffusion through the barrier and recombining. Given there is no external field $J_(i,t) (0)=J_(i,r) (0)$. The recombination current can be computed as
$
  J_(i, r) = J_(i,r) (0) exp((e V)/(k_B T))
$
then the total current is
$
  J_"tot" = (J_(n,t) + J_(p, t)) (exp[(e V)/(k_B T)]-1)
$
This shows clear rectifying behaviour. These devices act as diodes or solar cells.

#pagebreak()
= Fermi Surfaces
The Fermi surface is simply the surface of constant $epsilon.alt_F$ in $bold(k)$-space. This surface seperates unfilled and filled orbitals at $T = 0"K"$. All electrical properties are therefore determined by the shape and volume of the Fermi surface, since current is due to changes in the occupancy of states near the Fermi surface. For this reason we want to know how to construct the Fermi surface. This requires zone schemes.

== Zone schemes
=== Reduced zone scheme
We can always pick $bold(k)$ of any Bloch function to lie within the first Brillouin zone. This procedure is called _mapping_ the band in the _reduced zone scheme_.

Given a Bloch function $psi_bold(k)' = e^(i bold(k)'dot bold(r)) u_bold(k)'$ with $bold(k)'$ outside the first Brillouin zone we can always find a $bold(G)$ such that $bold(k) = bold(k)' + bold(G)$ is inside the first Brillouin zone. Then
$
  psi_bold(k)' = e^(i bold(k)' dot bold(r)) u_bold(k)' = e^(i bold(k) dot bold(r)) underbrace((e^(-i bold(G) dot bold(r)) u_bold(k)'), "periodic") = e^(i bold(k) dot bold(r)) u_bold(k) = psi_bold(k)
$
For free electrons the $epsilon.alt_bold(k)'$ with $bold(k)'$ outside the first Brillouin zone is equal to $epsilon.alt_bold(k)$ with $bold(k) = bold(k)' + bold(G)$ inside the first Brillouin zone. So for each band we only need to solve for the energies inside the first Brillouin zone. We may therefore find different energies, corresponding to different bands, for the same $bold(k)$. In fact the wavefunctions with different energies, but same $bold(k)$, will be independent! This is due to the coefficients $C(bold(k)+bold(G))$ differing for different bands, so we should really write the Bloch function as:
$
  psi_(n,bold(k)) = e^(i bold(k) dot bold(r)) u_(n,bold(k)) = sum_bold(G) C_n (bold(k)+bold(G)) e^(i (bold(k)+bold(G)) dot r)
$
for the $n$th band.

=== Periodic zone scheme
Given we can translate any band to the first Brillouin zone from other zones, then can also do the opposite and translate the first Brillouin zone into all other zones. The energy $epsilon.alt_bold(k)$ becomes a periodic function
$
  epsilon.alt_(bold(k)) = epsilon.alt_(bold(k)+bold(G))
$
These are understood to both refer to the same band. This is the _periodic zone scheme_.

The last common scheme is the _extended zone scheme_, here we essentially just plot $epsilon.alt_bold(k)$.

== Fermi surfaces
We go through an example construction of the Fermi surface for a square lattice:

We start by finding just the first three Brillouin zones. This is done by using $bold(G)_1$, $bold(G)_2$ and $bold(G)_3$ in Figure ...

The free electron Fermi surface is determined by
$
  epsilon.alt_F = hbar^2/(2 m) bold(k)_F^2 => "circle"
$
so the Fermi surface is a circle of radius $k_F$, meaning electrons occupy all states with $k_x^2 + k_y^2 <= k_F^2$. We map this region with the reduced zone scheme as in Figure ...

so parts of the Fermi surface fall in the second, third and fourth Brillouin zones, while the first Brillouin zone is completely occupied.

With the free electron Fermi surface in hand we can consider what happens in the nearly free electron case. Here the Fermi surface is generally distorted slightly.

=== Orbits
Previously we showed that electrons move according to $hbar dot(bold(k)) = bold(F)$. For a static magnetic field $bold(B)$ this can be written as
$
  dot(bold(k)) = - e/(hbar^2 c) nabla_bold(k) epsilon.alt times bold(B)
$
This tells us the electron moves in a direction normal to $nabla_bold(k) epsilon.alt$, meaning on a surface of constant energy! So an electron on the Fermi surface will move in a curve on the Fermi surface. This leads to three types of orbits as in Figure ...

Generally orbits that enclose filled states are electron orbits, and orbits that enclose empty states are hole orbits. Orbits that move from zone to zone without closing are open orbits.


