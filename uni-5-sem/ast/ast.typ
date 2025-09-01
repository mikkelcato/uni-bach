//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *astrophysics and cosmology*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on astrophysics and cosmology taken during the SDU course---follows three part structure; from very big to big to small. Based on Ryden's _Introduction to Cosmology_ and _An Introduction to Modern Astrophysics_ by Carroll and Ostie---supplemented with notes taken during lecture.
  ],
)

= Cosmology
Cosmology is the study of the universe, and given the universe is quite complicated we rely on many simplifications. Ignoring small things like galaxies, allows us to try to answer very fundamental questions about the universe as a whole.

Necessarily cosmology deals with very large scales, so we'll need an appropriate system of measurement. For length we use the megaparsec: $1 "Mpc" = 3.09 times 10^22"m"$. For mass we use Sun's mass: $1 "M"_dot.circle = 1.99 times 10^30"kg"$. For power we use the Sun's luminosity: $1 "L"_dot.circle = 3.83 times 10^26 "W"$. For time we use years, megayears or gigayears.

Cosmology also deals with very small scales, especially early in the history of our universe. Here it's common to use eV instead of J, with $1 "eV" = 1.6 times 10^(-19) "J"$. Using the fundamental constants $G, c, k_B "and" hbar$ we can construct the Planck scale:
$
  l_P & equiv ((G hbar)/c^3)^(1\/2) = 1.62 times 10^(-35) "m" \
  M_P & equiv ((hbar c)/G)^(1\/2) = 2.18 times 10^(-8) "kg" \
  t_P & equiv ((G hbar)/c^5)^(1\/2) = 5.39 times 10^(- 44) "s" \
  E_P & equiv M_P c^2 = 1.22 times 10^28 "eV" \
  T_P & equiv E_P \/k_B = 1.42 times 10^32 "K"
$
This is also known as the natural scale---since if distance, mass, time, and temperature are measured in the corresponding Planck units then $c = k_B = hbar = G = 1$ by definition.

#pagebreak()
== Fundamental observations
Cosmology is based, like many other things, on certain observations.

=== Olbers' paradox
Why is the night sky dark? We compute the expected brightness of the sky given an infinite universe. Let $n_*$ be the number density of stars in the universe, and let the typical size of such a star be $R_*$, we'll assume $R_* tilde R_dot.circle$. Considering a cylinder of radius $R_*$ about our line of sight then our view will be blocked if a stars center lies within this cylinder. If the length of the cylinder is $lambda$ then it has a volume $V = lambda pi R_*^2$. The average number of stars blocking our view is then $N = n_* V = n_* lambda pi R_*^2$. The typical distance before our view is blocked is the $lambda$ for which $N = 1$ implying $lambda^(-1) = n_* pi R_*^2$.

This is always finite, meaning in our infinite universe the sky should be completely full of stars. To estimate the brightness, consider a star of size $R_*$ at a distance $r >> R_*$. Its angular area is $Omega = pi R_*^2 \/ r^2$, if the luminosity is $L_*$ then the flux measured is $f = L_* \/ 4 pi r^2$. Then the surface brightness is $Sigma_* = f\/Omega = L_*\/4 pi^2 R_*^2$. This is independent of distance, so the surface brightness of the sky will be equal to the surface brightness of an individual star. We conclude that in an infinite universe the sky is as bright as the Sun. This is obviously wrong.

First we assumed that space is transparent over very long distances, which could be wrong. However, this wouldn't really change much, since whatever is absorbing the light would be heated.

Secondly we assumed that the universe is infinitely large. If the universe extends to some $r_"max" << lambda$, then only a fraction of the night sky would be covered with stars. Similarly if the universe is infinite, but contains no stars beyond some $r_"max"$.

Thirdly we assumed that the universe is infinitely old. The speed of light is finite, so looking further out in space corresponds to looking back in time. If the universe has some age $t_0 << lambda\/c$, then we can't see beyond $r tilde c t_0$ and only a fraction of the sky would be covered with stars. Similarly if the universe has only contained stars for a finite time $t_0$.

Fourthly we assumed that the surface brightness of a star is independent of distance. This is only true if the universe is static. But in an expanding unviverse the brightness of distance sources would be decreased, and vice versa for a contracting universe.

=== The cosmological principle
For large scales the universe is isotropic and homogeneous, so there is no preferred direction or location in the universe---large scales correspond to roughly $100 "Mpc"$ or more. This is not obvious, as an example consider a $3 "Mpc"$ sphere centered on our location, within this $tilde 90%$ of the luminosity would be within the Milky Way and M31, it is easy then to define a preferred direction---only on the scale of superclusters can the universe be treated as isotropic. Similarly the universe is lumpy or inhomogeneous at small scales, even the $3 "Mpc"$ sphere is denser by an order of magnitude compared to the entire universe.

Isotropicity does not imply homogeneity and vice versa, these things are very different. However, we adopt the Copernican principle---that there is nothing special about our location in the universe on large scales. The universe around us appears isotropic, so by the Copernican principle it is isotropic everywhere, and from this homogeneity necessarily follows. The statement that the universe is isotropic and homogeneous on large scales is the cosmological principle.

=== Hubble-Lemaître
Looking at emission or absorption lines we observe a redshift
$
  z equiv (lambda_"ob" - lambda_"em")/(lambda_"em")
$
we say a galaxy has redshift $z$. If $z < 0$ the galaxy is blueshifted, however most galaxies have $z > 0$.

Treating $z$ as being a result of the Doppler effect we can say $z = v \/c$. Since most galaxies are redshifted this implies that most are moving away from us. Lemaître was the first to say this would be explained by an expanding universe. Hubble followed this by plotting redshift against distance and fitted with
$
  z = H_0/c r
$
which is Hubble's law, intepreting these as Doppler shifts we get
$
  v = H_0 r
$
Hubble initially overestimated that $H_0 approx 500 "km""s"^(-1)"Mpc"^(-1)$, now it seems to be much smaller with a value around $H_0 = 68 plus.minus 2 "km""s"^(-1)"Mpc"^(-1)$. Hubble's law shows us that the universe is expanding homogeneously and isotropicly. To see this consider three galaxies forming a triangle with $r_12 equiv abs(arrow(r)_1 - arrow(r)_2)$ etc. We require the shape of this triangle to be preserved, for this we need $r_12 (t) = a(t) r_12 (t_0)$, and similarly for $r_23$ and $r_13$. We call $a(t)$ the scale factor, and we define $a(t_0) equiv 1$. The velocities take the form $ v_12 (t) = dot(a) r_12 (t_0) = dot(a)/a r_12 (t) $ this relationship is the same for all three galaxies. So in any universe experiencing a homogeneous and isotropic expansion we find a distance-velocity relationship of the form $v = H r$ with $H = dot(a)\/a$.

If everything is moving away from everything, then this implies that at some point in the past they were closer. The time since two galaxies seperated by $v = H_0 r$ were in contact is $ t_0 = r/v = H_0^(-1) = 14.38 plus.minus 0.42 "Gyr" $ assuming no forces acted to change their relative motion, this is known as the Hubble time. So one Hubble time ago everything would have been very close together. So the observation of redshifts naturally leads to a Big Bang model, i.e. a model where the universe expands from an initially very dense state. Naturally this number is not the true age of the universe, since gravity and other forces exist.

We can say the age of the universe is $t_0 tilde H_0^(-1)$, and the first light can only have travelled a distance $d tilde c t_0 tilde c\/H_0$. This therefore resolves Olbers' paradox, since light from past this wouldn't have had time to reach us.

#pagebreak()
= Structure formation

#pagebreak()
= Stars
