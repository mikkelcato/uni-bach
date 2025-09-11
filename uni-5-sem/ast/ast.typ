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

Cosmology deals with very large and very small scales. For large scales we generally use the megaparsec and the Sun as a reference---like we're used to. For very small scales we can use the fundamental constants $G, c, k_B "and" hbar$ to construct the Planck scale:
$
  l_P & equiv ((G hbar)/c^3)^(1\/2) = 1.62 times 10^(-35) "m" \
  M_P & equiv ((hbar c)/G)^(1\/2) = 2.18 times 10^(-8) "kg" \
  t_P & equiv ((G hbar)/c^5)^(1\/2) = 5.39 times 10^(- 44) "s" \
  E_P & equiv M_P c^2 = 1.22 times 10^28 "eV" \
  T_P & equiv E_P \/k_B = 1.42 times 10^32 "K"
$
this is also known as the natural scale---since if distance, mass, time, and temperature are measured in Planck units then $c = k_B = hbar = G = 1$ by definition.
#pagebreak()
== Fundamental observations
=== Why is the night sky dark?
This is known as Olbers' paradox. To see the problem we compute the brightness of the sky given an infinite universe.

Let $n_*$ be the number density of stars in the universe, and $R_*$ the typical size of a star. Now consider a cylinder of radius $R_*$ about our line of sight, our view will be blocked if the center of a star lies within this cylinder. If the length of the cylinder is $lambda$ then it has a volume $V = lambda pi R_*^2$. The average number of stars blocking our view can be found $ N = n_* V = n_* lambda pi R_*^2 $ the typical distance before our view is blocked is the $lambda$ for which $N = 1$ implying $ lambda = 1/(n_* pi R_*^2) < oo $

so the sky should be full of stars. To find the brightness consider a star of size $R_*$ at a distance $lambda >> R_*$. It takes up an angular area of $Omega = pi R_*^2 \/ lambda^2$ on the sky. If the luminosity is $L_*$ then the flux measured is $f = L_* \/ 4 pi lambda^2$. Then the brightness is $ Sigma_* = f/Omega = L_* /(4 pi^2 R_*^2) $ this does not depend on $lambda$, so the surface brightness of the sky will be equal to the surface brightness of any star---so as bright as the Sun. So what is wrong?

First we assumed that space is transparent over very long distances. However, this wouldn't really change much, since whatever is absorbing the light would be heated. Secondly we assumed that the universe is infinitely large. If the universe extends to some $r_"max" << lambda$, then only a fraction of the night sky would be covered with stars. Similarly if the universe is infinite, but contains no stars beyond some $r_"max"$. Thirdly we assumed that the universe is infinitely old. The speed of light is finite, so looking further out in space corresponds to looking back in time. If the universe has some age $t_0 << lambda\/c$, then we can't see beyond $r tilde c t_0$ and only a fraction of the sky would be covered with stars. Similarly if the universe has only contained stars for a finite time $t_0$. Fourthly we assumed that the surface brightness of a star is independent of distance. This is only true if the universe is static. But in an expanding unviverse the brightness of distance sources would be decreased, and vice versa for a contracting universe.

=== The cosmological principle
For large scales the universe is isotropic and homogeneous, so there is no preferred direction or location in the universe---large scales correspond to $tilde 100 "Mpc"$. This is not obvious, as an example consider a $3 "Mpc"$ sphere centered on our location, within this $tilde 90%$ of the luminosity would be within the Milky Way and M31, so we can define a preferred direction. Similarly the universe is lumpy or inhomogeneous at small scales, even the $3 "Mpc"$ sphere is denser by an order of magnitude compared to the entire universe.

Isotropicity does not imply homogeneity and vice versa. However, we adopt the Copernican principle---there is nothing special about our location in the universe. The universe around us appears isotropic, so by the Copernican principle it is isotropic everywhere---from this homogeneity follows. The statement that the universe is isotropic and homogeneous on large scales is the cosmological principle.

=== Hubble-Lemaître's law
Looking at emission or absorption lines we observe a redshift
$
  z equiv (lambda_"ob" - lambda_"em")/(lambda_"em")
$
we say a galaxy has redshift $z$.

Treating $z$ as being a result of the Doppler effect we can say $z = v \/c$. Since most galaxies have $z > 0$, this implies that most are moving away from us. Lemaître was the first to say this would be explained by an expanding universe. Hubble followed this by plotting redshift against distance and fitted with
$
  z = H_0/c r
$
which is Hubble's law, intepreting these as Doppler shifts we get
$
  v = H_0 r
$
Hubble initially overestimated that $H_0 approx 500 "km""s"^(-1)"Mpc"^(-1)$, now it seems to be much smaller with a value around $H_0 = 68 plus.minus 2 "km""s"^(-1)"Mpc"^(-1)$. Hubble's law shows us that the universe is expanding homogeneously and isotropicly. To see this consider three galaxies forming a triangle with $r_12 equiv abs(arrow(r)_1 - arrow(r)_2)$ etc. We require the shape of this triangle to be preserved, for this we need $r_12 (t) = a(t) r_12 (t_0)$, and similarly for $r_23$ and $r_13$. We call $a(t)$ the scale factor, and we define $a(t_0) equiv 1$. The velocities take the form $ v_12 (t) = dot(a) r_12 (t_0) = dot(a)/a r_12 (t) $ this relationship is the same for all three galaxies. So in any universe experiencing a homogeneous and isotropic expansion we find a distance-velocity relationship of the form $v = H r$ with $H = dot(a)\/a$.

If everything is moving away from everything, then at some point in the past everything was closer. The time since two galaxies seperated by $r$ were in contact---the Hubble time---is
$ t_0 = r/v = H_0^(-1) = 14.38 plus.minus 0.42 "Gyr" $
assuming no forces acted to change their relative motion. So the observation of redshifts naturally leads to a Big Bang model, i.e. a model where the universe expands from an initially very dense state. We can say the age of the universe is $t_0 tilde H_0^(-1)$---since this is obviously not exact---and the first light can only have travelled a distance $d tilde c t_0 tilde c\/H_0$. This therefore resolves Olbers' paradox, since light from past this wouldn't have had time to reach us.

#pagebreak()
== Cosmic Dynamics
=== Distance
We assume basic knowledge of (special) relativity and skip most of this. A useful metric in cosmology is the Robertson-Walker metric
$
  dd(s)^2 = -c^2 dd(t)^2 + a(t)^2 (dd(r)^2 + S_kappa (r)^2 dd(Omega)^2) "with" S_kappa (r) = cases(R_0 sin r\/R_0 #h(20pt)&(kappa=1), r& (kappa = 0), R_0 sinh r\/R_0& (kappa =-1))
$
which describes a spatially homogeneous and isotropic universe, wherein distances can expand or contract.

To find the proper distance to some thing we need the length of the spatial geodesic when the scale factor is fixed at some $a(t)$. Along the geodesic $(theta, phi)$ is constant, so at a fixed time $t$,
$
  dd(s) = a(t) dd(r) => d_p (t) = a(t) integral_0^r dd(r) = a(t) r
$
giving
$
  dot(d)_p = dot(a)/a d_p => v_p (t_0) = H_0 d_p (t_0)
$
with $H_0 equiv (dot(a)\/a)_(t=t_0)$ and $dot(d)_p equiv v_p$. This is what we had before, but now the change in distance is interpreted as being related to the expansion of space. If $ d_p > d_H (t_0) equiv c \/H_0 => v_p = dot(d)_p >c $ so objects further away than $d_H (t_0)$ move away from us faster than the speed of light.

We can't measure $d_p (t_0)$ easily, but we can measure the redshift $z$. Now consider some light emitted at $t_e$ and observed by us at $t_0$. The light travels with $dd(s)=0$ and constant $(theta, phi)$ so
$
  c^2 dd(t)^2 = a(t)^2 dd(r)^2 => c dd(t)/a(t) = dd(r)
$
now consider the wave crest emitted at $t_e$ and observed at $t_0$ and the wave crest emitted at $t_e + lambda_e\/c$ and observed at $t_0 + lambda_0 \/c$,
$
  c integral_(t_e)^t_0 dd(t)/a(t) = r = c integral_(t_e+lambda_e \/c)^(t_0 + lambda_0\/c) dd(t)/a(t)
$
now
$
  integral_(t_e)^(t_0) dd(t)/a(t) - integral_(t_e + lambda_e\/c)^(t_0) dd(t)/a(t) &= integral_(t_e + lambda_e\/c)^(t_0 + lambda_0 \/c) dd(t)/a(t) - integral_(t_e + lambda_e\/c)^t_0 dd(t)/a(t) \
  integral_(t_e)^(t_e + lambda_e\/c) dd(t)/a(t) &= integral_(t_0)^(t_0 + lambda_0 \/c) dd(t) /a(t)
$
taking the time between successive crests to be very small $a(t)$ is approximately constant
$
  1/a(t_e) integral_e = 1/a(t_0) integral_0 => lambda_e /(a(t_e)) = lambda_0/a(t_0) => 1 + z = a(t_0)/a(t_e) = 1/a(t_e)
$
which is nice.

=== Friedmann equation
We seek equations linking $a(t), kappa$ and $R_0$---the geometry---with $epsilon(t)$ and $P(t)$---the contents of the universe.

The first four are related by the Friedmann equation---see GR notes for proper derivation. This course doesn't use relativity proper so we'll derive a Newtonian version and just state the correct version.

We start with Newtons law for a mass $m$ on a sphere of mass $M_s$ which can contract or expand:
$
  F = - (G M_s m)/(R_s (t)^2) = m dv(R_s, t, 2) => 1/2 (dv(R_s, t))^2 = (G M_s)/R_s + U
$
with $U$ being a constant of integration. Roughly this states that $epsilon_"kin" + epsilon_"pot" = "constant"$. We can write
$
  M_s = (4 pi)/3 rho(t) R_s (t)^3
$
and let $R_s (t) = a(t) r_s$. Then we obtain
$
  (dot(a)/a)^2 = (8 pi G)/3 rho(t) + (2 U)/r_s^2 1/a(t)^2
$
but this is very wrong and assumes the universe is Euclidean, which it isn't.

The correct (first) Friedmann equation is
$
  (dot(a)/a)^2 = (8 pi G)/(3 c^2) epsilon(t) - (kappa c^2)/R_0^2 1/a(t)^2
$
and it can be rewritten in many ways to incorporate quantities that are nicer to work with. Defining $H(t) equiv dot(a)\/a$ we get
$
  H(t)^2 = (8 pi G)/(3 c^2) epsilon(t) - (kappa c^2)/R_0^2 1/(a(t)^2) => H_0^2 = (8 pi G)/(3 c^2) epsilon_0 - (kappa c^2)/R_0^2
$
in the case of a spatially flat universe $kappa = 0$ then
$
  H(t)^2 = (8 pi G)/(3 c^2) epsilon(t) =>^"critical density" epsilon_c (t) equiv (3 c^2)/(8 pi G) H(t)^2
$
if $epsilon > epsilon_c$ then $kappa = +1$ and if $epsilon < epsilon_c$ then $kappa = -1$---sometimes we use the equivalent mass density $rho_c equiv epsilon_c\/c^2$. We can also define the density parameter
$
  Omega(t) equiv epsilon(t)/(epsilon_c (t))
$
giving
$
  1 - Omega(t) = - (kappa c^2)/(R_0^2 a(t)^2 H(t)^2) => kappa/R_0^2 = H_0^2/c^2 (Omega_0 -1)
$
so if we know $Omega_0$ then we get the sign of $kappa$ and if we know $c\/H_0$ then we can find $R_0$---currently $Omega_0$ seems very close to $1$.

=== Fluid and Acceleration equations
The Friedmann equation is nice, but we still need something relating $a(t)$ and $epsilon(t)$. We start by deriving the fluid equation, again using a Newtonian approach---however, this time it turns out not to matter. We start with the first law of thermodynamics
$
  dd(Q) = dd(E) + P dd(V)
$
if the universe is homogeneous then $dd(Q) = 0$, i.e. there is no bulk flow of heat---no entropy increase $dd(S)=0$. So we have
$
  dot(E) + P dot(V) = 0
$
consider a sphere with $R_s (t) = a(t) r_s$, then
$
  V(t) = (4 pi)/3 r_s^3 a(t)^3 => dot(V) = (4pi)/3 r_s^3 (3 a^2 dot(a)) = V ((3 dot(a))/a)
$
and
$
  E(t) = V(t) epsilon(t) => dot(E) = V dot(epsilon) + dot(V) epsilon = V (dot(epsilon) + (3 dot(a) epsilon)/a)
$
so we obtain
$
  dot(epsilon) + (3dot(a))/a (epsilon + P) = 0
$
we can combine this with the Friedmann equation to get the acceleration equation (sometimes the second Friedmann equation)---we can write the Friedmann equation as
$
  dot(a)^2 &= (8 pi G)/(3 c^2) epsilon a^2 - (kappa c^2)/R_0^2 \
  dv(, t) => 2 dot(a) dot.double(a) &= (8 pi G)/(3 c^2) (dot(epsilon) a^2 + 2 epsilon a dot(a)) \
  dot.double(a)/a &= (4 pi G)/(3c^2) (dot(epsilon) a/dot(a) + 2 epsilon) \
  "fluid" => dot.double(a)/a &= - (4 pi G)/(3 c^2) (epsilon + 3 P)
$
which is the acceleration equation. If $epsilon > 0$ then it leads to a negative acceleration---it decreases the value of $dot(a)$ and reduces the relative velocity of any two points in the universe. It also includes the pressure associated with the content in the universe. Baryonic matter has positive $P$, as does photons---this also leads to a negative acceleration. But if something has $P < - epsilon\/3$ then it would lead to a positive acceleration.

=== Equations of state
What we have so far are three equations, two of which are independent. But we have three unknowns $a(t), epsilon(t) "and" P(t)$, this is a problem. What we need is some function $P equiv P(epsilon)$---an equation of state.

Many of these are very simple for our purposes---it can typically be written in linear form $P = w epsilon$ where $w$ is dimensionless. As an example take a non-relativistic gas which obeys the perfect gas law
$
  P_"nonrel" = rho/mu k_B T tilde.eq^(epsilon approx rho c^2) (k_B T)/(mu c^2) epsilon = expval(v^2)/(3 c^2) epsilon = w epsilon_"nonrel"
$
with $w << 1$. Most gases are non-relativistic, even in cosmology, so this is quite nice. For a gas of photons---or any other relativistic gas---we have
$
  P_"rel" = 1/3 epsilon_"rel"
$
just take $expval(v^2) tilde c^2$ in the previous.

Some $w$ are of course of more interest than others. We know the universe contains non-relativistic matter---so $w = 0$ is of interest---we refer to this component as matter. It also contains photons---so $w = 1\/3$ is of interest---we refer to this component as radiation. A component with $w < -1\/3$ is also of interest, since this would provide a positive acceleration---this is referred to as dark energy---evidence suggests that the universe contains a cosmological constant $Lambda$ with $w = -1 => P = - epsilon$.

Adding a cosmological constant changes our dynamics---which is why Einstein introduced it, he wanted a static universe. The Friedmann equation and acceleration equation change as
$
  (dot(a)/a)^2 &= (8 pi G)/(3 c^2) epsilon - (kappa c^2)/R_0^2 1/a^2 + Lambda/3 \
  dot.double(a)/a &= - (4 pi G)/(3 c^2) (epsilon + 3 P) + Lambda/3
$
from the Friedmann equation we see that adding $Lambda$ is equivalent to adding a component with
$
  epsilon_Lambda equiv c^2/(8 pi G) Lambda
$
if $Lambda$ is constant, then $epsilon_Lambda$ is constant---from the fluid equation this requires that
$
  P_Lambda = - epsilon_Lambda = - c^2/(8 pi G) Lambda
$
so we can treat $Lambda$ as a component. For a static universe with $dot.double(a) = 0$ we require $Lambda = 4 pi G rho$, in this case the universe would have $kappa = +1$ with
$
  R_0 = c/Lambda^(1\/2)
$
however Einstein discarded this after Hubble showed everything is moving apart---i.e. the universe is expanding. Currently it seems like $Lambda$ is needed, but its value should be greater than what Einstein believed, given that the universe is expanding---and the expansion is accelerating.

#pagebreak()
== Model Universes
