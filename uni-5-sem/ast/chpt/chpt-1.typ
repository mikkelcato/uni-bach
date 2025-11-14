//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

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
= Fundamental observations
== Why is the night sky dark?
This is known as Olbers' paradox. To see the problem we compute the brightness of the sky given an infinite universe.

Let $n_*$ be the number density of stars in the universe, and $R_*$ the typical size of a star. Now consider a cylinder of radius $R_*$ about our line of sight, our view will be blocked if the center of a star lies within this cylinder. If the length of the cylinder is $lambda$ then it has a volume $V = lambda pi R_*^2$. The average number of stars blocking our view can be found $ N = n_* V = n_* lambda pi R_*^2 $ the typical distance before our view is blocked is the $lambda$ for which $N = 1$ implying $ lambda = 1/(n_* pi R_*^2) < oo $

so the sky should be full of stars. To find the brightness consider a star of size $R_*$ at a distance $lambda >> R_*$. It takes up an angular area of $Omega = pi R_*^2 \/ lambda^2$ on the sky. If the luminosity is $L_*$ then the flux measured is $f = L_* \/ 4 pi lambda^2$. Then the brightness is $ Sigma_* = f/Omega = L_* /(4 pi^2 R_*^2) $ this does not depend on $lambda$, so the surface brightness of the sky will be equal to the surface brightness of any star---so as bright as the Sun. So what is wrong?

First we assumed that space is transparent over very long distances. However, this wouldn't really change much, since whatever is absorbing the light would be heated. Secondly we assumed that the universe is infinitely large. If the universe extends to some $r_"max" << lambda$, then only a fraction of the night sky would be covered with stars. Similarly if the universe is infinite, but contains no stars beyond some $r_"max"$. Thirdly we assumed that the universe is infinitely old. The speed of light is finite, so looking further out in space corresponds to looking back in time. If the universe has some age $t_0 << lambda\/c$, then we can't see beyond $r tilde c t_0$ and only a fraction of the sky would be covered with stars. Similarly if the universe has only contained stars for a finite time $t_0$. Fourthly we assumed that the surface brightness of a star is independent of distance. This is only true if the universe is static. But in an expanding unviverse the brightness of distance sources would be decreased, and vice versa for a contracting universe.

=== The cosmological principle
For large scales the universe is isotropic and homogeneous, so there is no preferred direction or location in the universe---large scales correspond to $tilde 100 "Mpc"$. This is not obvious, as an example consider a $3 "Mpc"$ sphere centered on our location, within this $tilde 90%$ of the luminosity would be within the Milky Way and M31, so we can define a preferred direction. Similarly the universe is lumpy or inhomogeneous at small scales, even the $3 "Mpc"$ sphere is denser by an order of magnitude compared to the entire universe.

Isotropicity does not imply homogeneity and vice versa. However, we adopt the Copernican principle---there is nothing special about our location in the universe. The universe around us appears isotropic, so by the Copernican principle it is isotropic everywhere---from this homogeneity follows. The statement that the universe is isotropic and homogeneous on large scales is the cosmological principle.

== Hubble-Lemaître's law
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
= Cosmic Dynamics
== Distance
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

== Friedmann equation
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

== Fluid and Acceleration equations
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

== Equations of state
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
= Model Universes
In a homogeneous and isotropic universe we can use
$
  (dot(a)/a)^2 & = (8 pi G)/(3 c^2) epsilon - (kappa c^2)/(R_0^2 a^2) \
             0 & = dot(epsilon) + 3 dot(a)/a (epsilon + P) \
             P & = w epsilon
$
to relate $epsilon(t)$, $P(t)$ and $a(t)$.

== Energy density
The energy density due to different components with different $w_i$ and $epsilon_i$ is additive
$
  epsilon = sum_i epsilon_i
$
similarly for the pressure
$
  P = sum_i w_i epsilon_i
$
assuming there is no interaction between different components then every component should obey a fluid equation
$
  dot(epsilon)_i +3 dot(a)/a (1 + w_i) epsilon_i = 0 => dd(epsilon_i)/epsilon_i = - 3 ( 1 + w_i ) dd(a)/a
$
assuming $w_i$ is constant we obtain
$
  epsilon_i (a) = epsilon_(i,0) a^(-3 (1 + w_i))
$
this tells us how $epsilon_i$ evolves---note that we never used the Friedmann equation. Take $epsilon_m$ with $w = 0$ then
$
  epsilon_m (a) = epsilon_(m,0)/a^3
$
or $epsilon_r$ with $w = 1\/3$ then
$
  epsilon_r (a) = epsilon_(r,0)/a^4
$
the difference is caused by
$
  epsilon_m & = n E = n (m c^2) prop a^(-3) \
  epsilon_r & = n E = n (h c\/lambda) prop a^(-4)
$
since $n prop a^(-3)$ and $lambda prop a$---during this we assume that photons are not created or destroyed, which they are, however $epsilon_"CMB"$ is so large that it dominates completely.

It turns out that the photon CMB and neutrino CMB together give
$
  Omega_(r,0) = Omega_("CMB",0)+Omega_(nu,0) = 9.00 times 10^(-5)
$
see pg. 72-73. This is very small in comparison to $Omega_(m,0) approx 0.31$ and $Omega_(Lambda,0) approx 0.69$---this is the Benchmark Model. In this model we have
$
  epsilon_(Lambda,0)/epsilon_(m,0) = Omega_(Lambda,0)/Omega_(m,0) approx 2.23
$
so today $Lambda$ dominates. In the past the ratio was
$
  (epsilon_Lambda (a))/(epsilon_m (a)) = epsilon_(Lambda, 0)/(epsilon_(m,0) a^(-3)) = Omega_(Lambda,0)/Omega_(m,0) a^3
$
since $epsilon_Lambda (a)$ is constant. The point when they were equal had
$
  a_(m Lambda) = (Omega_(m,0)/Omega_(Lambda,0))^(1\/3)
$
similarly
$
  (epsilon_m (a))/(epsilon_r (a)) = (epsilon_(m,0))/(epsilon_(r,0)) a => a_(r m) = (epsilon_(r,0))/(epsilon_(m,0))
$
the equation for $epsilon_i (a)$ tells us that in the limit $a -> 0$ components with big $w$ dominate, while for $a -> oo$ components with small $w$ dominate. This matches what observational evidence shows, namely that radiation ($w = 1\/3$) dominated, then matter $(w = 0)$, and then the cosmological constant $(w = -1)$. Given that $a$ is monotonically increasing with respect to $t$ then it is common to use it in place of time---similarly we can use the redshift $z$---especially since the conversion $a -> t$ is difficult.


== Empty universes
In the following sections we'll use the multi-component Friedmann equation
$
  dot(a)^2 = (8 pi G)/(3 c^2) sum_i epsilon_(i,0) a^(-1 - 3 w_i) - (kappa c^2)/R_0^2
$
where we've just plugged in the previous expression for $epsilon (a)$ and multipled by $a^2$.

We start with the simplest case, an empty universe. In this case the Friedmann equation becomes
$
  dot(a)^2 = - (kappa c^2)/R_0^2
$
one solution has $dot(a)=kappa=0$---so an empty, static, flat universe, i.e. the one described by the Minkowski metric. The other solution has $kappa = -1$ (Milne universe) giving
$
  dot(a) & = plus.minus c/R_0 \
         & => a(t) = t/t_0 "with" t_0 = R_0/c = H_0^(-1)
$
in this universe the light we observe at $t = t_0$ would be observed at some $t = t_e$ related by
$
  1 + z = 1/a(t_e) = t_0/t_e => t_e = (H_0^(-1))/(1 + z)
$
we also have
$
  d_p (t_0) = c integral_(t_e)^(t_0) dd(t)/a(t) = c t_0 integral_(t_e)^(t_0) dd(t)/t = c t_0 ln t_0/t_e = c/H_0 ln(1+z)
$
the distance at emission is smaller by a factor
$
  (a(t_e))/(a(t_0)) = 1/(1+z) => d_p (t_e) = c/H_0 ln(1+z)/(1+z)
$

== Single-component universes
Now we'll treat spacially flat universes with only one component having some $w$, in this case the Friedmann equation becomes
$
  dot(a)^2 = (8 pi G epsilon_0)/(3 c^2) a^(-(1+3 w))
$
to solve this we make the ansatz $a prop t^q$. In this case the LHS is $prop t^(2 q - 2)$ and the RHS is $prop t^(-(1+3 w)q)$ this is only true if
$
  q = 2/(3+ 3 w)
$
the scale factor is then
$
  a(t) = (t/t_0)^(2\/(3+3w)) "with" t_0 = 1/(1+w) (c^2/(6 pi G epsilon_0))^(1\/2)
$
and the Hubble constant would be
$
  H_0 equiv (dot(a)/a)_(t=t_0) = 2/(3(1+w)) t_0^(-1)
$
and the energy density would be
$
  epsilon(t) = epsilon_0 (t/t_0)^(-2)
$
letting $epsilon_0 = epsilon_(c,0)$ we can write
$
  epsilon(t) = 1/(6 pi (1+w)^2) c^2/G t^(-2)
$
Knowing $z$ we can find
$
  1 + z = (a(t_0))/(a(t_e)) = (t_0/t_e)^(2\/(3+3w)) => t_e = t_0/((1+z)^(3(1+w)\/2))
$
giving the proper distance
$
  d_p (t_0) = c/H_0 2/(1+3 w) [1 - (1+z)^(-(1+3w)\/2)]
$
in the limit $t_e = 0$ ($z = oo$) we get the horizon distance
$
  d_"hor" (t_0) = c integral_0^(t_0) dd(t)/a(t) = c/H_0 2/(1+3w)
$
the above is great for $w eq.not -1$, since then $q$ wouldn't be defined.

For a $Lambda$-dominated universe the Friedmann equation becomes
$
  dot(a)^2 = (8 pi G epsilon_Lambda)/(3 c^2) a^2 => dot(a) = H_0 a "with" H_0 = ((8 pi G epsilon_Lambda)/(3 c^2))^(1\/2)
$
clearly
$
  a(t) = e^(H_0 (t-t_0))
$
we can find
$
  d_p (t_0) = c/H_0 z => d_p (t_e) = c/H_0 z/(1+z)
$

== Multi-component universes
We can rewrite the Friedmann equation without explicit curvature as
$
  (H(t)^2)/H_0^2 = epsilon(t)/epsilon_(c,0) + (1-Omega_0)/a(t)^2
$
including matter ($w=0$), radiation ($w=1\/3$) and the cosmological constant ($w=-1$) we get
$
  H^2/H_0^2 = Omega_(r,0)/a^4 + Omega_(m,0)/a^3 + Omega_(Lambda,0) + (1-Omega_0)/a^2
$
in the Benchmark model $Omega_0 = Omega_(r,0) + Omega_(m,0) + Omega_(Lambda,0) = 1$---it is spatially flat.

Much of analyzing different component compositions consists of using seperation of variables to integrate the Friedmann equation giving $H_0 t = dots$, and then some expression for the scale factor pops out---it's also easy to find $a_"max"$. As an example take a universe with matter and a cosmological constant---then for $Omega_(Lambda, 0) < 0$ we get
$
  H^2/H_0^2 = Omega_(m,0)/a^3 + (1-Omega_(m,0)) => H_0 t = 2/(3 sqrt(Omega_(m,0)-1)) sin^(-1) [(a/a_"max")^(3\/2)]
$
where
$
  a_"max" = (Omega_(m,0)/(Omega_(m,0)-1))^(1\/3)
$
for $Omega_(Lambda, 0) > 0$ we get
$
  H_0 t = 2/(3 sqrt(1-Omega_(m,0))) ln[(a/a_(m Lambda))^(3\/2) + sqrt(1+(a/a_(m Lambda))^3)]
$
where
$
  a_(m Lambda) = (Omega_(m,0)/(1-Omega_(m,0)))^(1\/3)
$
this let's us calculate the approximate age of the universe with $Omega_(m,0) = 0.31$ and $Omega_(Lambda,0) = 0.69$ and
$
  t_0 = (2 H_0^(-1))/(3 sqrt(1-Omega_(m,0))) ln[(sqrt(1-Omega_(m,0))+1)/sqrt(Omega_(m,0))] = 13.74 plus.minus 0.40 "Gyr"
$
and $t_(m Lambda)$ or when $a=a_(m Lambda)$ is
$
  t_(m Lambda) = (2 H_0^(-1))/(3 sqrt(1-Omega_(m,0))) ln(1 + sqrt(2)) = 10.17 plus.minus 0.30 "Gyr"
$
which tells us that the cosmological constant has dominated for $approx 3.5 "Gyr"$. (with $Omega_(r,0) = 0$ and $kappa = 0$)

For a universe with just matter and radiation then
$
  t_(r m) = 4/3 (1 - 1/sqrt(2)) a^2_(r m)/sqrt(Omega_(r,0)) H_0^(-1) = 50000 "yr"
$
for the Benchmark model---thus radiation only dominated for a very short time and we can essentially ignore it.

== Benchmark model
\* skipped for now just a block of text.

Most of the important times ($t_(r m), t_(m Lambda)$ and $t_0$) have been listed, as well as the composition.

#pagebreak()
= Measuring Parameters
If we know the scale factor $a(t)$ we essentially know everything---so this is the quantity we want to measure. This is hard, so what we do is use a Taylor expansion
$
  a(t) = a(t_0) + dv(a, t)_(t=t_0) (t-t_0) + 1/2 dv(a, t, 2)_(t_t_0) (t-t_0)^2 + dots
$
assuming that $a$ isn't weird and fluctuates in a scuffed way we only really need the first couple of terms. We can write
$
  a(t) approx 1 + H_0 (t-t_0) - 1/2 q_0 H_0^2 (t-t_0)^2
$
with $H_0$ being the usual Hubble constant and $q_0$ being the deceleration parameter
$
  q_0 equiv - ((dot.double(a) a)/dot(a)^2)_(t=t_0) = - (dot.double(a)/(a H^2))_(t=t_0)
$
so we want to find $H_0$ and $q_0$---then we can approximate $a(t)$ near $t_0$.

Using the acceleration equation we can write
$
  -dot.double(a)/(a H^2) &= 1/2 [(8 pi G)/(3 c^2 H^2)] sum_(i=1)^N epsilon_i (1 + 3 w_i) \
  &= 1/2 sum_(i=1)^N Omega_i (1 + 3 w_i)
$
at the present moment then
$
  q_0 & = 1/2 sum_(i=1)^N Omega_(i,0) (1 + 3 w_i) \
      & = Omega_(r,0) + 1/2 Omega_(m,0) - Omega_(Lambda,0) approx^"b.mark" -0.53
$
In principle $H_0$ can be found using $c z = H_0 d$, but the distance is tough to work with. We can write
$
  1/a(t) approx 1 - H_0 (t-t_0) + (1 + q_0/2) H_0^2 (t-t_0)^2
$
then
$
  d_p (t_0) approx c (t_0-t_e) + (c H_0)/2 (t_0 - t_e)^2
$
here the first term is what the proper time would be in a static universe, and the second term is a correction since the universe expands. We don't know $t_0 - t_e$ but
$
  z = 1/a(t_e) -1 => z approx H_0 (t_0 - t_e) + (1+q_0/2) H_0^2 (t_0 - t_e)^2
$
so
$
  t_0-t_e approx H_0^(-1) [z-(1+q_0/2) z^2] => d_p (t_0) approx (c z)/H_0 [1-(1+ q_0)/2 z]
$
but $d_p (t_0)$ is not measurable.

So we need some way to measure distance from stuff we can measure---one way is using the luminosity distance
$
  d_L equiv (L/(4 pi f))^(1\/2)
$
which can give us the distance to some object (if the universe were static and flat) with a known luminosity $L$ using a measured flux $f$. In a universe described by the RW-metric the surface area of a sphere is given by
$
  A_p (t_0) = 4 pi S_kappa (r)^2
$
and the observed flux will be decreased by a factor $(1+z)^(-2)$ since the energy of a photon will be decreased by $(1+z)^(-1)$ (due to wavelength stretching) and the time between photons will be lengthened by $(1+z)$. The result is
$
  f = L/(4 pi S_kappa (r)^2 (1+z)^2) => d_L = S_kappa (r) (1+z)
$
given space seems flat $d_L = r(1+z) = d_p (t_0) (1+z)$, this gives
$
  d_L approx c/H_0 z (1 + (1-q_0)/2 z)
$
we can measure $d_L$, $q_0$ and $z$ so this gives us $H_0$ in terms of measurable quantities.

Another distance we need later when talking CMB is the angular-diameter distance defined by
$
  d_A equiv l/dd(theta, d: delta)
$
the distance between two ends of a body on the sky can be written as
$
  dd(s) = a(t_e) S_kappa (r) dd(theta, d: delta) = l => d_A = (S_kappa (r))/(1+z) = d_L/(1+z)^2 =^"kappa = 0" d_p (t_e)
$
so for $kappa = 0$ we have
$
  d_A (1+z) = d_p (t_0) = d_L/(1+z)
$
importantly for $z-> oo$ and $d_p (t_0) -> d_"hor" (t_0)$ we find
$
  d_L tilde z d_"hor" (t_0)\
  d_A tilde (d_"hor" (t_0))/z
$

#pagebreak()
= The Cosmic Microwave Background
The sky is uniformly bright at a $T_0 = 2.7255 "K"$ due to photons from the Big Bang (note that $T prop a^(-1)$). Currently
$
  epsilon_(gamma,0) = alpha T_0^4 = 0.2606 " MeVm"^(-3)
$
which is small, but since $h f_"mean"$ is small
$
  n_(gamma,0) = beta T^3 = 4.107 times 10^8 "m"^(-3)
$
which is large. Whereas for baryons
$
  epsilon_("bary",0) tilde 234 "MeVm"^(-3) -> n_("bary",0) tilde 0.25 "m"^(-3)
$
giving
$
  eta = n_("bary",0)/n_(gamma,0) tilde 6.1 times 10^(-10)
$
== Observing the CMB
The CMB was first observed as a constant isotropic noise---it is now known to be a blackbody, so the entire universe is itself a blackbody, around the microwave range.

We have three big observational results due to e.g. WMAP and COBE; at any angular position $(theta, phi.alt)$ on the sky, the spectrum of the CMB is a blackbody to within $10^(-4)$. The CMB also has a dipole distortion caused by our relative motion to it, blueshifting half the sky. When the dipole distortion is accounted for the remaining temperature fluctuation are minute in amplitude---at any point on the sky let the temperature be $T(theta,phi.alt)$ then
$
  expval(T) = 1/(4 pi) integral T(theta,phi.alt) sin theta dd(theta) dd(phi.alt) = 2.7255 "K"
$
the fluctuations are characterized by
$
  dd(T, d: delta)/T (theta,phi.alt) equiv (T(theta,phi.alt) - expval(T))/expval(T)
$
the deviations are of order
$
  expval((dd(T, d: delta)/T)^2)^(1\/2)
$
this is what we'd expect in a hot Big Bang model.

== Recombination and Decoupling
We want to know how $"ionized plasma" -> "gas of neutral atoms"$ in the early universe---and the related process of $"opaque" -> "transparent"$. We distinguish between three epochs; recombination, photon decoupling and last scattering.

We assume the only baryonic component is hydrogen, either in the form of neutral hydrogen or protons $p$. For charge neutrality in this universe we require $n_p = n_e$. The ionization can be characterized by
$
  X equiv n_p/(n_p+n_H) = n_p/n_"bary" = n_e/n_"bary"
$
in this universe the relevant energy scale is $Q = 13.6 "eV"$ a photon with energy $h f > Q$ can photoionize a hydrogen atom---$p$ and $e^-$ can also recombine
$
  H + gamma ->^"ionization" p + e^- ->^"recombination" H + gamma
$
the balance between these processes determine $X$. At a time when $a tilde 10^(-5)$ we had $T tilde 3 times 10^5 "K" -> h f_"mean" tilde 2.7 k_B T tilde 60 "eV"$. With $1.6$ billion photons per baryon all hydrogen atoms that would form by recombination are almost instantly ionized$-> X = 1$. At this time photons interacted with $e^-$ through Thomson scattering
$
  gamma + e^- -> gamma + e^-
$
this has $sigma_e = 6.65 times 10^(-29) "m"^2$ giving a mean distance before scattering
$
  lambda = 1/(n_e sigma_e) =>^"rate" Gamma = c/lambda = n_e sigma_e c
$
when the universe is fully ionized $n_e = n_p = n_"bary"$ and $n_"bary" = n_("bary",0) a^(-3)$ so the scattering rate for photons is
$
  Gamma = (n_("bary",0) sigma_e c)/a^3 = (5 times 10^(-21) "s"^(-1))/a^3
$
so for $a tilde 10^(-5) -> Gamma = 5 times 10^(-6) "s"^(-1)$ this is not very high---but photons are coupled to electrons as long as $Gamma > H$ or equivalently $lambda < c\/H$ so the mean free path is shorter than the Hubble distance. When the inequality flips the $e^-$ dilute faster than photons can interact with them---at this point the universe becomes transparent, and the baryonic matter is no longer compelled to have the same temperature as the CMB. For $a < a_(r m) tilde 2.9 times 10^(-4)$ we have
$
  H^2/H_0^2 = Omega_(r,0)/a^4 -> H = (H_0 Omega_(r,0)^(1\/2))/a^2
$
giving for $a tilde 10^(-5) -> H = 2.1 times 10^(-10) "s"^(-1)$ this is way smaller than the corresponding $Gamma$ so photons we strongly coupled to $e^-$. If hydrogen had remained ionized then photons would have remained coupled to protons and electrons until recently. One can find that if this were the case then decoupling would happen at $a tilde 0.025$ with $T tilde 110 "K"$, but at this temperature the CMB photons don't have enough energy to keet hydrogen ionized. For this reason decoupling is not gradually caused just by dilution of $e^-$---instead it is a sudden process caused by recombination.

== Recombination physics
Naively one can say that when the mean photon energy falls below $Q$ recombination happens---this gives
$
  T_"rec" tilde Q/(2.7 k_B) tilde 60000"K"
$
this is very crude, since the CMB is a blackbody with an exponential tail. We expect the true $T_"rec"$ be dependent on $eta$ and $Q$.

Consider the reaction
$
  H + gamma harpoons.rtlb p + e^-
$
this determines $X$ in the early universe. We assume thermal equilibrium and kinetic equilibrium---so the temperature is the same for everything and the particles adhere to either the Fermi-Dirac or Bose-Einstein distribution. Then the number density of a particle type $x$ in the momentum range $p -> p+dd(p)$ is given by (abuse of notation)
$
  n_x (p) dd(p) = g_x (4 pi)/h^3 (p^2 dd(p))/(exp([E-mu_x]/(k_B T))plus.minus 1)
$
with $+$ for fermions, and $-$ for bosons---$g_x$ is the statistical weight $gamma, e^- "and" p$ have $g_x = 2$. For photons with $E = p c = h f$ and $mu_gamma = 0$ we find
$
  n_gamma (f) dd(f) = (8pi)/c^3 (f^2 dd(f))/(exp[(h f)/(k_B T)]-1)
$
this is just the blackbody spectrum, integrating gives
$
  n_gamma = (2.4041)/pi^2 ((k_B T)/(hbar c))^3
$
note from lecture; assuming $E = p c$ and $mu_x << E(p)$ (radiation) we can find
$
  epsilon_x &= integral_0^oo dd(p) n_x (p) E \
  &= g_x c/(hbar^3 2 pi^2) integral_0^oo (dd(p) p^3)/(exp[p c\/k_B T] plus.minus 1) \
  &= g_x (k_B^4 T^4)/(hbar^3 c^3 2 pi^2) integral_0^oo dd(y) y^3/(exp(y) plus.minus 1) \
  &= g_x pi^2/30 (k_B^4 T^4)/(hbar^3 c^3) cases(1 ": boson", 7/8 ": fermion")
$
which is the familiar Stefan-Boltzmann law: $epsilon_gamma prop alpha T^4$, since $epsilon_r prop a^(-4)$ this implies $T prop a^(-1)$ as mentioned. We define (since over time different species stop being relativistic)
$
  epsilon_r (t) = sum_("ultra rel." x) epsilon_x (t) equiv g_* (t) pi^2/30 (k_B^4 T^4)/(hbar^3 c^3)
$
with $g_*$ being the number of (effective) bosonic relativistic degrees of freedoms:
$
  g_* (t) = sum_("rel. bosons" x) g_x + 7/8 sum_("rel. fermions" y) g_y
$

When recombination happened $e^-, p "and" H$ atoms all had non-relativistic speeds (no longer radiation) $m c^2 >> k_B T$ so $p tilde m_x v$ and
$
  E tilde m_x c^2 + p^2/(2 m_x)
$
so for $m_x c^2 - mu_x >> k_B T$ we find
$
  n_x (p) dd(p) = g_x (4 pi)/h^3 exp[(-m_x c^2 + mu_x)/(k_B T)] exp[- p^2/(2 m_x k_B T)] p^2 dd(p)
$
this is just the Maxwell-Boltzmann distribution, integrating gives
$
  n_x = g_x ((m_x k_B T)/(2 pi hbar^2))^(3\/2) exp[(-m_x c^2 + mu_x)/(k_B T)]
$
this gives
$
  epsilon_x tilde.eq m_x c^2 n_x tilde exp(-(m_x c^2 - mu_x)/(k_B T)) << 1 => epsilon_x << epsilon_r
$
this is consistent with the universe being radiation dominated during this time---until $mu_x tilde m_x c^2$. Note from lecture; we observe $n_overline(x) << n_x$ since we have matter, this implies that $mu_x eq.not mu_overline(x)$ since otherwise both $n_i$ would be exponentially supressed---see baryon-antibaryon assymmetry.

Now we assume that there was chemical equilibrium at the time of recombination. Given $mu_gamma = 0$ this means that
$
  mu_H = mu_p + mu_e
$
this gives
$
  n_H/(n_p n_e) = g_H/(g_p g_e) (m_H/(m_p m_e))^(3\/2) ((k_B T)/(2 pi hbar^2))^(-3\/2) exp[((m_p + m_e - m_H)c^2)/(k_B T)]
$
we now let $m_H\/m_p = 1$ and $Q = (m_p+m_e-m_H)c^2$ by definition and with $g_p=g_e=2$ and $g_H = 4$ we find
$
  n_H/(n_p n_e) = ((m_e k_B T)/(2 pi hbar^2))^(-3\/2) exp[Q/(k_B T)]
$
which is the Saha equation. We make the substitution
$
  n_H = (1-X)/X n_p
$
and with $n_e = n_p$ and $eta equiv n_"bary"\/n_gamma = n_p \/X n_gamma$ we find
$
  n_p = 0.2436 X eta ((k_B T)/(hbar c))^3
$
plugging everything into the Saha equation
$
  (1-X)/X^2 = 3.84 eta ((k_B T)/(m_e c^2))^(3\/2) exp[Q/(k_B T)]
$
this has a solution
$
  X = (-1 + sqrt(1+4S))/(2 S)",  " S(T,eta) = 3.84 eta ((k_B T)/(m_e c^2))^(3\/2) exp[Q/(k_B T)]
$
we define recombination by $X = 1\/2$ giving
$
  k_B T_"rec" = 0.324 "eV" = Q/42 -> T_"rec" = 3760 "K" tilde 250000"yr"
$
this was obviously not instantaneous but $X = 1\/2$ is reasonable. We can write
$
  Gamma (z) = n_e (z) sigma_e c = X(z) (1+z)^3 n_("bary",0) sigma_e c = 5 times 10^(-21) "s"^(-1) X(z) (1+z)^3
$
at this time the universe is matter-dominated so
$
  H^2/H_0^2 = Omega_(m,0)/a^3 = Omega_(m,0) (1+z)^3 => H(z) = 1.23 times 10^(-18) "s"^(-1) (1+z)^(3\/2)
$
so we find at $Gamma = H$
$
  1+z_"dec" = 39.3/(X(z_"dec")^(2\/3)) => z_"dec" = 1120
$
in reality the value is smaller---since the reaction is not in equilibrium and the real value is $z_"dec" tilde 1090 <-> T_"dec" = 2970 "K" <-> t_"dec" = 371000 "yr"$.

We define the optical depth by
$
  tau(t) = integral_t^(t_0) Gamma(t) dd(t)
$
this quantity is the expected number of scatterings a detected CMB photon at $t_0$ has undergone since $t$. The time of last scattering was when $tau = 1$. We can obtain
$
  tau(a) &= integral_a^1 Gamma (a) dd(a)/dot(a) = integral_a^1 (Gamma(a))/H(a) dd(a)/a \
  tau(z) &= integral_0^z (Gamma(z))/H(z) dd(z)/(1+z) = 0.0041 integral_0^z X(z) (1+z)^(1\/2) dd(z)
$
for our purposes $z_"ls" tilde z_"dec" tilde 1090$.

== Temperature fluctuations
The angular size $dd(theta, d: delta)$ of a temperature fluctuation in the CMB is related to a physical size $l$ of the last scattering surface by
$
  d_A = l/dd(theta, d: delta)
$
with $d_A$ being the angular-diameter distance to the last scattering surface. We have $z_"ls" >> 1$ so
$
  d_A tilde (d_"hor" (t_0))/z_"ls"
$
this gives
$
  d_A tilde 12.8 "Mpc" => l = 3.7 "kpc" (dd(theta, d: delta)/(1 "arcmin"))
$
the smallest WMAP has resolved are of size $dd(theta, d: delta) tilde 5 "arcmin"$ corresponding to $l tilde 18 "kpc"$ at $t_"ls"$ or $l(1+z_"ls") tilde 20 "Mpc"$ today.

To analyze fluctuations it's typically smart to do a Fourier transform---on a sphere this corresponds to using spherical harmonics
$
  dd(T, d: delta)/T (theta,phi.alt) = sum_(l=0)^oo sum_(m=-l)^l a_(l m) Y_(l m) (theta,phi.alt)
$
we are interested in the statistical properties of hot and cold spots---we use the correlation function $C(theta)$. We consider two points on the sky, relative to some observer they have directions $hat(n)$ and $hat(n)'$, these are seperated by $cos theta = hat(n) dot hat(n)'$. We define
$
  C(theta) equiv expval(dd(T, d: delta)/T (hat(n)) dd(T, d: delta)/T (hat(n)'))_(hat(n)dot hat(n)' = cos theta)
$
so we average over all points seperated by $theta$---so it gives us a way to characterize fluctuations in terms of scales. We can find
$
  C(theta) = 1/(4 pi) sum_(l=0)^oo (2 l + 1) C_l P_l (cos theta)
$
the $P_l$ are Legendre polynomials and the $C_l$ are the multiple moments of $C(theta)$. The $l= 0$ (monopole) term should vanish if the mean temperature is defined properly. The $l=1$ (dipole) term corresponds to the dipole distortion. For bigger $l$ the $C_l$ are a measure of temperature fluctuations on an angular scale $theta tilde 180 degree\/l$---so they are interchangeable---with bigger $l$ corresponding to fluctuations on smaller and smaller scales. It is typical to plot the power-spectrum given by
$
  Delta_T equiv ((l(l+1))/(2 pi) C_l)^(1\/2) expval(T)
$
this tells us the contribution per logarithmic interval in $l$ to the total temperature fluctuation $dd(T, d: delta)$.

== Cause of Fluctuations
We are interested in the horizon distance
$
  d_"hor" (t_"ls") = a(t_"ls") c integral_0^(t_"ls") dd(t)/(a(t))
$
we ignore dark energy since we are so early, so we take the universe to just contain matter and radiation. Then
$
  d_"hor" (t_"ls") = 2.24 c t_"ls" = 0.251 "Mpc"
$
this correspond to
$
  theta_"hor" = (d_"hor" (t_"ls") )/d_A tilde 1.1 degree -> l_"hor" tilde 160
$
this angle _splits_ the power-spectrum in two parts---on larger angular scales $theta > 4theta_"hor"$ it levels off at a nearly constant $Delta_T tilde 30 mu"K"$ on smaller angular scales $theta < theta_"hor"$ we see multiple peaks at $l tilde 220, 520, 800$.

For $theta > theta_"hor"$ the fluctuations are caused by gravitational effects due to density fluctuations in the distribution of nonbaryonic dark matter---since at this point $epsilon_"dm" > epsilon_gamma > epsilon_"bary"$, meaning the distribtuion of dark matter dominated the gravitational potential at the time of last scattering. Assuming the distribution was not perfectly homogeneous we can write
$
  epsilon(arrow(r)) = overline(epsilon) + dd(epsilon(arrow(r)), d: delta)
$
according to Newtonian mechanics the spatially varying deviation $dd(epsilon(arrow(r)), d: delta)$ leads to a varying potential $dd(Phi, d: delta)$ related by
$
  nabla^2 dd(Phi, d: delta) = (4 pi G)/c^2 dd(epsilon, d: delta)
$
so at the time of last scattering there would be gravitational fluctuations unless the distribution was perfectly homogeneous. CMB photons could be at either potential wells or hills (or somewhere in between) at the time of last scattering---this necessarily changes their energy since they either have to climb out of a potential well losing energy (redshift) or fall down a potnetial hill gaining energy (blueshift). So looking at the CMB hot spots correspond to maxima in $dd(Phi, d: delta)$ while cool spots correspond to minima---precisely
$
  dd(T, d: delta)/T = 1/3 dd(Phi, d: delta)/c^2
$
this is called the Sachs-Wolfe effect after the two guys who derived it. Given that $Delta_T$ is nearly constant for $l tilde 2 "to" l tilde 40$ tells us that $dd(Phi, d: delta)$ was constant across a range of scales.

For $theta < theta_"hor"$ the fluctuations are more complicated given that baryons and photons exist. Prior to decoupling these lived together in a photon-baryon fluid with energy density being $tilde 40%$ of $epsilon_"dm"$---so this fluid moved primarily due to gravity by dark matter. If this fluid happens to be in a well then it will start to collapse due to gravity, but as this happens the pressure in the fluid will increase due to photon-baryon interactions---this effect eventually wins and the fluid will begin to expand, until the pressure decreases sufficiently and it will collapse again. These cycles are called acoustic oscillations. If the fluid in a well is at maximum compression at the time of photon decoupling then its density will on average be higher, and since $T prop epsilon_gamma^(1\/4)$ the photons will be hotter than average. Similarly if the fluid is a maximum expansion then the photons will be cooler than average. If the fluid is in the process of expanding or collapsing then the photons will be cooler or hotter due to Doppler shifting. The first peak of $Delta_T$ at $l tilde 220 "or" theta tilde 0.8 degree$ corresponds to potential wells where the fluid had just reached maximum compression at the time of last scattering---these have size comparable to the sound horizon distance for the fluid at the time of last scattering. $d_s$ is the maximum proper distance that a sound wave in the fluid could have traveled since the Big Bang
$
  d_s (t_"ls") = a(t_"ls") integral_0^(t_"ls") (c_s (t) dd(t))/a(t)
$
we can approximate $c_s tilde c\/sqrt(3)$---the speed of sound in a pure photon gas. This gives
$
  theta_s tilde (d_s (t_"ls"))/d_A tilde 0.7 degree
$
this assumes the universe is flat, since the peak would shift given curvature. The observed first peak is consistent with $kappa = 0$. The amplitude of the peak is dependent on $c_s$ with lower speeds giving higher amplitudes---by definition $c_s = sqrt(w_"pb") c$ with $w_"pb"$ dependent on $eta$ the photon-baryon ratio. So the first peak gives us information about both curvature and composition of the universe.

#pagebreak()
= Particle Physics
== The observation
An easy observation to make is that our universe contains different stuff. For us the most significant difference between stuff is what elementary particles make them up. Baryonic matter corresponds to protons, neutrons, and electrons---since the electrons weigh so little---most of this is found as hydrogen and helium (we also have dark matter, which we won't discuss). We also have three types of neutrinos $nu$ and three mass states---these have very little mass and are not very reactive.

The important mass-less particle for our purposes is the photon, which unlike the neutrino is very reactive---for blackbody radiation we know
$
  epsilon(f) dd(f) = (8 pi h)/c^3 (f^3 dd(f))/(exp[h f\/k_B T]-1)
$
with a peak at $h f_"peak" approx 2.82 k T$, and
$
  epsilon_gamma = alpha T^4
$
the number density is just
$
  n(f) dd(f) = (epsilon (f) dd(f))/(h f) = (8 pi)/c^2 (f^2 dd(f))/(exp[h f\/k_B T]-1)
$
so $n_gamma = beta T^3 => E_"mean" approx 2.70 k_B T$---the blackbody is important since our universe is essentially just one big blackbody (see CMB).

== Nucleosynthesis and the Early Universe
Before the time of the last scattering surface $t_(1s) approx 0.37 "Myr"$ the universe was opaque, meaning we can't see what the universe was like.

In the very early universe radiation dominated at times $t << t_(r m) approx 50000 "yr"$ here $a(t) prop t^(1\/2)$ and the temperature of blackbody photons in the universe with $T prop a^(-1)$ is given by
$
  T(t) prop t^(-1\/2) -> E_"mean" (t) = 2.7 k_B T(t)
$
using actual numbers we find that the LHC can achieve energies corresponding to a time of $t tilde 10^(-13) "s"$.

As the universe expanded and cooled the energy dropped from
$
  E_"mean" (t_P) tilde E_P tilde 10^(28)"eV" -> E_"mean" (t_0) tilde 10^(-3)"eV"
$
this is a very wide range---and obviously some scales are of more interest according to which process is being studied. Recombination and photoionization occur at energies $tilde 10"eV"$, while fission and fusion occur at much higher energies.

When talking atomic nuclei we use the mass number $A = Z + N$, where $Z$ is the amount of protons, and $N$ is the amount of neutrons---two of the most importan are $p$ (proton, hydrogen) and $D$ (deuterium). The binding energy $B$ of a nucleus is the energy required to pull it apart into its component protons and neutrons---for example $B_D = 2.22 "MeV"$
$
  p + n harpoons.rtlb "D" + 2.22 "MeV"
$
the relative binding energy can be written as $B\/A$, the usual energies here are $tilde 8 "MeV"$. This shows us that before recombination at $tilde 10 "eV"$ there was a period of nucleosynthesis (BBN) when protons and neutrons began to fuse and form deuterium, which could then fuse and form heavier nuclei.

The building blocks for nucleosynthesis are $n$ and $p$---note
$
  Q_n = (m_n - m_p )c^2 = 1.29 "MeV"
$
a free neutron is unstable and decays
$
  n -> p + e^- + overline(nu)_e
$
with decay time $tau_n = 880 "s"$---given that $Q_n > m_e c^2 = 0.51 "MeV"$ some energy is carried away by the kinetic energy of the electron and the antineutrino---this decay time is short on cosmological scales so free neutrons quickly went extinct. When $E_"mean" tilde 10 "MeV"$ at $t = 0.1"s"$ we had positrons and electrons since this energy is much larger then their rest energies, these were created by pair production
$
  gamma + gamma harpoons.rtlb e^- + e^+
$
also at this time neutrinos were still coupled to neutrons and protons so
$
  n + nu_e & harpoons.rtlb p + e^- \
   n + e^+ & harpoons.rtlb p + overline(nu)_e
$
in the universe everything was in kinetic equilibrium with $k_B T << m_p c^2$---using this eventually leads to
$
  n_n/n_p tilde.eq exp[- Q_n/(k_B T)]
$
but at some point the neutrons and protons won't be in equilibrium---the interaction between them requires a neutrino or antineutrino. These interact through the weak force which has a relatively small cross-section and $sigma_w prop t^(-1)$ combined with $n_nu prop t^(-3\/2)$ the interaction rate is
$
  Gamma prop n_nu sigma_w = t^(-5\/2)
$
when $Gamma tilde H$ the neutrino decouples from the protons and neutrons and hence there is no longer a conversion---the ratio is frozen---it can be shown that
$
  n_n/n_p = exp[- Q_n/(k_B T_"freeze")] approx 0.2
$
and $t_"freeze" tilde 1"s"$ with $k_B T tilde.eq 0.8 "MeV"$---the ratio would be valid for $t_"freeze" < t << tau_n$. This is one reason why BNN is quite inefficient since $p + n -> "D"$ is easy but $p + p -> "D"$ is a two-step process
$
  p + p harpoons.rtlb isotope("He", a: 2) -> "D" + e^+ + nu_e
$
note that $isotope("He", a: 2)$ is really, really unstable and decays into two $p$ with $tau_"split" tilde 10^(-23) "s"$---so fusing two $p$ is really rare. For this reason we can say that BNN stops when every neutron is bound---with the remaining protons just being alone.

After freezeout at $t tilde 2"s"$ the first step is as mentioned production of deuterium
$
  p + n harpoons.rtlb "D" + gamma
$
with $B_D$ being carried away by a gamma photon---to analyze this we can use the nucleosynthetic Saha equation
$
  n_D/(n_p n_n) = 6 ((m_n k_B T)/(pi hbar^2))^(-3\/2) exp[B_D/(k_B T)]
$
the time of nucleosynthesis defined by $n_D\/n_n = 1$ can be found to be $t_"nuc" tilde 200"s"$.

After deuterium is formed many other reactions are possible, most leading to the formation of helium. We have
$
  D + p & harpoons.rtlb isotope("He", a: 3)+gamma \
    D+n & harpoons.rtlb isotope("H", a: 3)+gamma
$
or
$
  D + D harpoons.rtlb isotope("He", a: 4) + gamma
$
more likely
$
  D + D & harpoons.rtlb isotope("H", a: 3) + p \
  D + D & harpoons.rtlb isotope("He", a: 3) + n
$
which then quickly fuse
$
   isotope("H", a: 3) + p & harpoons.rtlb isotope("He", a: 4) + gamma \
  isotope("He", a: 3) + n & harpoons.rtlb isotope("He", a: 4) + gamma \
   isotope("H", a: 3) + D & harpoons.rtlb isotope("He", a: 4) + n \
  isotope("He", a: 3) + D & harpoons.rtlb isotope("He", a: 4) + p \
$
when helium is reached however nucleosynthesis has trouble since it is very tightly bounded relative to its mass number. Further, there exist no stable nuclei with $A=5$ meaning it can't fuse with $p$ or $n$ and hope it achieves something. It can form small amount of lithium through
$
  isotope("He", a: 4) + D &harpoons.rtlb isotope("Li", a: 6) + gamma \
  isotope("He", a: 4) + isotope("H", a: 3) &harpoons.rtlb isotope("Li", a: 7) + gamma
$
and beryllium
$
  isotope("He", a: 4) + isotope("He", a: 3) &harpoons.rtlb isotope("Be", a: 7) + gamma
$
again there are no stable nuclei for $A=8$ so we reach another roadblock---for obvious reasons this is way more complicated, e.g. temperature dependent cross-sections.

One question still being researched today can be expressed as
$
  n_"antibary" << n_"bary" << n_gamma
$
why is there a baryon-antibaryon asymmetry? In the early universe quarks and antiquarks were constantly formed through pair-production and quickly annihilated. However, suppose there was a tiny asymmetry $delta_q << 1$ leading to
$
  n_q/n_gamma tilde delta_q
$
If we had $1003$ quarks per $1000$ antiquarks then for every three surviving quarks there'd be $2000$ photons---the three quarks would when the universe cooled down form a single baryon giving a very low $eta$. In reality we just need $800000003$ quarks per $800000000$ antiquarks to get the $eta$ we observe.

