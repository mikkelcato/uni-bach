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
In a homogeneous and isotropic universe we can use
$
  (dot(a)/a)^2 & = (8 pi G)/(3 c^2) epsilon - (kappa c^2)/(R_0^2 a^2) \
             0 & = dot(epsilon) + 3 dot(a)/a (epsilon + P) \
             P & = w epsilon
$
to relate $epsilon(t)$, $P(t)$ and $a(t)$.

=== Energy density
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


=== Empty universes
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

=== Single-component universes
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

=== Multi-component universes
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

=== Benchmark model
\* skipped for now just a block of text.

Most of the important times ($t_(r m), t_(m Lambda)$ and $t_0$) have been listed, as well as the composition.

#pagebreak()
== Measuring Parameters
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

#pagebreak()
== Particle Physics
=== The observation
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

=== Nucleosynthesis
