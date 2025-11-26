//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/mannot:0.3.0": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Friedmann-Robertson-Walker metric
In the preceding part of the course we developed all the essential tools of relativity and applied them to the canonical example---the black hole. The second part of this course is focused on the application of relativity to cosmology.

=== Basic assumptions
At _large enough_ scales clusters of galaxies etc. can be treated like pointlike particles in a continuous cosmic fluid. This combined with the assumption that there is no preferred point in the Universe leads to the cosmological principle:

_The universe is homogeneous and isotropic._

We also assume there exists some cosmological time, implying the Universe has a rest frame---i.e. at some fixed time $t$ matter is at rest. The motion of matter (the cosmic gas) would then follow this cosmological time. This also implies we can define a set of co-moving coordinates.

These two assumptions: the cosmological principle and the existence of cosmological time, make up the cosmological standard model.

== The metric
Using the above assumptions we will now attempt to construct the co-moving metric.

By the second assumption we consider two spatial slices at different fixed cosmic times $t$, denoted by $S_0$ and $S^*$ respectively---any point $P_0 (x^1,x^2,x^3) in S_0$ evolves in time to become a point $P^* (x^1,x^2,x^3) in S^*$. In the rest-frame $dd(t^2) = dd(tau^2)$ implying $g_(0 0) = -1$. Further we should not experience any mixing of coordinates $g_(0 i) = 0$ since we should be able to move along $t$ without changing our spatial position. We therefore consider a metric of the form
$
  dd(tau^2) = dd(t^2) - g_(i j) dd(x^i, x^j) "with" x^i in S_0
$
before trying to determine $g_(i j)$ we check if this is consistent. Consider the geodesic equation
$
  dv(x^mu, tau, 2) + tensor(Gamma, mu, -lambda rho) dv(x^lambda, tau) dv(x^rho, tau) = 0
$
on $S_0$ this implies
$
  tensor(Gamma, i, -00) = 0 => 2 pdv(g_(0 i), t) = underbrace(pdv(g_00, x^i), 0)
$
so the $g_(0 i)$ will be zero for all time, and therefore the metric is consistent. Now consider the metric at some spatial slice
$
  dd(s^2) = g_(i j) dd(x^i, x^j)
$
since $g_(i j)$ can depend on $t$ then $dd(s^2)$ can change over time, but the coordinates $x^i$ are still fixed. This implies a time-dependent rescaling of $dd(s^2)$.

By the cosmological principle we write
$
  dd(tau^2) = dd(t^2) - f_1 (t,r) dd(r^2) - f_2 (t,r) r^2 dd(Omega^2)
$
we decompose $f_1$ and $f_2$ to write
$
  dd(tau^2) = dd(t^2) - f(t) L(r) dd(r^2) - g(t) H(r) r^2 dd(Omega^2)
$
redefining $r$ we can absorb $H$ to obtain
$
  dd(tau^2) = dd(t^2) - f(t) L(r) dd(r^2) - g(t) r^2 dd(Omega^2)
$
to determine $f, L$ and $g$ we use the EFE. Again by the cosmological principle it is reasonable to assume that the Universe is a perfect fluid meaning
$
  T_(mu nu) = p g_(mu nu) + (p + rho) U_mu U_nu
$
in the co-moving frame $U^i = 0$ and $U^t = 1$ so
$
  R_(t r) = 0 => dot(g)/g = dot(f)/f
$
implying $g = f$ up to a constant---so the metric becomes
$
  dd(tau^2) = dd(t^2) - f(t) [L(r) dd(r^2)+r^2 dd(Omega^2)]
$
To determine $L(r)$ notice that $tensor(R, r, -r)$ cannot depend on $r$ so the spatial part must be a constant giving
$
  - L'/(r L^2) = 2 k => L = 1/(tilde(k)-k r^2)
$
similarly for $tensor(R, theta, -theta)$ giving
$
  1/r^2 - 1/(r^2 L) + L'/(2 r L^2) = 2 k
$
since the curvature is spatially constant these are the same. Then we can solve to find
$
  (1 - tilde(k))/r^2 + 2k = 2 k => tilde(k) = 1
$
putting everything together and defining the scale factor $a equiv sqrt(f)$ we find the FRW metric
$
  dd(tau^2) = dd(t^2) - a^2 (t) [1/(1-k r^2) dd(r^2)+r^2 dd(Omega^2)]
$
to see the meaning of $k$ one can consider the spatial curvature by
$
  tensor(R, i, -i) = g^(i k) R_(k i) = - (6 k)/a^2 = - 6 K
$
where $K = k a^(-2)$ is the Gaussian curvature---importantly we see that $k$ determines the spatial curvature.

=== Geometric interpretation of $bold(k)$
We construct an embedding diagram by
$
  dd(s^2) = (dd(bold(x)))^2 + (dd(z))^2 => bold(x)^2 + z^2 = R^2
$
this describes a hypersphere. Fixing any angular coordinate would then reduce the dimensionality and give us a normal three-sphere. By the above:
$
  dd(z^2) = 2 z dd(z) = - dd(bold(x)^2)
$
meaning
$
  dd(s^2) & = (dd(bold(x)))^2 + (dd(bold(x)^2))^2/(4(R^2-bold(x)^2)) \
          & =^"polar" dd(r^2) + r^2 dd(Omega^2) + (r^2 dd(r^2))/(R^2-r^2) \
          & = 1/(1-r^2\/R^2) dd(r^2) + r^2 dd(Omega^2)
$
defining $rho equiv r\/R$ this becomes
$
  dd(s^2) = R^2 (dd(rho^2)/(1-rho^2) + rho^2 dd(Omega^2))
$
if we define $R equiv a(t)$ at some fixed $t$ then this is the FRW metric if $rho^2 = k r^2$ meaning $k > 0$. So for $k > 0$ the FRW describes a sphere! A similar calculation shows that if one uses
$
  - bold(x)^2 + z^2 = R^2
$
then $k <^! 0$, and obviously for $k = 0$ the FRW metric is spatially flat.

#pagebreak()
== Hubble's law
Above we found
$
  dd(tau^2) = dd(t^2) - a^2 (t) [dd(r^2)/(1-k r^2) + r^2 dd(Omega^2)]
$
Consider an observer standing far away from a galaxy. Light follows null geodesics $dd(tau)= 0$. We can choose coordinates such that $dd(Omega^2) = 0$ so
$
  dd(t^2) = a^2 (t) dd(r^2)/(1-k r^2) => integral dd(t)/a(t) = integral dd(r)/sqrt(1-k r^2)
$
We take the light to leave the galaxy at radial distance $r_1$ and time $t_1$ and arrive at $r = 0$ and time $t_0$ (so observer is at the origin). We can then integrate to find
$
  integral_(t_1)^(t_0) dd(t)/a(t) = integral_0^(r_1) dd(r)/sqrt(1-k r^2) = cases(sin^(-1) r_1 #h(2em) &"for" k = +1, r_1 &"for" k = 0, sinh^(-1) r_1 &"for" k=-1)
$
We now consider two crests of light emitted at $(t_1, r_1)$ and $(t_1 + dd(t_1, d: delta), r_1)$, the observer would observe these at $(t_0, r=0)$ and $(t_0+dd(t_0, d: delta), r=0)$. So we have
$
  integral_(t_1+dd(t_1, d: delta))^(t_0+dd(t_0, d: delta)) dd(t)/(a(t)) = integral_(t_1)^(t_0) dd(t)/a(t)
$
taking $dd(t_i, d: delta)$ to be small so we treat $a(t_i)$ as a constant (our main assumption!),
$
  (t_0+dd(t_0, d: delta))/a(t_0) - (t_1 + dd(t_1, d: delta))/a(t_1) = t_0/a(t_0) - t_1/a(t_1) => dd(t_0, d: delta)/a(t_0) = dd(t_1, d: delta)/a(t_1)
$
so we have a redshift! In terms of frequencies
$
  lambda_1/lambda_0 = nu_0/nu_1 = dd(t_1, d: delta)/dd(t_0, d: delta) = a(t_1)/a(t_0)
$
and if the Universe is expanding $a(t_1) > a(t_0)$ so $lambda_1 > lambda_0$. We define the redshift
$
  z & equiv (lambda_0 - lambda_1)/lambda_1 \
  & = a(t_0)/a(t_1) -1 \
  & tilde.eq^"small time" underbrace((dot(a) (t_0))/a(t_0), equiv H_0) (t_0-t_1)
$
with $H_0$ being the Hubble rate today. We can relate this to a distance
$
  r_1 = integral_(t_1)^(t_0) dd(t)/a(t) tilde.eq (t_0-t_1)/a(t_0)
$
The physical distance is defined as $L = a(t_0) r$, so we can write
$
  z = H_0 L
$
taking the derivative of $L$ then gives
$
  v_r = dot(a) (t_0) r_1 = (dot(a) (t_0))/a(t_0) underbrace(a(t_0) r_1, t_0-t_1) = z
$
this is Hubble's law
$
  z = v_r = H_0 L
$
the observational verification of this essentially proves the Universe is expanding.

#pagebreak()
== Friedmann equations
We want an equation for the scale factor $a(t)$. To do this we use the EFEs
$
  tensor(R, mu, - nu) = 8 pi G (tensor(T, mu, - nu) - 1/2 T tensor(delta, mu, -nu))
$
with
$
  T_(mu nu) = p g_(mu nu) + (p + rho) U_mu U_nu
$
and $U^t = 1$, $U^i = 0$. We can find
$
  tensor(T, t, -t) = p-(p+rho)";  " tensor(T, r, -r) = p
$
and
$
  tensor(R, t, -t) = -3 dot.double(a)/a";  " tensor(R, r, -r) = - dot.double(a)/a - 2 dot(a)/a - (2 k)/a^2
$
we find two equations
$
                     3 dot.double(a) & = - 4 pi G (rho+3 p) a \
  a dot.double(a) + 2 dot(a)^2 + 2 k & = 4 pi G ( rho -p) a^2
$
these can be rewritten slightly
$
     (dot(a)/a)^2 & = (8 pi G)/3 rho - k/a^2 \
  dot.double(a)/a & = - 4 pi G (p + 1/3 rho)
$
these are the Friedmann equations. We can also use $D_mu T^(mu nu) = 0$ (but does not give new information) to find
$
               a^3 dv(p, t) & = dv(, t) [a^3 (p+rho)] \
          dv(, a) (rho a^3) & = - 3 p a^2 \
  dot(a) pdv(, t) (rho a^3) & = - 3 p a^2
$
The Friedmann equations are used to find how a Universe with a given composition will evolve.

We still need an equation of state $p = p(rho)$:

1. $p << rho$ (pressure-less matter). By the above we immediately find
$
  rho prop 1/a^3
$
this is also called non-relativistic matter.

2. $p = 1/3 rho$ (radiation). By the above
$
  rho prop 1/a^4
$
this is also called relativistic matter.

3. $p = - rho$. By the above
$
  rho prop "const"
$
as we would expect for a cosmological constant $Lambda$.

All of the above have the form $p = w rho$ (ideal fluid) with $w = \{0, 1/3, -1}$.

=== Age of the Universe
Consider $p >= 0$ then from
$
  dot.double(a)/a = - 4 pi G (p + 1/3 rho)
$
we can see $dot.double(a)\/a < 0$ and $dot(a)\/a < 0$. Meaning at some point $a$ was zero. We define the time of the beginning of the Universe (the Big Bang) as $a(t=0)=0$.

We can then compute the age of the Universe by
$
  t_0 = integral dv(t, a) dd(a) = integral dd(a)/dot(a) < integral_0^(a(t_0)) dd(a)/(dot(a) (t_0))
$
since $dot(a) (t_0) < dot(a) (t)$. Then we find
$
  t_0 < 1/H_0
$
so the Hubble rate is the upper-bound for the age of the Universe (if we have non-relativistic and relativistic matter with $p >= 0$)---using modern values we find $t_0 lt.tilde 13.8 "Gyr"$.

Consider
$
  dv(, a) (rho a^3) = - 3 p a^2 " for " p >= 0
$
so $rho a^3$ decreases, meaning for $a -> oo$ we have
$
  rho a^2 -> 0
$
Then
$
  (dot(a)/a)^2 = (8 pi G)/3 rho - k/a^2 => dot(a)^2 = (8 pi G)/3 rho a^2 - k
$
so $dot(a)^2 -> - k$. Which is physical for $k = - 1$ (keeps expanding) or $k = 0$ (asymptotic slow down). For $k = +1$ we have $dot(a) -> 0$ for $a -> a_"max"$ with
$
  rho a_"max"^2 = 3/(8 pi G)
$
and then the Universe would collapse on itself.

=== Critical density
For $k=0$ we have
$
  H^2 =(dot(a)/a)^2 = (8 pi G)/3 rho => rho_c equiv (3 H^2)/(8 pi G)
$
so the critical density is the density of energy in the universe if it were flat. We then define
$
  Omega = rho/rho_c
$
for $k eq.not 0$ we can then write the curvature as
$
  k/a^2 = H^2 (Omega-1)
$
meaning for $k = 0$ we have $Omega = 1$.

We can write the Friedmann equation as (see cosmology notes)
$
  H^2/H_0^2 = Omega_(m,0)/a^3 + Omega_(r,0)/a^4 + Omega_(Lambda,0) + (1- Omega_0)/a^2
$

#pagebreak()
= Evolution of the Universe
Above we had
$
  D_mu T^(mu nu) = 0 => dv(, a) (rho a^3) = - 3 p a^2
$
and we consider three types of matter: non-relativistic, relativistic and the cosmological constant. For these we found
$
       rho_m & prop a^(-3) \
       rho_r & prop a^(-4) \
  rho_Lambda & prop "const"
$
then plotting $(ln a, ln rho)$ gives the usual diagram highlighting periods where one type of matter dominated the rest. At some time $rho_m = rho_r$ and at some later time $rho_m = rho_Lambda$ by interpolating the above---periods of equality. At late enough times $rho_Lambda$ dominates completely.

== Big Bang
Going backwards in time it is clear that the energy density of the Universe increases, so early in the Universe everything was very hot and dense. Then it is natural to assume that particles through interactions would be in thermal equilibrium at early times.

The number density of bosons (fermions) would then be given by
$
  dd(n) = (4 pi g)/(h^3 c^3) (E^2 dd(E))/(e^(E\/k_B T) plus.minus 1)
$
and this is naturally related to $rho$. We compute
$
  n & = (4 pi g)/(h^3 c^3) integral (E^2 dd(E))/(e^(E\/k_B T) plus.minus 1) \
  n & prop cases(3/2 T^3 zeta(3) "for fermions", 2 T^3 zeta(3) "for bosons")
$
or we can write
$
  n_b = 4/3 n_f = 2.404 g/(2 pi^2) ((k_B T)/(h c))^3
$
so in the non-relativistic case when $rho prop n$ we find $T prop a^(-1)$. We can also compute the energy density
$
  rho = integral E dd(n) prop T^4
$
so in the relativistic case we also have $T prop a^(-1)$. In this case one can show by solving the Friedmann equation that $a prop sqrt(t)$.

Combining these we find for a radiation-dominated universe
$
  k_B T tilde.eq 0.46 E_(p l) (t/t_(p l))^(-1\/2)
$
with $E_(p l)$ and $t_(p l)$ being the Planck energy and time respectively. At these scales relativity breaks down, our description only makes sense for time scales longer than $t_(p l)$. The point $t = t_(p l)$ is what we mean by Big Bang.

=== BBN
At temperatures $T > "MeV"$ we had equilibrium with
$
  p^+ + n <--> D^+ + dots
$
at lower temperatures $T < "MeV"$ the formed nuclei could no longer be ripped apart so we have reactions like
$
  p^+ + n --> D^+ --> "He" --> "Li"
$
etc. This is the formation of the light elements.

#pagebreak()
== The cosmic microwave background
At temperatures $T gt.tilde 13.6 "eV"$ we had
$
  e^- + p^+ <--> H + gamma
$
at lower temperatures $T lt.tilde "eV"$ the free electrons and protons combine to form neutral hydrogen
$
  e^- + p^+ --> H + gamma
$
this process is called recombination. (recombination $z tilde 1370$ $->$ photon decoupling $->$ last scattering $z tilde 1100$ or $t tilde 380 "kyr"$)

Before recombination photons scatter with free electrons and protons forming a coupled plasma. At this point the Universe is opaque to photons. After recombination photons can move more freely since they do not scatter with neutral hydrogen making the Universe transparent to photons. These photons are what we call the CMB. The temperature of the CMB photons today is $T = 2.73 "K"$, this is found by fitting to a blackbody, since the CMB is a blackbody. Since we know the energy at which recombination happens we can also predict the temperature, and this is in agreement with the measured value.

=== Anisotropies
When we observe the CMB we observe the last scattering surface, through some angle. One of the primary observations are temperature anisotropies on the scale of $delta T\/T tilde 10^(-5)$. This in part confirms that the Univserse is isotropic, but at small scales we observe anisotropies.

This is done by essentially doing a Fourier transform on a sphere, where the Fourier basis are the spherical harmonics $Y_l^m$. The multipole $l$ is what we really care about in this case. For the monopole $l = 1$ we essentially just have the mean $T = T_"mean"$. For the dipole $l = 2$ we observe $delta T\/T tilde 10^(-3)$ which is caused by our relative motion to the CMB. For higher multipoles we observe the very small fluctuations. The higher multipoles correspond to smaller angular separations on the last scattering surface $delta theta$.

Considering a length $cal(l)$ on the last scattering surface we can write
$
  sin(delta theta) = cal(l)/d_A => delta theta = cal(l)/d_A
$
By the FRW-metric we have
$
  dd(s) = a (t_e) r delta theta =>^"small angles" cal(l) tilde.eq a (t_e) r delta theta
$
we can rewrite this as
$
  cal(l) & = a(t_e)/a(t_0) a(t_0) r delta theta \
         & = (a(t_0) r delta theta)/(1+z)
$
We define the horizon $d_"hor"$ to be the distance a photon could have travelled since the Big Bang
$
  d_"hor" & = a(t_0) integral_0^(t_0) dd(t)/(a(t)) \
          & = a(t_0) r
$
using this we can rewrite $cal(l)$
$
  cal(l) = (d_"hor" delta theta)/(1+z) =>^(z >> 1) d_A = d_"hor"/z
$
we can also compute $d_"hor"$ giving
$
  d_"hor" tilde.eq 1.4 times 10^(4) "Mpc"
$
this is approximately what one would find ignoring the expansion of the Universe. Then we can find the physical distance at last scattering by
$
  d_A tilde.eq 13 "Mpc" => cal(l) = d_A delta theta tilde.eq 0.22 "Mpc" ((delta theta)/(1 degree))
$
Fluctuations with $delta theta > 7 degree$ have been observed. At last scattering these correspond to length scales $cal(l) > 1.6 "Mpc"$. Today these correspond to the length scale $cal(L) > 1700 "Mpc"$ which is comparable to the size of the Universe. Therefore the fluctuations appear to be an initial condition of the Universe.

As written before we can write the fluctuations as a Fourier transform
$
  (delta T)/T (theta,phi) = sum_(l=0)^oo sum_(m=-l)^(m=l) a_(l m) Y_(l m) (theta, phi)
$
we are interested in the correlation function
$
  C(theta) = evaluated(expval((delta T)/T (hat(n)_1) (delta T)/T (hat(n)_1)))_(hat(n)_1 dot hat(n)_2 = cos theta)
$
this is essentially the average of the correlation between all pairs of points seperated by an angle $theta$. Given the fluctuations are Gaussian one can show that
$
  expval(a_(l m))=0";  " expval(a_(l m)^* a_(l' m')) = overbrace(C_(l), "power") delta_(l l') delta_(m m')
$
then we can write
$
  C(theta) = 1/(4 pi) sum_(l=0)^oo (2 l + 1) C_l P_l (cos theta)
$
so we have removed the $m$-dependence. From this one can construct the power-spectrum by plotting $l(l+1) C_l$ against $l$. The angle is related to a length scale since for smaller angles (or higher $l$) we can measure smaller lengths. As an example for $l = 2$ we can only split the sky in two, and we only have one Universe. For this reason the uncertainty for low $l$ is massive since we just have less points (cosmic variance). For higher multipoles we can take many different measurements, allowing us to decrease or eliminate the variance.

\* $-dd(f) = underbrace(f(r), "fraction uncovered") dd(Omega)$
