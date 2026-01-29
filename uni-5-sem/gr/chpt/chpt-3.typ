//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/mannot:0.3.0": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Friedmann-Robertson-Walker metric
== Basic assumptions
At _large enough_ scales clusters of galaxies etc. can be treated like pointlike particles in a continuous cosmic fluid. This combined with the assumption that there is no preferred point in the Universe leads to the cosmological principle:

_The universe is homogeneous and isotropic._

We also assume there exists some cosmological time implying the Universe has a rest frame. So at some fixed time $t$ matter is at rest. The motion of the cosmic fluid would then follow this cosmological time.

== The metric
We describe the above by considering Gaussian coordinates (the synchronous gauge). Consider two spatial slices at fixed cosmic times $t$ denoted by $S_0$ and $S_*$. We let every point $P_0 (x^i)$ evolve along a geodesic perpendicular to $S_0$ to become points $P_* (x^i)$. Along such a geodesic the $x^i$ are fixed (comoving) and we take $t = tau$, meaning
$
  dd(tau^2) = dd(t^2)
$
implying $g_(00) = -1$. We also require $g_(i 0) = g_(0 i) = 0$ if the geodesics are to remain perpendicular. Then
$
  dd(tau^2) = dd(t^2) - g_(i j) dd(x^i, x^j) "with" x^i in S_0
$
We check this is consistent. Consider
$
  dv(x^mu, tau, 2) + tensor(Gamma, mu, -lambda rho) dv(x^lambda, tau) dv(x^rho, tau) = 0
$
For a point on $S_0$
$
  tensor(Gamma, i, -00) = 0
$
implying
$ 2 pdv(g_(0 i), t) = underbracket(pdv(g_00, x^i), 0) $
so $g_(0 i)$ will be zero for all time. Consider a spatial slice ($dd(t)^2 = 0$)
$
  dd(s^2) = g_(i j) dd(x^i, x^j)
$
Since $g_(i j)$ can depend on $t$ then $dd(s^2)$ can change over time but the $x^i$ are still fixed. This implies a time-dependent rescaling of $dd(s^2)$.

Using isotropy we require
$
  dd(tau^2) = dd(t^2) - f_1 (t,r) dd(r^2) - f_2 (t,r) r^2 dd(Omega^2)
$
we can decompose $f_1$ and $f_2$
$
  dd(tau^2) = dd(t^2) - f(t) L(r) dd(r^2) - g(t) H(r) r^2 dd(Omega^2)
$

#proof[
  Consider two points at the same $t$ with $dd(theta) = dd(phi) = 0$. The line elements of these points are $f_1 (t,r_1) dd(r^2)$ and $f_1 (t,r_2) dd(r^2)$. The ratio of these is the same at all times since only the overall scale can change
  $
    f_1 (t,r_1) = f_1 (t,r_2) F(r_1,r_2)
  $
  or taking $r_2 = "constant"$ we find
  $
    f_1 (t,r) = f(t) L(r)
  $
  Similarly $f_2$ can be decomposed.
]

By redefining $r$ we can absorb $H$
$
  dd(tau^2) = dd(t^2) - f(t) L(r) dd(r^2) - g(t) r^2 dd(Omega^2)
$
To determine $f, L$ and $g$ we use the EFE. We assume the Universe is a perfect fluid
$
  T_(mu nu) = p g_(mu nu) + (p + rho) U_mu U_nu
$
The comoving frame has $U^i = 0$ and $U^t = 1$. Then $R_(t r) = 0$ implies
$
  dot(g)/g = dot(f)/f
$
so $g = f$ up to a constant. Then
$
  dd(tau^2) = dd(t^2) - f(t) [L(r) dd(r^2)+r^2 dd(Omega^2)]
$
To determine $L(r)$ we consider $tensor(R, r, -r)$. This can not depend on $r$ due to homogeneity meaning the spatial part must be constant
$
  - L'/(r L^2) = 2 k
$
implying
$ L = 1/(tilde(k)-k r^2) $
Similarly for $tensor(R, theta, -theta)$
$
  1/r^2 - 1/(r^2 L) + L'/(2 r L^2) = 2 k
$
Then
$
  (1 - tilde(k))/r^2 + 2k = 2 k
$
so $tilde(k) = 1$. We define the _scale factor_ $a equiv sqrt(f)$ giving the _Friedmann-Robertson-Walker metric_
$
  dd(tau^2) = dd(t^2) - a^2 (t) [1/(1-k r^2) dd(r^2)+r^2 dd(Omega^2)]
$

== Geometric interpretation of $bold(k)$
To see the meaning of $k$ consider the spatial curvature
$
  tensor(R, i, -i) = g^(i k) R_(k i) = - (6 k)/a^2 = - 6 K
$
where $K = k a^(-2)$ is the _Gaussian curvature_. So $k$ determines the spatial curvature!

To determine the shape we construct embedding diagrams. We embed it using
$
  dd(s^2) = (dd(bold(x)))^2 + (dd(z))^2
$
Consider a hypersphere
$
  bold(x)^2 + z^2 = R^2
$
implying
$
  dd((z^2)) = 2 z dd(z) = - dd((bold(x)^2))
$
Then
$
  dd(s^2) & =^("remove" z) (dd(bold(x)))^2 + (dd((bold(x)^2)))^2/(4(R^2-bold(x)^2)) \
  & =^"polar" dd(r^2) + r^2 dd(Omega^2) + (r^2 dd(r^2))/(R^2-r^2) \
  & = 1/(1-r^2\/R^2) dd(r^2) + r^2 dd(Omega^2)
$
defining $rho equiv r\/R$
$
  dd(s^2) = R^2 (dd(rho^2)/(1-rho^2) + rho^2 dd(Omega^2))
$
And if we define $R equiv a(t)$ for some fixed $t$ then this is the FRW metric if $rho^2 = k r^2$ meaning $k > 0$. So for $k > 0$ the FRW describes a sphere! A similar calculation shows that using
$
  - bold(x)^2 + z^2 = R^2
$
then $k <^! 0$. And for $k = 0$ the FRW metric is spatially flat.

== Hubble's law
Consider an observer standing far away from a galaxy. Light has $dd(tau^2)= 0$ (with $dd(Omega^2)=0$) meaning
$
  dd(t^2) = a^2 (t) dd(r^2)/(1-k r^2)
$
implying
$
  integral dd(t)/a(t) = underbracket(-, "incoming" #linebreak() "light") integral dd(r)/sqrt(1-k r^2)
$
We assume the light travels from $(t_1,r_1)$ to $(t_0, 0)$ with $t_0 > t_1$. Then
$
  integral_(t_1)^(t_0) dd(t)/a(t) = - integral_(r_1)^0 dd(r)/sqrt(1-k r^2) = cases(sin^(-1) r_1 #h(2em) &"for" k = +1, r_1 &"for" k = 0, sinh^(-1) r_1 &"for" k=-1)
$
Consider two crests of light emitted at $(t_1, r_1)$ and $(t_1 + dd(t_1, d: delta), r_1)$. The observer would see these at $(t_0, 0)$ and $(t_0+dd(t_0, d: delta), 0)$ respectively. Then
$
  integral_(t_1+dd(t_1, d: delta))^(t_0+dd(t_0, d: delta)) dd(t)/(a(t)) = integral_(t_1)^(t_0) dd(t)/a(t)
$
We assume $dd(t_i, d: delta)$ is very small meaning $a(t_i) tilde "constant"$. Then
$
  (t_0+dd(t_0, d: delta))/a(t_0) - (t_1 + dd(t_1, d: delta))/a(t_1) = t_0/a(t_0) - t_1/a(t_1)
$
implying
$ dd(t_0, d: delta)/a(t_0) = dd(t_1, d: delta)/a(t_1) $
so we have a redshift! With frequencies
$
  lambda_1/lambda_0 = nu_0/nu_1 = dd(t_1, d: delta)/dd(t_0, d: delta) = a(t_1)/a(t_0)
$
We define the redshift $z$ by
$
  z & equiv (lambda_0 - lambda_1)/lambda_1 = a(t_0)/a(t_1) -1 \
    & tilde.eq^"small time" [1 + (dot(a) (t_0))/a(t_0) (t_1-t_0)]^(-1) \
    & tilde.eq^"Taylor" underbracket((dot(a) (t_0))/a(t_0), H_0) (t_0-t_1)
$
with $H_0$ being the _Hubble rate_ at $t_0$.

Consider the _comoving distance_
$
  r equiv integral_(t_1)^(t_0) dd(t)/a(t) tilde.eq (t_0-t_1)/a(t_0)
$
The _proper distance_ at $t_0$ is defined by $d equiv a(t_0) r$ so
$
  z tilde.eq H_0 d
$
Taking the derivative of $d$ we find
$
  v_r = dot(a) (t_0) r = (dot(a) (t_0))/a(t_0) underbracket(a(t_0) r, t_0-t_1) = z
$
This is _Hubble's law_
$
  v_r = H_0 d = z
$
which _proves_ the Universe is expanding!

= The Friedmann equations
== The equations
Consider the EFE
$
  tensor(R, mu, - nu) = 8 pi G (tensor(T, mu, - nu) - 1/2 T tensor(delta, mu, -nu))
$
with
$
  T_(mu nu) =^"perfect fluid" p g_(mu nu) + (p + rho) U_mu U_nu
$
where $U^t = 1$ and $U^i = 0$. We compute
$
  tensor(T, t, -t) = -rho";  " tensor(T, i, -i) = p
$
meaning
$
  T = - rho + 3 p
$
Similarly we compute
$
  tensor(R, t, -t) = 3 dot.double(a)/a";  " tensor(R, r, -r) = dot.double(a)/a +2 dot(a)^2/a^2 + (2 k)/a^2
$
The $t t$-component gives
$
  (dot.double(a))/a = -(4 pi G)/3 (rho + 3 p)
$
this is the _acceleration equation_ or _second Friedmann equation_! The $r r$-component gives
$
  4 pi G (rho-p) & = dot.double(a)/a + 2 (dot(a)^2)/a^2 + (2 k)/a^2
$
We rewrite this using the acceleration equation
$
  H^2 equiv (dot(a)^2)/a^2 & = (8 pi G rho)/3 - k/a^2
$
which is the _first Friedmann equation_! We have defined the _Hubble rate_ $H = dot(a)\/a$. These are very important and can be solved to find $a(t)$. Then solving these give the metric meaning we do not need to solve the full EFE!

We can also use $D_mu T^(mu nu) = 0$ to find
$
        dot(rho) + 3 dot(a)/a (rho + p) & = 0 \
  a^3 dot(rho) + 3 a^2 dot(a) (rho + p) & = 0 \
                      dv(, t) (rho a^3) & = -3 p a^2 dot(a) \
                      dv(, a) (rho a^3) & = - 3 p a^2
$
which is the _fluid equation_.

We also need an equation of state $p = p(rho)$. We assume $p = w rho$ as typical for an _ideal fluid_. Then by the above
$
  dot(rho)/rho & = - 3 (1+p/rho) dot(a)/a \
               & = -3 (1+ w) dot(a)/a
$
implying
$
  rho prop a^(-3 (1+ w))
$

1. $p << rho$ (pressure-less matter). Then $w tilde 0$ giving
$
  rho_m prop 1/a^3
$
this is also called non-relativistic matter.

2. $p = 1/3 rho$ (radiation). Then $w = 1/3$ giving
$
  rho_r prop 1/a^4
$
this is also called relativistic matter.

3. $p = - rho$. Then $w = -1$ giving
$
  rho_Lambda prop "constant"
$
as we would expect for a cosmological constant.

== Consequences
Consider $p >= 0$. Then
$
  dot.double(a)/a = - (4 pi G)/3 (rho+ 3p) < 0";  " dot(a)/a >^"today" 0
$
This implies that $a = 0$ at some point. We define the beginning of the Univeres by $a(t=0) =^! 0$. This is called the _Big Bang_.

The age of the Universe can be computed by
$
  t_0 = integral dv(t, a) dd(a) = integral dd(a)/dot(a) < integral_0^(a(t_0)) dd(a)/(dot(a) (t_0))
$
We find
$
  t_0 < H_0^(-1)
$
for $p >= 0$. Using modern values we have $t_0 lt.tilde 13.8 "Gyr"$.

Consider $p >= 0$. Then
$
  dv(, a) (rho a^3) = - 3 p a^2
$
implies $rho a^3$ is decreasing. Then for $a -> oo$ we have $rho a^2 -> 0$. We have
$
  dot(a)^2 = (8 pi G)/3 rho a^2 - k
$
so $dot(a)^2 -> - k$ as $a -> oo$. For $k = - 1$ the Universe expands forever. For $k = 0$ the expansion slows down asymptotically. For $k = +1$ we have $dot(a) = 0$ for some $a = a_"max"$ given by
$
  rho a_"max"^2 = 3/(8 pi G)
$
after $a_"max"$ the Universe would begin collapsing.

Assume $k=0$. Then
$
  H^2 = (8 pi G)/3 rho
$
We define the _critical density_ by
$
  rho_c equiv (3 H^2)/(8 pi G)
$
so this is the energy density in the Universe if it were flat. We define the _density parameter_
$
  Omega = rho/rho_c
$
Then assuming $k eq.not 0$ we can write the curvature as
$
  k/a^2 = H^2 (Omega-1)
$
meaning $k = 0$ implies $Omega = 1$.

Using the density parameter we can write the first Friedmann equation as
$
  H^2/H_0^2 = Omega_(m,0)/a^3 + Omega_(r,0)/a^4 + Omega_(Lambda,0) + (1- Omega_0)/a^2
$
As an example of why we would use this form consider a _radiation-dominated_ Universe. Then (assuming $k=0$)
$
  H^2/H_0^2 = Omega_(r,0)/a^4
$
or
$
  a dd(a) prop dd(t)
$
integrating gives
$
  a prop t^(1\/2)
$
Similarly a _matter-dominated_ Universe would have $a prop t^(2\/3)$.

= Big Bang cosmology
== The Big Bang
Above we found
$
       rho_m & prop a^(-3) \
       rho_r & prop a^(-4) \
  rho_Lambda & prop "const"
$
These are different implying there are periods where different $rho_i$ dominate. As an example when at early times $rho_r tilde rho_"tot"$ we say the Universe was radiation-dominated. At some early time $rho_r tilde rho_m$ and at some later time $rho_m tilde rho_Lambda$. We call these _periods of equality_. At very late times we expect a Universe completely dominated by $rho_Lambda$ since all other $rho_i tilde 0$ eventually.

Travelling backwards in time $rho$ clearly increases. This implies the Universe has expanded from an initially very hot and dense state. We call this the _hot Big Bang_. We assume particles were in thermal equilibrium due to inter-particle interactions in the early Universe. Then the number density $n$ of bosons ($+$) and fermions ($-$) would be given by
$
  dd(n) = (4 pi g)/(h^3 c^3) (E^2 dd(E))/(e^(E\/k_B T) plus.minus 1)
$
We compute
$
  n & = (4 pi g)/(h^3 c^3) integral (E^2 dd(E))/(e^(E\/k_B T) plus.minus 1) \
  n & prop cases(3/2 T^3 zeta(3) &"   for fermions", 2 T^3 zeta(3) &"   for bosons")
$
implying
$
  n_b = 4/3 n_f = 2.404 g/(2 pi^2) ((k_B T)/(h c))^3
$
Then in the non-relativistic case with $rho prop n$ we find $T prop a^(-1)$. We compute
$
  rho = integral E dd(n) prop T^4
$
Then in the relativistic case we also find $T prop a^(-1)$. With $a_r prop t^(1\/2)$ we find
$
  k_B T tilde.eq 0.46 E_(p l) (t/t_(p l))^(-1\/2)
$
At these scales general relativity breaks down and the formalism we have developed only makes sense for $t > t_(p l)$. So being pedantic the point $t = t_(p l)$ and not $t = 0$ is the Big Bang.

== Big Bang nucleosynthesis
At temperatures $T > "MeV"$ we had equilibrium
$
  p^+ + n <--> D^+ + dots
$
At temperatures $T < "MeV"$ the formed nuclei could no longer be ripped apart so we had reactions like
$
  p^+ + n --> D^+ --> "He" --> "Li"
$
This is _Big Bang nucleosynthesis_.

== The cosmic microwave background
At temperatures $T gt.tilde 13.6 "eV"$ we had equilibrium
$
  e^- + p^+ <--> H + gamma
$
At temperatures $T lt.tilde "eV"$ free electrons and protons form neutral hydrogen
$
  e^- + p^+ --> H + gamma
$
This is called _recombination_. Before recombination photons scatter with free electrons and protons forming a coupled plasma. At this point the Universe is _opaque_. After recombination photons can move freely since they do not scatter with hydrogen making the Universe _transparent_. At some point we had _last scattering_. These photons constitute the _cosmic microwave background_ or CMB. The temperature of CMB photons today is $T_("mm, CMB") = 2.73 "K"$ found by measurement. We can also compute $T_("CMB")$ since we know $z_"rec"$ and $rho_"rec"$. These values agree to high precision!

To summarize recombination happened at $z_"rec" tilde 1370$ after which photons decoupled leading to last scattering at $z_"LS" tilde 1100$.

== CMB anisotropies
When observing the _CMB_ we observe the _last scattering surface_. The primary observation are temperature anisotropies on the scale of
$
  dd(T, d: delta)/T tilde 10^(-5)
$
These are very small confirming the isotropy and homogeneity of the Universe.



Consider a physical distance $cal(l)$ on the last scattering surface. Then
$
  dd(theta, d: delta) tilde.eq cal(l)/d_A
$
with $d_A$ being the _angular diameter distance_ and $dd(theta, d: delta)$ being the _angular separation_. Using the FRW metric we have
$
  dd(s) = a (t_e) r delta theta tilde.eq cal(l)
$
Then
$
  cal(l) & = a(t_e)/a(t_0) a(t_0) r delta theta \
         & = (a(t_0) r delta theta)/(1+z)
$
We define the _particle horizon_ $d_H$ by
$
  d_H & = a(t_0) underbracket(integral_0^(t_0) dd(t)/(a(t)), "comoving horizon") = a(t_0) r
$
Then
$
  cal(l) = (d_H delta theta)/(1+z)
$
implying
$
  d_A = d_H/(z+1) tilde.eq d_H/z
$
We can compute
$
  d_H tilde.eq^"last scattering" 1.4 times 10^(4) "Mpc"
$
Then
$
  cal(l) tilde.eq (d_H delta theta)/z tilde.eq^"last scattering" 0.22 "Mpc" ((delta theta)/(1 degree))
$
We have observed anisotropies with $delta theta tilde 7 degree$. These correspond to $cal(l)_"LS" tilde 1.6 "Mpc"$ and $cal(l)_"today" tilde 1700 "Mpc" tilde$ size of the Universe. This implies the anisotropies are the initial conditions of the Universe.

As is typical when we have some signal with noise we can do a Fourier transform. When done on a sphere the Fourier basis are the _spherical harmonics_ $Y_l^m$. Then
$
  (delta T)/T (theta,phi) = sum_(l=0)^oo sum_(m=-l)^(m=l) a_(l m) Y_(l m) (theta, phi)
$
We want the _correlation function_
$
  C(theta) = evaluated(expval((delta T)/T (hat(n)_1) (delta T)/T (hat(n)_1)))_(hat(n)_1 dot hat(n)_2 = cos theta)
$
This is the average correlation between all pairs of points seperated by the angle $theta$. Assuming Gaussian fluctuations we have
$
  expval(a_(l m))=0";  " expval(a_(l m)^* a_(l' m')) = underbracket(C_(l), "power") delta_(l l') delta_(m m')
$
Then
$
  C(theta) = 1/(4 pi) sum_(l=0)^oo (2 l + 1) C_l P_l (cos theta)
$
with $P_l$ being the _Legendre polynomials_. This leads to the _power-spectrum_ as in ...

The _multipole_ $l$ is related to length since we can probe smaller angles with higher multipoles. As an example for $l = 2$ (dipole) we split the sky in two and for $l = 1$ (monopole) we simply measure the average $T_"mean"$. This also leads to _cosmic variance_.



$l < 100$ corresponds to large $dd(theta, d: delta)$ and is called the _Sachs-Wolfe plateau_. These scales are larger than the _sound horizon_ at last scattering. We define the sound horizon by
$
  lambda = integral_0^(t_0) c_s/(a(t)) dd(t)
$
Then
$
  alpha tilde.eq lambda/d_A
$
is the corresponding angular separation. Assuming a matter-dominated Universe with $a prop t^(2\/3)$ we find
$
  alpha =^"last scattering" (c_s (1+z_gamma)^(-1\/2))/(c[(1+z_0)^(-1\/2)-(1+z_gamma)^(-1\/2)]) tilde.eq 1 degree
$
which corresponds to $l tilde 200$! We used
$
  c_s tilde.eq sqrt(p/rho) = 1/sqrt(3)
$
as the speed of sound. We see $l tilde 200$ corresponds to the first peak of the power-spectrum. We say this is the _fundamental mode_. This implies that for $l < 200$ the initial density and pressure fluctuations are fixed since they have no way to reach equilibrium. At larger $l$ initial overdensities lead to potential wells with higher $rho_"baryon"$ and $rho_"dm"$. Then gravity tries to pull matter together. But baryons and photons are tightly coupled in the early Universe giving rise to radiation pressure. Then the photons try to escape the potential wells pulling baryons along with them. The dark matter resists the baryons escaping and acts to pull them back. This leads to damped oscillations and the production sound waves. The modes of these oscillations, the standing waves, become peaks in the power-spectrum. At very large $l$ we see exponential damping due to the last scattering surface having a _thickness_.

= Inflation
== Problems with the Big Bang
Big Bang cosmology has multiple problems. The first problem is the _flatness problem_. The model assumes the Universe is flat. We see from the first Friedmann equation
$
  H^2 = (8 pi G)/3 rho - k/a^2
$
that $rho_"curvature" prop a^(-2)$. We should be able to measure this but we simply do not. This is explainable if $rho_"curvature"$ is very small.

The second problem is the _horizon problem_ which is unexplainable. We previously defined the particle horizon
$
  d_H = a(t) integral_0^t dd(t')/(a(t'))
$
We also define the _Hubble radius_ $H^(-1)$. Taking $H^(-1)$ as constant then any light from $d > H^(-1)$ would never reach us since space would recede faster that the speed of light. Then each volume $(H^(-1))^3$ would be _causally disconnected_. Assuming $a prop t^p$ with $0 < p < 1$ then
$
  d_H tilde H^(-1)
$
for a matter-dominated Universe we would have $d_H = 2 H^(-1)$. Consider
$
  lambda_H (t_"LS") & = d_H (t_0) (a_"LS"/a_0) \
                    & = d_H (t_0) (T_0/T_"LS")
$
This is the size of the observable Universe at last scattering if we scaled the current size. Consider also
$
  H_"LS"^(-1) & tilde.eq^"matter-dominated" H_0^(-1) (a_"LS"/a_0)^(3\/2) tilde d_H (t_0) (T_0/T_"LS")^(3/2)
$
This is the Hubble radius at last scattering. We compute
$ lambda_H^3/d_H^3 tilde 10^6 $
so the Universe consisted of $tilde 10^6$ causally disconnected volumes at last scattering! This is a problem since observations of the CMB imply isotropy and homogeneity. There needs to be some mechanism allowing the Universe to _thermalize_. This is where _inflation_ comes in. We can plot $(lambda_H, a)$ as in ...

Assuming the Universe is radiation-dominated we have
$
  H^(-1) prop a^2
$
while
$
  lambda_H prop a
$
by definition. We see today $lambda_H < H^(-1)$ but at early times $lambda_H > H^(-1)$ which is the problem described above. We need $lambda_H < H^(-1)$ at some even earlier time. To do this inflation assumes a period with $H^(-1) tilde.eq "constant"$ and $a prop e^(H t)$. This is the same as requiring a period with
$
  0 < dv(, t) (lambda/(abs(H)^(-1))) = dv(, t) (a abs(dot(a)/a)) = dv(, t) abs(dot(a))
$
Assuming $dot(a) > 0$ this implies $dot.double(a) > 0$.

== Slow-roll inflation
Consider a scalar field described by the action
$
  S & = integral dd(x, 4) overbracket(sqrt(-g), "minimal coupling") cal(L) \
    & = integral dd(x, 4) sqrt(-g) [1/2 partial_mu phi partial^mu phi + V(phi)]
$
Then solving the Euler-Lagrange equations give
$
  dot.double(phi) + underbracket(3 H dot(phi), "friction") - (nabla^2 phi)/a^2 + V_phi (phi) = 0
$
Assuming homogeneity the energy-momentum tensor is
$
  T_(mu nu) = partial_mu phi partial_nu phi - g_(mu nu) cal(L)
$
implying
$
  rho = dot(phi)^2/2 + V(phi)";  " p= dot(phi)^2/2 - V(phi)
$
When $V(phi) >> dot(phi)^2$ then $p = - rho$ which is just the condition for a cosmological constant. This implies an accelerated expansion! We define the number of _$e$-folds_ by $a(t_R) = a(t_i) e^N$ or
$
  underbracket(N, "number of "e"-folds") = ln((a(t_R))/(a(t_i))) gt.tilde underbracket(70, "length of inflation")
$
with $t_R$ being the time at _reheating_. The $tilde 70$ comes from requiring $lambda_H (t_i) lt.tilde H_"inf"^(-1)$
$
  lambda_H (t_i) & = d_H (t_0) (a_i/a_R) (a_R/a_0) \
                 & tilde.eq H_0^(-1) e^(-N) T_0/T_R lt.tilde H_"inf"^(-1)
$
so
$
  e^N gt.tilde T_0/H_0 H_"inf"/T_R => N gt.tilde ln(T_0/H_0) - ln(T_R/H_"inf") tilde 70
$

This restricts $V(phi)$ since $phi$ has to be _stuck for long enough_ determined by $+3 H dot(phi)$.

With _slow-roll inflation_ we assume $dot(phi)^2 << V(phi)$ and small $dot.double(phi)$. This is enough for the above assumptions to hold. Then
$
  underbracket(H^2 tilde.eq (8 pi G)/3 V(phi), "first Friedmann equation")";  " 3 H dot(phi) = - V_phi
$
which are the _slow-roll equations_. We define the _slow-roll parameters_
$
  epsilon.alt = 4 pi G dot(phi)^2/H^2";  " eta = 1/(8 pi G) (V''/V)
$
both are $<< 1$ since
$
  eta - epsilon.alt = - dot.double(phi)/(H phi)
$
Then
$
  N & = integral_(t_i)^(t_R) H dd(t) \
    & tilde 8 pi G integral_(phi_R)^(phi_i) V/V_phi dd(phi) gt.tilde^! 70
$
which is now a condition on $V(phi)$!

From the slow-roll equations we see $H tilde$ constant if $V(phi) tilde$ constant. When $V(phi)$ eventually begins _rolling quickly_, $H$ stops being constant, and inflation stops. This is quantified by $epsilon.alt tilde 1$ or equivalently $dot(phi)^2 tilde V(phi)$. When this happens $p tilde 0$ and $rho > 0$ implying $dot.double(a) < 0$ by the fluid equation! While $V(phi)$ is _rolling slowly_ the _effective cosmological constant_ $Lambda_"eff"$ dominates leading to a large repulsive pressure causing rapid expansion with $a prop e^(H t)$. After slow-roll inflation we have $Lambda_"eff" -> Lambda_"vacuum"$ with $Lambda_"vacuum"$ being dominated by $rho_r$ and $rho_m$.

== Quantum fluctuations
Assuming the field is inhomogeneous we have
$
  phi = phi_c (t) + underbracket(dd(phi(t, bold(x)), d: delta), "small")
$
We can find an equation for $dd(phi, d: delta)$ by
$
  dd(phi_k, d: delta) = integral dd(x, 3)/(2 pi)^(3\/2) dd(phi, d: delta) e^(-i bold(k) dot bold(x))
$
We find
$
  dd(dot.double(phi)_k, d: delta) + 3 H dd(dot(phi)_k, d: delta) + k^2/a^2 dd(phi_k, d: delta) + V_(phi phi) dd(phi_k, d: delta) = 0
$
Assuming $k < a H$ then
$
  dd(phi_k, d: delta) tilde H/(2 pi) (k/k_*)^(n_s-1)";  " n_s - 1 = -6 epsilon.alt + 2 eta
$
We can show $dd(phi_k, d: delta)$ gives rise to $dd(T, d: delta)\/T$. So they determine the initial conditions of the Universe! Using the CMB for $l < 100$ we can measure these perturbations.

= Gravitational waves
== Linearized field equations
We write the metric as
$
  g_(mu nu) = eta_(mu nu) + h_(mu nu)
$
with $abs(h_(mu nu)) << 1$. The Riemann tensor is
$
  R_(alpha mu beta nu) prop partial^2 g underbracket(+ dots + Gamma^2, "vanishing")
$
with the linearized metric we find (trivially $R^((0)) = 0$)
$
  R^((1))_(alpha mu beta nu) = 1/2 (partial_alpha partial_nu h_(mu beta) + partial_mu partial_beta h_(alpha nu) - partial_alpha partial_beta h_(mu nu) - partial_mu partial_nu h_(alpha beta))
$
since $partial eta = 0$. Then the Ricci tensor becomes
$
  R^((1))_(mu nu) = eta^(alpha beta) R^((1))_(alpha mu beta nu) = 1/2 ( partial_alpha partial_nu tensor(h, alpha, -mu) + partial_mu partial_alpha tensor(h, alpha, -nu) - partial^2 h_(mu nu) - partial_mu partial_nu h )
$
with $h$ being the trace of $h_(mu nu)$. Then the Ricci scalar becomes
$
  R^((1)) = eta^(mu nu) R_(mu nu)^((1)) = partial_mu partial_nu h^(mu nu) - partial^2 h
$
The Einstein field equations become
$
  R_(mu nu)^((1)) - 1/2 eta_(mu nu) R^((1)) = - 8 pi G T_(mu nu)^((0))
$
we take the trace
$
  R^((1)) = -8 pi G tensor(T^((0)), mu, -mu)
$
and we find
$
  partial_mu partial_nu h^(mu nu)-partial^2 h = 8 pi G tensor(T^((0)), mu, -mu)
$
which looks like a wave equation!

== Gauge freedom
The above $h$ has $10$ degrees of freedom. We want to reduce this number.

Consider
$
  x^mu -> x'^mu = x^mu + xi^mu (x)
$
We know how the metric transforms
$
  g'_(alpha beta) &= pdv(x^mu, x'^alpha) pdv(x^nu, x'^beta) g_(mu nu) \ &= pdv(x^mu, x'^alpha) pdv(x^nu, x'^beta) (eta_(mu nu) + h_(mu nu)) \
  &= (tensor(delta, mu, -alpha) tensor(delta, nu, -beta) - tensor(delta, mu, -alpha) partial_beta xi^nu - partial_alpha xi^mu tensor(delta, nu, -beta) + cal(O) (xi^2)) (eta_(mu nu) + h_(mu nu)) \
  &= eta_(alpha beta) + h_(alpha beta) - partial_beta xi_alpha - partial_alpha xi_beta
$
we see
$
  h_(alpha beta) -> h'_(alpha beta) = h_(alpha beta) - partial_alpha xi_beta - partial_beta xi_alpha
$
This gives us gauge freedom. We can pick
$
  partial_mu h^(mu nu) = 1/2 partial^nu h
$
we can do this since we have four degrees of freedom in $xi^mu$. This is called the _Lorentz gauge_ or the _de Donder gauge_ if written in the form
$
  g^(mu nu) tensor(Gamma, lambda, -mu nu) = 0
$
Using this we have in the Lorentz gauge
$
  partial^2 h = - 16 pi G tensor(T^((0)), mu, -mu)
$
We define
$
  overline(h)_(mu nu) equiv h_(mu nu) - 1/2 eta_(mu nu) h
$
then
$
  partial^mu overline(h)_(mu nu) = 0
$
giving
$
  partial^2 overline(h)_(mu nu) = - 16 pi G T^((0))_(mu nu)
$
Consider
$
  h -> h' = h - 2 partial_mu xi^mu
$
then
$
  partial^mu h'_(mu nu) &= partial^mu h_(mu nu) - partial^2 xi_nu - partial_mu partial_nu xi^mu \
  &=^"Lorentz" 1/2 partial_nu h' = 1/2 partial_nu h - 2 partial_nu partial_mu xi^mu
$
so
$
  1/2 partial_nu h = partial^mu h_(mu nu) underbracket(- partial^2 xi_nu, =^! 0)
$
meaning if $partial^2 xi_nu = 0$ the Lorentz gauge is satisfied! This is called residual gauge degree of freedoms. And we have now fixed $8$ out of $10$ degrees of freedom! So we are left with two dynamical degrees of freedom.

== Plane waves in vacuum
We consider
$
  T_(mu nu) = 0
$
meaning
$
  partial^2 h = 0";  " partial^2 overline(h)_(mu nu) = 0
$
implying
$
  partial^2 h_(mu nu) = 0
$
with the condition
$
  partial_mu h^(mu nu) = 1/2 partial^nu h
$
We make the ansatz
$
  h_(mu nu) (lambda) = underbracket(epsilon_(mu nu), "polarization") e^(i k_alpha x^alpha)
$
and $k^alpha = vecrow(omega, bold(k))$. We find
$
  k^2 epsilon_(mu nu) e^(i k_alpha x^alpha) = 0
$
meaning
$
  k^2 = - omega^2 + bold(k)^2 = 0
$
So the ansatz works given the normal dispersion relation for a plane wave holds. We also find
$
  partial_mu h^(mu nu) = 0 => k^mu epsilon_(mu nu) = 0
$
which is called the _transverse gauge condition_.

The condition $k^mu epsilon_(mu nu) = 0$ still has degrees of freedom we need to fix. We impose
$
  epsilon_(mu 0) = epsilon_(0 mu) = 0
$
and
$
  tensor(epsilon, -mu, mu) = 0
$
As an example take
$
  k^alpha = vecrow(omega, 0, 0, underbracket(k_z, =^! omega))
$
then
$
  k^mu epsilon_(mu nu) = 0 => omega epsilon_(3 nu) = 0
$
giving (with $epsilon_(nu 0) = 0$ and $tensor(epsilon, mu, -mu) = 0$ and $epsilon_(mu nu) = epsilon_(nu mu)$)
$
  h_(mu nu) (z,t) = mat(0, 0, 0, 0; 0, h_+, h_times, 0; 0, h_times, - h_+, 0; 0, 0, 0, 0) e^(i omega (z-t))
$
We can define the basis
$
  epsilon_+^(mu nu) = h_+ mat(0, 0, 0, 0; 0, 1, 0, 0; 0, 0, -1, 0; 0, 0, 0, 0)";  " epsilon_times^(mu nu) = h_times mat(0, 0, 0, 0; 0, 0, 1, 0; 0, 1, 0, 0; 0, 0, 0, 0)
$
Consider a test particle at rest with $U^mu = vecrow(1, 0, 0, 0)$. The geodesic equation is
$
  dv(U^mu, tau) + tensor(Gamma, mu, -nu lambda) U^nu U^lambda = 0 => dv(U^mu, tau) = - tensor(Gamma, mu, -00)
$
where
$
  tensor(Gamma, mu, -00) = 1/2 eta^(mu nu) (partial_0 h_(nu 0) + partial_0 h_(0 nu) - partial_nu h_00) = 0
$
So the effect on a test particle is
$
  dv(U^mu, tau) = 0
$
meaning there is no effect!

Consider two test particles with $x^mu = vecrow(0, 0, 0, 0)$ and $x^mu + dd(x^mu) = vecrow(0, xi, 0, 0)$. Then
$
  dd(s) & = sqrt(g_(mu nu) dd(x^mu, x^nu)) \
        & = sqrt(g_11) xi \
        & =^"binom" (eta_11 + 1/2 h_11) xi
$
for $+$ polarization we find
$
  dd(s) & = (1 + 1/2 h_+ sin(t-z)) xi
$
for separation along $y$ we find
$
  dd(s) & = (1 - 1/2 h_+ sin(t-z)) xi
$
So we find a stretch along $x$ and a squeeze along $y$ (initially) with it oscillating as one would expect.

For $times$ polarization one would find a _diagonal_ stretching and squeezing.

