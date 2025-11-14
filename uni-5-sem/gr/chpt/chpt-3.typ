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
