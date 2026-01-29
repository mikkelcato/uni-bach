//**** init-ting
#import "@preview/physica:0.9.8": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Foundational observations
== The sky is dark
Let $n_*$ be the number density of stars in the Universe and $R_*$ the typical size of a star. Then consider a cylinder of radius $R_*$ about our line of sight. Given the center of a star is within the cylinder then our view will be blocked. Let the length of the cylinder be $lambda$ then it has volume $V = lambda pi R_*^2$. The average number of stars blocking our view can then be computed as
$ N = n_* V = n_* lambda pi R_*^2 $
The typical distance before our view is blocked will be the $lambda$ for which $N = 1$ meaning
$ lambda = 1/(n_* pi R_*^2) < oo $
implying the sky should be full of stars! To determine the brightness of the sky we consider a star of size $R_*$ with luminosity $L_*$ at a distance $lambda >> R_*$. This star takes up an area
$ Omega = (pi R_*^2 )/ lambda^2 $
The flux we would measure is then
$ f = L_* /(4 pi lambda^2) $
Giving the brightness
$ Sigma_* = f/Omega = L_* /(4 pi^2 R_*^2) $
This does not depend on $lambda$! So the brightness of the sky will be equal to the brightness of a typical star, i.e. the Sun. This is referred to as Olbers' paradox since we observe the sky is dark.

Later we will see that the Universe is finite which resolves the paradox.

== The cosmological principle
At very large scales $(~ 100 "Mpc")$ the Universe is isotropic and homogeneous. This is the cosmological principle. We believe this since the Universe appears isotropic around us and we have no reason to think our place in the Universe is special (the Copernican principle).

== Hubble's law
When observing galaxies we see a redshift $z$ defined by
$
  z equiv (lambda_0 - lambda_e)/(lambda_e)
$
We interpret this as a Doppler shift implying
$
  v = c z
$
Since most galaxies have $z > 0$ everything appears to move away from us. Hubble showed by observation that
$
  c z tilde.eq H_0 r tilde.eq v
$
This is Hubble's law.

This result is expected given the Universe is expanding homogeneously and isotropically. To see this consider three galaxies and let $r_(i j) equiv abs(bold(r)_i - bold(r)_j)$. Assuming the Universe expands in the above manner, then the shape of the triangle they form must be preserved. This implies
$ r_(i j) (t) = a(t) r_(i j) (t_0) $
with $a(t)$ being the scale factor. Then
$
  v_(i j) (t) = dot(a) r_(i j) (t_0) = dot(a)/a r_(i j) (t) equiv H r_(i j) (t)
$
where we define the Hubble parameter
$
  H equiv dot(a)/a
$
This implies a Big Bang since if everything is moving away from everything then at some point everything was clumped. We can actually estimate the age of the Universe by the Hubble time
$
  t_0 = r/v = H_0^(-1) tilde 14 "Gyr"
$

= The FRW cosmology
== The FRW metric
We assume the Universe is described by the FRW metric
$
  dd(s)^2 = -c^2 dd(t)^2 + a(t)^2 (dd(r)^2/(1-k r^2) + r^2 dd(Omega)^2)
$
and typically we assume a flat Universe ($k = 0$).

We define the proper distance by
$
  d_p (t) = a(t) underbracket(integral_0^r dd(r')/sqrt(1-k r'^2), "comoving distance") = a(t) chi
$
With a redefinition of $r$ (where $chi ->^tilde r$) we can write
$
  d_p (t) & = a(t) r \
  v_p (t) & = H d_p (t)
$
so $d_p > d_H equiv c H_0^(-1)$ implies $v_p > c$.

Under this redefinition the metric becomes
$
  dd(s)^2 = - c^2 dd(t^2) + a(t)^2 (dd(r^2) + S_k^2 (r) dd(Omega^2))
$
where the dependence on $k$ is hidden in $S_k$. Then
$
  c^2 dd(t)^2 =^"for photons" a(t)^2 dd(r)^2
$
implying
$
  (c dd(t))/a = dd(r)
$
We integrate to find
$
  r & = c integral_(t_e)^t_0 dd(t)/a(t) \
  r & =^"and" c integral_(t_e+lambda_e \/c)^(t_0 + lambda_0\/c) dd(t)/a(t)
$
Then
$
  integral_(t_e)^(t_e + lambda_e\/c) dd(t)/a(t) &= integral_(t_0)^(t_0 + lambda_0 \/c) dd(t) /a(t)
$
Assuming $dd(t, d: Delta) tilde 0$ we have $a(t) tilde "constant"$ which eventually gives
$
  1 + z = a(t_0)/a(t_e) =^"convention" 1/a(t_e)
$
This is nice since we can measure $z$.

== The FRW dynamics
The first equation we need is the Friedmann equation
$
  H^2 = (8 pi G)/(3 c^2) epsilon - (k c^2)/R_0^2 1/a^2
$
which relates $epsilon, k$ and $a$. Defining
$
  epsilon_c equiv^(k=^!0) (3 c^2 H^2)/(8 pi G)";  " Omega equiv (epsilon)/(epsilon_c)
$
We have
$
  k/R_0^2 = H_0^2/c^2 underbracket((Omega_0 -1), tilde 0)
$

The second equation we need is the fluid equation
$
  dot(epsilon) + (3dot(a))/a (epsilon + P) = 0
$
which does not depend on $k$. We can combine this with the Friedmann equation to obtain the acceleration equation
$
  dot.double(a)/a & = - (4 pi G)/(3 c^2) (epsilon + 3 P)
$
As the name implies this tells us whether the expansion of the Universe is speeding up or slowing down.

We also need an equation of state $P = P(epsilon)$. Typically  we assume everything is a perfect fluid with $P = w epsilon$. For non-relativistic matter we have $w tilde 0$ while relativistic matter has $w tilde 1/3$. Aside from these observations hint at a component with constant $epsilon$. This is referred to as dark energy or the cosmological constant $Lambda$. With the addition of $Lambda$ we have
$
  H^2 &= (8 pi G)/(3 c^2) epsilon - (kappa c^2)/R_0^2 1/a^2 + Lambda/3 \
  dot.double(a)/a &= - (4 pi G)/(3 c^2) (epsilon + 3 P) + Lambda/3
$
Then $Lambda$ corresponds to a density $epsilon_Lambda$ with
$
  epsilon_Lambda equiv (c^2 Lambda)/(8 pi G)
$
For $Lambda =^! "constant"$ we require $P_Lambda = - epsilon_Lambda$ by the fluid equation implying $w = -1$.

== Solving the fluid equation
The $epsilon_i$ and $P_i$ are additive. We assume different components do not interact, meaning each component satisfies the fluid equation
$
  dot(epsilon) +3 dot(a)/a (1 + w) epsilon = 0
$
Which we can solve to find
$
  epsilon = epsilon_0 a^(-3 (1 + w))
$
We find
$
       epsilon_m & = epsilon_(m,0)/a^(3) \
       epsilon_r & = epsilon_(r,0)/a^(4) \
  epsilon_Lambda & = epsilon_(Lambda,0)
$
Today $epsilon_Lambda$ dominates, but due to the different scalings the Universe has had different periods of domination. Consider
$
  epsilon_Lambda/epsilon_m = epsilon_(Lambda, 0)/(epsilon_(m,0)) a^3 = Omega_(Lambda,0)/Omega_(m,0) a^3
$
The scale factor at matter-dark energy equality $a_(m Lambda)$ can be found by requiring $epsilon_Lambda = epsilon_m$
$
  a_(m Lambda) = (Omega_(m,0)/Omega_(Lambda,0))^(1\/3)
$
Similarly the scale factor at radiation-matter equality $a_(r m)$ is
$
  a_(r m) = (Omega_(r,0))/(Omega_(m,0))
$
Since we have now derived how $epsilon_i$ scales we can write the Friedmann equation as
$
  H^2/H_0^2 = Omega_(m,0)/a^3 + Omega_(r,0)/a^4 + Omega_(Lambda,0) + Omega_(k,0)/a^2
$
which is simple to solve. And everything we would care about basically follows from this. We typically use the _Benchmark Model_ where
$
  Omega_0 = Omega_(m,0) + Omega_(r,0) + Omega_(Lambda,0) tilde.eq 1";  " Omega_(k,0) = 1 - Omega_0
$

== The empty universe
As an example of a model universe we consider an empty universe. The Friedmann equation becomes
$
  dot(a)^2 = - (k c^2)/R_0^2
$
We find two solutions. The first has $dot(a) = 0$ meaning $k = 0$. This describes an empty, static and flat universe. The second has $k = -1$ (a Milne universe) with
$
  dot(a) = plus.minus c/R_0
$
This implies
$
  a = t/t_0";  " t_0 = R_0/c = H_0^(-1)
$
Then light we observe at $t_0$ is emitted at
$
  t_e = (H_0^(-1))/(1 + z)
$
Similarly
$
  d_p (t_0) & = c integral_(t_e)^(t_0) dd(t)/a(t) \
            & = c/H_0 ln(1+z)
$
Using
$
  (d_p (t_e))/a(t_e) = (d_p (t_0))/a(t_0)
$
we find
$
  d_p (t_e) = c/H_0 ln(1+z)/(1+z)
$
For any other model universe everything follows as above, the results are just more complicated. As an example one would consider a three-component universe with radiation, matter and dark energy to compute the age of the Universe.


== Finding $a(t)$
We Taylor expand the scale factor as
$
  a(t) &= a(t_0) + dv(a, t)_(t=t_0) (t-t_0) + 1/2 dv(a, t, 2)_(t_t_0) (t-t_0)^2 + dots \
  & approx 1 + H_0 (t-t_0) - 1/2 q_0 H_0^2 (t-t_0)^2
$
with $q_0$ being the deceleration parameter
$
  q_0 equiv - ((dot.double(a) a)/dot(a)^2)_(t=t_0) = - (dot.double(a)/(a H^2))_(t=t_0)
$
Then knowing $H_0$ and $q_0$ we can approximate $a(t)$.

Using the acceleration equation we have
$
  -dot.double(a)/(a H^2) &= 1/2 ((8 pi G)/(3 c^2 H^2)) sum_(i=1)^N epsilon_i (1 + 3 w_i) \
  &= 1/2 sum_(i=1)^N Omega_i (1 + 3 w_i)
$
Then
$
  q_0 & = 1/2 sum_(i=1)^N Omega_(i,0) (1 + 3 w_i) \
  & = Omega_(r,0) + 1/2 Omega_(m,0) - Omega_(Lambda,0) approx^"benchmark" -0.53
$
For $H_0$ we consider
$
  1/a(t) tilde.eq 1 - H_0 (t-t_0) + (1 + q_0/2) H_0^2 (t-t_0)^2
$
Then
$
  d_p (t_0) tilde.eq c (t_0-t_e) + (c H_0)/2 (t_0 - t_e)^2
$
Using
$
  z & = 1/a(t_e) -1 \
    & tilde.eq H_0 (t_0 - t_e) + (1+q_0/2) H_0^2 (t_0 - t_e)^2
$
we find
$
  t_0-t_e tilde.eq H_0^(-1) [z-(1+q_0/2) z^2]
$
Then
$
  d_p (t_0) tilde.eq (c z)/H_0 [1-(1+ q_0)/2 z]
$
So we just need $d_p (t_0)$. We define the luminosity distance
$
  d_L equiv (L/(4 pi f))^(1\/2)
$
With the FRW metric the area of a sphere is
$
  A_p (t_0) = 4 pi S_kappa (r)^2
$
The flux will be decreased by a factor $(1+z)^(-2)$, meaning we have
$
  f = L/(4 pi S_kappa (r)^2 (1+z)^2)
$
Then
$
  d_L & = S_kappa (r) (1+z) \
      & =^(k=0) r (1+ z) \
      & = d_p (t_0) (1+z)
$
Using this we find
$
  d_L tilde.eq c/H_0 z (1 + (1-q_0)/2 z)
$
Which can be used to measure $H_0$.

We also define the angular-diameter distance
$
  d_A equiv cal(l)/dd(theta, d: delta)
$
with $cal(l)$ being the size of what we observe and $dd(theta, d: delta)$ being the angle it subtends. Then
$
  cal(l) & = dd(s) \
         & = a(t_e) S_k (r) dd(theta, d: delta)
$
We find
$
  d_A & = (S_k (r))/(1+z) \
      & = d_L/(1+z)^2 \
      & =^(k = 0) d_p (t_e)
$
Then we have
$
  d_p (t_0) = d_A (1+z) = d_L/(1+z)
$
Taking $z -> oo$ where $t_e -> 0$ meaning $d_p (t_0) -> d_"hor" (t_0)$ we find
$
  d_L tilde.eq z d_"hor" (t_0)\
  d_A tilde.eq (d_"hor" (t_0))/z
$

= The cosmic microwave background
== The observation
We observe the sky being uniformly bright with a temperature $expval(T)_0 = 2.7255 "K" equiv T_0$. We call the radiation giving rise to this temperature the cosmic microwave background or CMB. The CMB is famously a blackbody up to very small fluctuations. These are characterized by
$
  dd(T, d: delta)/T (theta,phi) equiv (T-expval(T))/expval(T)
$
with
$
  expval((dd(T, d: delta)/T)^2)^(1/2) tilde 10^(-4)
$

For future reference we note $T prop a^(-1)$ and define the baryon-photon ratio
$
  eta equiv n_("bary",0)/n_(gamma,0) tilde 6.1 times 10^(-10)
$

== The physics
We assume all baryons consist of neutral hydrogen $H$ and protons $p$. By charge neutrality we require $n_p = n_e$. We define
$
  X equiv n_p/(n_p+n_H) = n_p/n_"bary" = n_e/n_"bary"
$
We would like to determine $X$. We care about
$
  H + gamma ->^"ionization"_(hbar omega > Q) p + e^- ->^"recombination" H + gamma
$
At $a tilde 10^(-5)$ any hydrogen would be ionized instantly meaning $X = 1$ since $h f_"mean" tilde 60 "eV"$. At this time the Universe was opaque since $gamma$ scatter with $e^-$ through
$
  gamma + e^- -> gamma + e^-
$
This happens at a rate
$
  Gamma & = c/lambda \
        & = n_e sigma_e c \
        & =^(X=1) (n_("bary",0) sigma_e c)/a^3 tilde 10^(-6) "s"^(-1)
$
As long as $Gamma > H$ we say $gamma$ and $e^-$ are coupled. When $Gamma < H$ they decouple and the Universe becomes transparent. For $a < a_(r m) tilde 10^(-4)$ we have
$
  H = (H_0 Omega_(r,0)^(1\/2))/a^2 tilde 10^(-10) "s"^(-1)
$
So $gamma$ and $e^-$ were strongly coupled.

For recombination we expect
$
  T_"rec" tilde Q/(2.7 k_B) tilde 6 times 10^4 "K"
$
which is very crude. We would like to better understand
$
  H + gamma harpoons.rtlb p + e^-
$
Assuming thermal equilibrium the number density $n_x$ of a species $x$ is determined by
$
  n_x (p) dd(p) = g_x (4 pi)/h^3 (p^2 dd(p))/(exp((E-mu_x)/(k_B T))plus.minus 1)
$
For $gamma$ with $E = p c = h f$ and $mu_gamma = 0$ we have
$
  n_gamma (f) dd(f) = (8pi)/c^3 (f^2 dd(f))/(exp((h f)/(k_B T))-1)
$
Then
$
  n_gamma = (2.4041)/pi^2 ((k_B T)/(hbar c))^3
$
Similarly with $E = p c$ and $mu_x << E$ we have
$
  epsilon_x &= integral_0^oo dd(p) n_x (p) E \
  &= g_x (k_B^4 T^4)/(hbar^3 c^3 2 pi^2) integral_0^oo dd(y) y^3/(exp(y) plus.minus 1) \
  &= g_x pi^2/30 (k_B^4 T^4)/(hbar^3 c^3) cases(1 ": boson", 7/8 ": fermion")
$
This implies $T prop a^(-1)$. Since different species stop being relativistic over time we write
$
  epsilon_r (t) = sum_("ultra rel." x) epsilon_x (t) equiv g_* (t) pi^2/30 (k_B^4 T^4)/(hbar^3 c^3)
$
where $g_*$ is defined by
$
  g_* (t) equiv sum_("rel. bosons" x) g_x + 7/8 sum_("rel. fermions" y) g_y
$
This is the effective number of relativistic bosonic degrees of freedom.

At recombination $e^-, p "and" H$ were non-relativistic. Then we have $m c^2 >> k_B T$ and $p tilde m_x v$ with
$
  E tilde m_x c^2 + p^2/(2 m_x)
$
We find
$
  n_x (p) dd(p) = g_x (4 pi)/h^3 exp((-m_x c^2 + mu_x)/(k_B T)) exp(- p^2/(2 m_x k_B T)) p^2 dd(p)
$
Then
$
  n_x = g_x ((m_x k_B T)/(2 pi hbar^2))^(3\/2) exp((-m_x c^2 + mu_x)/(k_B T))
$
and
$
  epsilon_x tilde.eq m_x c^2 n_x prop exp(-(m_x c^2 - mu_x)/(k_B T)) << 1
$
Which is consistent with the Universe being radiation dominated at recombination.

We observe $n_overline(x) << n_x$ since we have matter. This implies $mu_x eq.not mu_overline(x)$, otherwise both $n_i$ would be exponentially supressed.

At recombination we assume equilibrium meaning (since $mu_gamma = 0$)
$
  mu_H = mu_p + mu_e
$
This implies
$
  n_H/(n_p n_e) = g_H/(g_p g_e) (m_H/(m_p m_e))^(3\/2) ((k_B T)/(2 pi hbar^2))^(-3\/2) exp[((m_p + m_e - m_H)c^2)/(k_B T)]
$
We have $1 tilde.eq m_H\/m_p$ and $Q equiv (m_p+m_e-m_H)c^2$ and $g_p=g_e=2$ and $g_H = 4$. Then
$
  n_H/(n_p n_e) = ((m_e k_B T)/(2 pi hbar^2))^(-3\/2) exp(Q/(k_B T))
$
which is the Saha equation. By definition
$
  n_H = (1-X)/X n_p
$
Then using $n_e = n_p$ and
$ eta equiv n_"bary"/n_gamma = n_p/(X n_gamma) $
we find
$
  n_p = 0.2436 X eta ((k_B T)/(hbar c))^3
$
Then the Saha equation becomes
$
  (1-X)/X^2 = 3.84 eta ((k_B T)/(m_e c^2))^(3\/2) exp(Q/(k_B T))
$
The solution is
$
  X = (-1 + sqrt(1+4S))/(2 S)";  " S(T,eta) = 3.84 eta ((k_B T)/(m_e c^2))^(3\/2) exp(Q/(k_B T))
$
Defining recombination by $X = 1/2$ we find
$
  k_B T_"rec" = 0.324 "eV" = Q/42
$
so
$
  T_"rec" = 3760 "K" tilde 2.5 times 10^5 "yr"
$

To find the time of decoupling consider
$
  Gamma (z) & = n_e (z) sigma_e c \
            & = X(z) (1+z)^3 n_("bary",0) sigma_e c
$
Assuming the Universe was matter-dominated we find
$
  H = H_0 Omega_(m,0)^(1/2) (1+z)^(3/2)
$
At $Gamma = H$
$
  1+z_"dec" = 39.3/(X z_"dec"^(2\/3))
$
giving $z_"dec" = 1120$. A more complicated computation gives $z_"dec" tilde.eq 1090$.

The last thing we need is the time of last scattering $z_"ls"$. This differs slightly from $z_"ls"$. We consider the optical depth
$
  tau(t) & = integral_t^(t_0) Gamma(t) dd(t)
$
This is the average number of scatterings a $gamma$ from the CMB has undergone since $t$. We can rewrite this as
$
  tau(a) & = integral_a^1 Gamma(a)/dot(a) dd(a) \
         & = integral_a^1 Gamma(a)/H(a) dd(a)/a \
  tau(z) & = integral_0^z Gamma(z)/H(z) dd(z)/(1+z) \
         & prop integral_0^z X(z) sqrt(1+z) dd(z) =^"last scattering" 1
$
This is annoying to solve so we just take $z_"ls" tilde z_"dec" tilde 1090$.

== Nucleosynthesis briefly
The very early universe was radiation dominated meaning $a(t) prop t^(1\/2)$ so
$
  T(t) prop t^(-1\/2)
$
By $E_"mean" (t) tilde.eq 2.7 k_B T(t)$ we can then find the energy. We can probe energies corresponding to $tilde 10^(-13) "s"$. As the Universe has expanded the energy has dropped from
$
  E_P tilde 10^(28)"eV" -> E_"mean" (t_0) tilde 10^(-3)"eV"
$
Above we say that the CMB physics happened at $E tilde 10 "eV"$. We now consider nucleosynthesis in the early Universe which happened at $E tilde "MeV"$.

We define the mass number $A equiv Z + N$ and the binding energy $B$ which is the energy required to pull a nucleus apart into its component protons and neutrons. As an example
$
  p + n harpoons.rtlb "D" + underbracket(2.22 "MeV", B_D)
$
We define
$
  Q_n equiv (m_n - m_p )c^2 = 1.29 "MeV"
$
All nucleosynthesis starts with $n$ and $p$ fusing to form heavier elements. This process is naturally limited since any free neutron is unstable and will decay
$
  n -> p + e^- + overline(nu)_e
$
with $tau_n = 880 "s"$.

When $E_"mean" tilde 10 "MeV"$ at $t = 0.1"s"$ we had
$
  gamma + gamma harpoons.rtlb e^- + e^+
$
and $n$ and $p$ were coupled by $nu_e$
$
  n + nu_e & harpoons.rtlb p + e^- \
   n + e^+ & harpoons.rtlb p + overline(nu)_e
$
With $k_B T << m_p c^2$ and assuming thermal equilibrium we find
$
  n_n/n_p tilde.eq exp(- Q_n/(k_B T))
$
At some time $n$ and $p$ stopped being in equilibrium since they are coupled through $nu$. This interaction has
$
  Gamma tilde n_nu sigma_w = t^(-5/2)
$
When $Gamma tilde H$ they decouple and the $n$ and $p$ ratio is frozen. We can find
$
  n_n/n_p = exp(- Q_n/(k_B T_"freeze")) tilde.eq 0.2
$
with $t_"freeze" tilde 1"s" << tau_n$ and $k_B T_"freeze" tilde.eq 0.8 "MeV"$. This is why nucleosynthesis is very inefficient!

After $t_"freeze"$ we have
$
  p + n harpoons.rtlb "D" + gamma
$
with $B_D$ being carried away by a $gamma_gamma$. We use the nucleosynthetic Saha equation
$
  n_D/(n_p n_n) = 6 ((m_n k_B T)/(pi hbar^2))^(-3\/2) exp[B_D/(k_B T)]
$
We define nucleosynthesis by when $n_D\/n_n = 1$. This gives $t_"nuc" tilde 200"s"$.

After $"D"$ is formed they fuse to form $isotope("He", a: 3)$ and $isotope("H", a: 3)$ by
$
  D + p & harpoons.rtlb isotope("He", a: 3)+gamma \
    D+n & harpoons.rtlb isotope("H", a: 3)+gamma \
  D + D & harpoons.rtlb isotope("He", a: 4) + gamma \
  D + D & harpoons.rtlb isotope("H", a: 3) + p \
  D + D & harpoons.rtlb isotope("He", a: 3) + n
$
Which quickly fuse to form $isotope("He", a: 4)$ by
$
   isotope("H", a: 3) + p & harpoons.rtlb isotope("He", a: 4) + gamma \
  isotope("He", a: 3) + n & harpoons.rtlb isotope("He", a: 4) + gamma \
   isotope("H", a: 3) + D & harpoons.rtlb isotope("He", a: 4) + n \
  isotope("He", a: 3) + D & harpoons.rtlb isotope("He", a: 4) + p \
$
After this nucleosynthesis becomes stuck since $isotope("He", a: 4)$ is very stable and because no stable nuclei have $A=5$. We do have $ isotope("He", a: 4) + D &harpoons.rtlb isotope("Li", a: 6) + gamma \
isotope("He", a: 4) + isotope("H", a: 3) &harpoons.rtlb isotope("Li", a: 7) + gamma \
isotope("He", a: 4) + isotope("He", a: 3) &harpoons.rtlb isotope("Be", a: 7) + gamma $
But no stable nuclei have $A=8$ so nucleosynthesis effectively stops.

== Reionization
We observe
$
  rho_("bary",0) tilde 4.2 times 10^(-28)" kg"" m"^(-3)
$
However, in our solar neighbourhood we find
$
  delta_"sn" = (rho_"sn"-rho_("bary",0))/rho_("bary",0) tilde 2 times 10^7
$
and the Sun has $delta_dot.circle tilde 3 times 10^30$.

Most of the observed $rho_("bary",0)$ is intergalactic gas ($tilde 85%$). Half of this is diffuse with $T < 10^5 "K"$ and $delta <= 0$ while the rest is hot with $10^5 "K" < T < 10^7 "K"$ and $3 < delta < 300$. The final $tilde 15%$ is in stars, galaxies, clusters and the gas between these. Most of the intergalactic gas has $X tilde 1$. This is unexpected since $X tilde 0.1$ around $z tilde 1200$. This implies something has reionized the hydrogen in the Universe.

The reionized gas will scatter $gamma_"CMB"$ at a rate
$
  Gamma = n_e sigma_e c
$
Letting reionization start at $t_*$ we have
$
  tau_* = integral_(t_*)^t_0 Gamma (t) dd(t) = c sigma_e integral_(t_*)^(t_0) n_e (t) dd(t)
$
When $tau_* >> 1$ all $gamma_"CMB"$ would have scattered many times giving no information, the spectrum would be smeared. Since we observe detail in the CMB we know $tau_* < 1$. The smearing we do observe is consistent with $tau_* = 0.066 plus.minus 0.016$.

We assume $n_H + n_p = n_"bary" = n_("bary",0) a^(-3)$. Taking everything to ionize instantly at $t=t_*$ then
$
  n_e =^(t > t_*) n_p = n_("bary",0) a^(-3)
$
while $n_e=0$ for $t < t_*$. Then
$
  tau_* &= underbracket(c sigma_e n_("bary",0), Gamma_0) integral_(t_*)^(t_0) dd(t)/(a(t)^3) \ &= Gamma_0 integral_(a(t_*))^1 dd(a)/(dot(a)a^3) \ &= Gamma_0 integral_(a(t_*))^1 dd(a)/(H(a) a^4) \ &= Gamma_0 integral_0^z_* ((1+z)^2 dd(z))/(H(z))
$
With a two-component universe with matter and dark energy we have
$
  H(z) = H_0 (Omega_(m,0) (1+z)^3 + Omega_(Lambda,0))^(1\/2)
$
Then we find
$
  tau_* = 2/(3 Omega_(m,0)) Gamma_0/H_0 ((Omega_(m,0) (1+z_*)^3 + Omega_(Lambda,0))^(1\/2)-1)
$
We can then compute $z_* = 7.8 plus.minus 1.3$ or $t_* tilde 650 "Myr"$.

The reionization is mostly caused by high energy $gamma$ generated by galaxies and AGN. These would both be present at $t_*$.
