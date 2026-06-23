#import "../../temp.typ": *
#show: chpt-note.with()

= Gravitational Waves
The existence of GWs was proposed by Poincaré in 1905. Using relativity Einstein and Rosen concluded in 1936 that GWs did not exist. Robertson#footnote[Of Friedmann-Robertson-Walker fame.] refereed their paper and found a mistake changing the conclusion. The first direct detection of GWs was published in 2015. However, searches has been done since the 1960s and the first indirect observion was done through the _Hulse-Taylor system_. This system is a NS binary with one NS being a pulsar. The measured change in orbital period $dot(P)$ was consistent with the emission of GWs.

== Basics of GWs
According to relativity space is a four-dimensional manifold with gravity being caused by the curvature of the manifold. The object of interest in relativity is the metric $g_(mu nu)$ defined by
$
  dd(s^2) = g_(mu nu) (x) dd(x^mu, x^nu)
$
Minkowski spacetime has $g_(mu nu) = eta_(mu nu)$ with $eta_(mu nu)$ being the Minkowski metric. When spacetime is not Minkowski the metric can be much more complicated. An example of such a metric is the _Friedmann-Robertson-Walker metric_ describing a homogeneous and isotropic Universe
$
  dd(s^2) = dd(t^2) - (a^2(t))/c^2 dd(Sigma^2)
$
with
$
  dd(Sigma^2) = dd(r^2)/(1-k r^2) + r^2 dd(Omega^2)
$
The metric is related to the energy-momentum tensor $T_(mu nu)$ by the _Einstein field equations_
$
  cal(R)_(mu nu) - 1/2 g_(mu nu) cal(R) = (8 pi G)/c^4 T_(mu nu)
$
We consider a non-static source. Away from the source spacetime should be $tilde$ Minkowski with small distortions. Then
$
  g_(mu nu) (x) = eta_(mu nu) + h_(mu nu) (x)
$
with $h_(mu nu)$ small. We plug this metric into the EFE with $T_(mu nu) = 0$. Then ignoring higher order terms we have
$
  square h_(mu nu) = 0
$
which is a wave-equation! We make the ansatz
$
  h_(mu nu) = h_(mu nu)^0 e^(i k_mu x^mu) = h_(mu nu)^0 e^(i bold(k) dot bold(x)) e^(-i omega t)
$
where $h_(mu nu)^0$ is a symmetric $4 times 4$ matrix and $omega = k c$. Let the propagation be along $hat(z)$. Then in the _transverse-traceless gauge_ we have
$
  h_(mu nu)^0 = mat(0, 0, 0, 0; 0, h_+, h_times, 0; 0, h_times, -h_+, 0; 0, 0, 0, 0)
$
with $h_+$ and $h_times$ being constant amplitudes. With $T_(mu nu) eq.not 0$ we find the inhomogeneous wave-equation
$
  square h_(mu nu) = (16 pi G)/c^4 T_(mu nu)
$
which can be solved using Green's functions
$
  h_(mu nu) (x) = (4 G)/c^4 integral dd(x', 3)/abs(bold(x)-bold(x)') T_(mu nu) (t',x')
$
where we evaluate $T_(mu nu)$ at the _retarded time_
$
  t' = t - abs(bold(x)-bold(x)')/c
$

We will consider sources of size $R$ that vary harmonically over time with frequency $nu_s$. Also, we assume $T_(mu nu)$ is dominated by the rest mass of the rotating objects whose density is $rho_m$. We also need some approximations
$
  "long wavelength approx.  " & lambda >> R \
        "far field approx.  " & r >> R
$
With these we have
$
  h_(mu nu) (t,bold(x)) tilde.eq (4 G)/(r c^4) integral dd(x', 3) T_(mu nu) (t- r/c, bold(x)')
$
Then
$
  h_(i j) tilde.eq (4 G)/(r c^4) dv(Q_(i j), t, 2)
$
with $Q_(i j)$ being the _mass quadropole moment_
$
  Q_(i j) = integral dd(x, 3) (x_i x_j - 1/3 r^2 delta_(i j)) rho_m (bold(x))
$
== Energy of GWs
Any GW will distort the distance $L$ between masses by $Delta L$. We quantify the distortion using the _strain_ $h$ defined by
$
  h = (Delta L)/L
$
which will oscillate with the frequency of the GW.

The GW radiation is quadropolar to lowest order. There is no monopole- or dipole radiation due to energy conservation and momentum conservation respectively.#footnote[Since dipole radiation is associated with movement of the CM.] Also, GW radiation is not spherically symmetric just like EM radiation. The GW flux in a volume $V tilde L^3$ with $L >> lambda$ is given by
$
  F = 1/(32 pi) dot(abs(h))^2 c^3/G
$
and the total luminosity is given by
$
  L = G/(5 c^5) sum_(i,j=1)^3 (dv(Q_(i j), t, 3))^2
$
Note, as opposed to EM radiation GWs with $nu lt.tilde "kHz"$ come from sources with sizes $R lt.tilde lambda$ implying the signal is _coherent_. Also, the scaling is more favorable as
$
  "EM flux" tilde 1/r^2
$
while
$
  "strain" tilde 1/r
$
Then GWs are detectable at $nu lt.tilde "kHz"$, however, _graviton_ energies $hbar omega$ are small and their detection therefore seems unfeasible. GWs also propagate almost freely allowing us to in principle detect GWs from the very early Universe. Thereby we can probe beyond the CMB.

== Binary mergers
We consider two masses $m_1$ and $m_2$ in a binary orbit with orbital frequency $omega_s$.#footnote[We will take their velocities to be non-relativistic.] The angular velocity is related to the orbit size $R$ by
$
  omega_s^2 = (G M)/R^3
$
with
$
  M = m_1 +m_2
$
The geometry we consider is shown in @binarysystem.

#let dbinary = diagram(
  node-fill: black,
  label-sep: 0.25em,
  node((-4, 0), name: <0>),
  edge(<0>, <1>, "->", label: $x$, label-pos: 100%),
  node((4, 0), name: <1>),
  node((0, 2), name: <2>),
  edge(<2>, <3>, "->", label: $y$, label-pos: 100%),
  node((0, -2), name: <3>),
  node((0, 0), name: <CEN>),
  node((3, -1), name: <M2>, radius: .2em),
  node((-3, 1), name: <M1>, radius: .2em, label: $m_2$),
  edge(<CEN>, <M2>, "--", label: $m_2$, label-pos: 100%),
  edge(<CEN>, <M1>, "--", label: $m_1$, label-pos: 100%),
  edge(
    (1, 0),
    (1, -0.33),
    bend: -45deg,
    label: $omega_s t$,
    label-pos: 55%,
  ),
)

#figure(
  scale(dbinary, 100%),
  caption: [Geometry for the binary system.],
)<binarysystem>

Also, denote the distance from $(0,0)$ to $m_i$ by $r_i$. We have
$
  rho_m (bold(x)) = sum_(alpha=1)^2 m_alpha delta(bold(x)-bold(x)^alpha)
$
implying
$
  Q_(i j) &= integral dd(x, 3) (x_i x_j - 1/3 r^2 delta_(i j)) rho_m (bold(x)) \
  &= sum_(alpha=1)^2 m_alpha (x_i^alpha x_j^alpha - 1/3 r_alpha^2 delta_(i j))
$
We have $z=0$ for $m_alpha$ so
$
  Q_33 = sum_(alpha=1)^2 m_alpha (-1/3 r_alpha^2)
$
and
$
  Q_(3i) = Q_(i 3) = 0
$
Also,
$
  Q_11 & = sum_(alpha) m_alpha ((x_1^alpha)^2 - 1/3 r_alpha^2) \
       & = sum_alpha m_alpha (2/3 (x_1^alpha)^2 - 1/3 (y_2^alpha)^2)
$
and
$
  Q_22 & = sum_alpha m_alpha (2/3 (y_1^alpha)^2 -1/3 (x_2^alpha)^2) \
  Q_12 & = Q_21 = sum_alpha m_alpha x^alpha y^alpha
$
With $R = r_1+r_2$ we can write
$
  Q_(i j) = 1/2 mu R^2 J_(i j)
$
with
$
  mu =(m_1 m_2)/(m_1+m_2)
$
and
$
  J_(i j) = mat(cos 2 omega_s t + 1/3, sin 2 omega_s t, 0; sin 2 omega_s t, 1/3 - cos 2 omega_s t, 0; 0, 0, -2/3)
$
We can approximate
$
  J_(i j) tilde cos 2 omega_s t";  " Q_(i j) tilde Q";  " h_0 tilde h_+ tilde h_times
$
so
$
  dv(Q, t, 2) = 1/2 mu R^2 (4 omega_s^2) cos 2 omega_s t
$
implying
$
  h(t) & tilde (4 G)/(r c^4) 2 mu R^2 omega_s^2 cos omega_("GW") t \
       & = h_0 cos omega_"GW" t
$
with $omega_"GW" = 2 omega_s$. We can rewrite $h_0$ as
$
  h_0 & tilde (4 G)/(r c^4) 2 mu R^2 omega_s^2 \
      & = (4 G)/(r c^4) 2 mu R^2 (G M)/R^3 \
      & = 2 ((2 G M)/(c^2 r))((2 G mu)/(c^2 R)) \
      & = 2/(r R) ((2 G m_1)/c^2) ((2 G m_2)/c^2) \
      & = (2 R_1 R_2)/(r R)
$
where $R_i$ are the masses Schwarzschild radii! As an example consider a binary NS system with $m_1 tilde m_2 tilde 1.4 M_dot.o$ so $R_i tilde 4"km"$. With $r = 40"Mpc"$ and $R=100"km"$ we find
$
  h_0 tilde 10^(-22) tilde (Delta L)/L
$
which is tiny.#footnote[This is why constructing GW detectors is hard. One needs to take extreme care to remove noise.]

We can find the luminosity by computing $dot.triple(Q)_(i j)$ and we find
$
  L = (32 G)/(5 c^5) (mu R^2 omega_s^3)^2
$
and the flux is
$
  F tilde 10^(-5) "W"/"m"^2
$
for $omega_s tilde 400 "Hz"$. This is a huge amount of flux!

== The inspiral stage
The energy of the binary system is given by
$
  E & = 1/2 mu omega_s^2 R^2 - (G m_1 m_2)/R \
    & = -1/2 m_1 m_2 G/R \
    & = - (G mu M)/(2 R)
$
During the insprial stage GW are emitted leading to orbital decay. The energy loss is
$
  -dv(E, t) & = -1/2 (G mu M)/R dot(R)/R \
            & = 1/3 G mu M dot(omega)_s/omega_s \
            & =^! L = (32 G)/(5 c^5) mu^2 R^4 omega_s^6
$
We eventually find
$
  dot(omega)_s^3 = (96/5)^3 omega_s^11/c^15 (G m_c)^5
$
where we define the _chirp mass_ $m_c$ by
$
  m_c = (mu^3 M^2)^(1\/5) = (m_1 m_2)^(3\/5)/(m_1+m_2)^(1\/5)
$
which is the characteristic scale of the inspiral stage. We can use $omega_s = pi nu_"GW"$ to obtain
$
  m_c = c^3/G [5/96 pi^(-8\/3) nu_"GW"^(-11\/3) dot(nu)_"GW"]^(3\/5)
$
Then we can compute $m_c$ by knowing $nu_"GW"$ and $dot(nu)_"GW"$! This is done in practice by using
$
  nu_"GW" = 1/(Delta t)";  " dot(nu)_"GW" = (Delta nu_"GW")/(Delta t)
$

== The merger stage
When _coalescence_ or merging begins our Newtonian approximation breaks down horribly. We expect this happens when the orbital radius corresponds to the Schwarzschild radius of the combined masses
$
  R_s = (2 G M)/c^2
$
implying
$
  omega_"max" = 1/sqrt(8) c^3/(G M) tilde 2 times 10^5 "Hz" (M_dot.o/M)
$
This frequency should correspond to the maximally observed GW frequency
$
  omega_"max" tilde.eq pi nu_"GW"^"max"
$
We can then measure $nu_"GW"^"max"$ to compute the total mass $M$. With the chirp mass $m_c$ we can now determine $m_1$ and $m_2$. We define the _mass ratio_ $q$ by
$
  m_1 = q M";  " m_2 = (1-q) M
$
implying
$
  M = m_c/(q(1-q))^(3\/5)
$
Then $q$ can be inferred!

== Luminosity distance
The _luminosity distance_ $d_L$ is defined by
$
  F = L/(4 pi d_L^2)
$
Using the above equations for $F$ and $L$ we find
$
  d_L = (8 G)/(sqrt(5) c^4) 1/h_0 (mu R^2 omega_s^4)
$
We can use this to probe the _Hubble constant_ $H_0$ since
$
  d_L (z_0) = (1+z_0) c integral_0^(z_0) dd(z)/(H(z))
$
Then knowing $z_0$ allows us to probe these cosmological parameters using GWs. This is likely not viable for BH mergers since there is no EM counterpart to measure $z_0$. However, we can apply this to NS mergers which have been observed!

Also, GWs are redshifted over cosmological distances as
$
  nu_"GW"^"obs" = nu_"GW"^"emit"/(1+z)
$
and
$
  dot(nu)_"GW"^"obs" &= dv(nu^"obs", t^"obs") = underbracket(dv(nu^"obs", nu^"emit"), (1+z)^(-1)) dv(nu^"emit", t^"emit") underbracket(dv(t^"emit", t^"obs"), (1+z)^(-1)) \
  &= dot(nu)_"GW"^"emit"/(1+z)^2
$
implying the chirp mass is also redshifted as
$
  m_c^"obs" = (1+z) m_c^"emit"
$
and similarly for $m_1^"obs"$ and $m_2^"obs"$.

== Spinning black holes
The spin of BHs leads to velocity-dependent interaction during inspiral which are annoying to deal with analytically. Therefore such effects are typically handled numerically and by fitting data.

We consider a massive object of mass $m$ and spin $bold(S)$. Then we define a dimensionless _spin parameter_
$
  chi = (c abs(bold(S)))/(G m^2)
$
Spinning or Kerr BHs have smaller horizons as compared to Schwarzschild BHs which have $chi = 0$. When $chi = 1$ one can show
$
  R_"horizon" = (G m)/c^2 = R_"Schwarzschild"/2
$
The $chi$ of the final BH can be determined analytically by modelling the _ringdown stage_ as a damped oscillator. As an example in GW150914 the more massive BH had $chi < 0.7$ and it was found that the final BH had
$
  chi tilde.eq 0.67
$
in agreement with the spin of the _primary_ BH before the merger.

== Pulsar timing arrays
The interferometers we use to measure GWs are not sensitive to GWs with frequencies below $tilde 10"Hz"$. These frequencies would occur for mergers of SMBHs with $M gt.tilde M_dot.o$. However, such mergers could cause a _stochastic GW background_ which we can measure using _pulsar timing arrays_.

Using PTAs is based on the observation that millisecond pulsars have extremely stable rotation periods
$
  nu (t) = nu_0 + dot(nu)_0 t
$
with
$
  dot(nu)_0/nu_0 tilde 10^(-20)"Hz"
$
As EM radiation propagates through spacetime where the metric is perturbed due to GWs it accumulates a time delay $Delta T$. We consider a pulsar a distance $L$ away from us. Then
$
  Delta T(t) = 1/2 hat(p)^i hat(p)^j integral_0^L dd(s) h_(i j) (t(s),bold(x)(s))
$
where $hat(p)$ is "along" $L$ and
$
  t(s) = t-(L-s)";  " bold(x) = (L-s) hat(p)
$
with $c = 1$. We consider a binary merger a distance $d$ away. We define $hat(Omega)$ to be along $d$ toward us and denote the wavelength of the produced GWs by $lambda$. We can rewrite the above in terms of redshift
$
  z(t) &= dv(Delta T(t), t) \ &= 1/2 hat(p)^i hat(p)^j integral_0^L dd(s) pdv(h_(i j), t) (t(s),bold(x)(s))
$
We consider pairs of pulsars call each member of a pair $a$ and $b$. The time averaged product over an observation time $T$ is then
$
  rho_(a b) & = overline(z_a (t) z_b (t)) \
            & = 1/T integral_(-T\/2)^(T\/2) dd(t) z_a (t) z_b (t)
$
with
$
  z_a (t) = 1/2 hat(p)_a^i hat(p)_a^j integral_0^L_a dd(s) pdv(h_(i j), t) (t-(L_a-s), (L_a-s) hat(p)_a)
$
and similarly for $z_b (t)$. We average over all pairs of pulsars $a$ and $b$ separated by an angle $gamma$ on the sky
$
  Gamma(gamma) & = expval(rho_(a b))_(a, b in gamma) \
               & = [h^2 mu_u (gamma) (1+delta_(a b))]_(a,b in gamma)
$
where $mu_u$ is the _Hellings and Downs function_. Then comparing $Gamma(gamma)$ with observations could show the existence of a GW background. The explanations for the GW background is usually mergers of SMBH or sources from the very early Universe such as inflation!

