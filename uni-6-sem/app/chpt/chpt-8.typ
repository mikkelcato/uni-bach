#import "chpt-temp.typ": *
#show: chpt-note.with()

= Evidence for Dark Matter
== The Coma cluster
Fritz Zwicky showed in 1933 by applying the virial theorem that visible matter was not enough to account for the observed dynamics of the Coma cluster.

We retrace Zwicky's calculation. By the virial theorem
$
  2 expval(E_"kin") = - expval(E_"pot")
$
We assume the Coma cluster is spherical with
$
  expval(E_"pot") = - alpha (G M^2)/r
$
where $alpha tilde$ density profile and $M$ is the total mass. Computing $expval(E_"kin")$ requires knowing the average velocity relative to the COM in the Coma cluster. We use the dispersion
$
  sigma^2 = expval(v^2) - expval(v)^2
$
However, we can only measure $sigma_parallel$ so we need to assume the Coma cluster is isotropic implying
$
  sigma = sqrt(3) sigma_parallel
$
Also, assuming the Coma cluster is uniform we have
$
  alpha = 3/5
$
Then by the virial theorem
$
  sigma_parallel^2 = (G M)/(5 R)
$
with $R$ being the radius of the Coma cluster. Zwicky measured seven galaxies in the Coma cluster and found
$
  sigma_parallel tilde 1000 "km"/"s"
$
Likewise, he measured
$
  R tilde 1 "Mpc"
$
Then
$
  M_"virial" tilde 10^14 M_dot.o
$
Also, taking the luminosity of the Coma cluster to be
$
  L_"Coma" tilde [10^12,10^13] L_dot.o
$
Assuming the Coma cluster has $N_"Coma" = 800$ galaxies with $M_"Galaxy" tilde 10^9 M_dot.o$ we find
$
  M_"count" tilde 10^12 M_dot.o
$
This is much less than $M_"virial"$ and we can define the _mass-to-light-ratio_
$
  M_"virial"/L_"Coma" tilde 5 dot 10^2
$
We can make much more conservative estimates and still miss $tilde 80%$ of the mass.

This was the $tilde$ first evidence  for dark matter.

== Rotation curves
Vera Rubin studied the velocities of stars in galaxies in the 70s and 80s deriving their rotaiton curves $v(r)$
$
  v(r) = sqrt((G M(r))/r)
$
We observe almost all the luminous matter in a galaxy to be toward the center bulge. Therefore at large $r$ we would expect $M(r) tilde "constant"$ and $v tilde r^(-1\/2)$. However, Rubin found $v(r) tilde "constant"$ for large $r$ implying $M(r) tilde r$ and $rho tilde r^(-2)$. Then galaxies appear to be surrounded by a halo of non-luminous dark matter. This is shown to very large $r$ using the $21"cm"$ line.

One of the most commonly used density profiles for the dark matter halo is the Navarro-Frenk-White profile
$
  rho_"NFW" (r) = rho_s r_s/r (1+r_s/r)^(-2)
$
Just to give an example.

== N-body simulations
Modern _proof_ of dark matter. Using simulations people have shown we need dark matter to achieve proper structure formation due to the clumping nature of dark matter.

== Gravitational lensing
We know from relativity that matter can deform spacetime and therefore bend light. Then a large enough mass will act as a lens. This could be a quasar, etc. or more interestingly dark matter. We can use this to learn things about the source or attempt to reconstruct the lens.

An example where this has been used is the Bullet cluster where dark matter is needed to explain the amount of gravitational lensing. The Bullet cluster consists of two clusters colliding. When this happens the dark matter halos of the cluster pass through eachother uninhibited while clouds of gas are slowed due to friction. This acts to separate the luminous matter and dark matter from which one can show the lensing happens mostly where there is no visible matter.

== MOND
Modified Newtonian dynamics or MOND attempts to explain what dark matter does without dark matter. MOND assumed Newtonian dynamics fail at small accelerations. This is modelled by
$
  mu(a/a_0) a = -abs(grad phi.alt_N)
$
where $mu$ is a function going to $1$ as $a>>a_0$ with $a_0 tilde 10^(-10) "m" "s"^(-2)$ and
$
  phi.alt_N = - (G M)/r
$
With rotation curves we have $a << a_0$ at large $r$ so
$
  mu(a/a_0) tilde a/a_0
$
and
$
  a^2/a_0 = (G M)/r^2 => a^2 = a_0 (G M)/r^2
$
Using $a = v^2\/r$ we find
$
  v = (a_0 G M)^(1\/4) tilde "constant"
$
so MOND leads to $v tilde "constant"$ at large $r$, however, there is little else it explains and it has no real motivation. Briefly, it is stupid.

= $Lambda$CDM
== Hubble's law
Hubble's law is the observation that galaxies move away from us with a velocity given by
$
  v = H_0 d
$
where
$
  H_0 tilde 68 "km" "s"^(-1) "Mpc"^(-1)
$
is the Hubble constant.

Within relativity this is a consequence of the Universe and spacetime expanding isotropically and homogeneously. Under these assumptions any observer would see Hubble like behaviour. We write
$
  d = a(t) chi
$
where $chi$ is the comoving distance between objects measured at $t = t_0$ and $a(t)$ is the scale factor normalised as $a(t_0) = 1$. Then
$
  v equiv dot(d) = dot(a) chi = dot(a)/a d = H(t) d
$
with
$
  H(t) equiv dot(a)/a
$
being the Hubble parameter where $H_0 = H(t_0)$. We see $H(t)$ is a measure of the expansion rate of the Universe.

== The FRW cosmology
We would like equations describing the expansion of the Universe more qualitatively. We make three assumptions:

1. We can use relativity.

2. The cosmological principle holds.#footnote[This is just isotropy and homogeneity at large scales.]

3. There exists a cosmic time.

Under these we can derive the Friedmann-Robertson-Walker metric
$
  dd(s^2) = dd(t^2) - a^2 (t) [dd(r^2)/(1-k r^2) + r^2 dd(Omega^2)]
$
where $k$ is related to spatial curvature.

With this metric and assuming the Universe is a perfect fluid one finds the two Friedmann equations
$
      H^2 + k/a^2 & = (8 pi G)/3 rho + Lambda/3 \
  dot.double(a)/a & = - (4 pi G)/3 (rho + 3 p) + Lambda/3
$
which can be combined to give
$
  dv(, t) (rho a^3) = -p dv(, t) (a^3)
$
We need an equation of state $P(rho)$ to solve the Friedmann equations. We will typically assume
$
  P = w rho
$
implying
$
  rho tilde a^(-3(w+1))
$
We see for non-relativistic matter where $w = 0$ we find
$
  rho_m tilde 1/a^3
$
while for relativistic matter where $w = 1/3$ we find
$
  rho_r tilde 1/a^4
$
and for $rho_Lambda tilde "constant"$ we need $w = -1$. Also, we can show
$
  a(t) tilde t^(2\/3(w+1))
$
Consider a de-Sitter Universe i.e. one dominated by $Lambda >0$. Then
$
  a(t) tilde e^(H t)
$
with
$
  H = sqrt(Lambda/3)
$
Because different types of matter evolve differently the Universe has experienced various periods of domination and equality. This is obvious by considering the above scalings.

== Typical form
Assume $Lambda = 0$ then
$
  k/a^2 = (8 pi G)/3 rho - H^2
$
we define
$
  rho_"crit" = (3 H^2)/(8 pi G)
$
which corresponds to $rho_"Universe"$ if $k =^! 0$. We define density parameters
$
  Omega_i = rho_i/rho_"crit"
$
and
$
  Omega_k = -k/(H^2 a^2)
$
Then the $1^"st"$ Friedmann equation becomes
$
  Omega_m + Omega_r + Omega_Lambda + Omega_k = 1
$
We could also write
$
  H^2/H_0^2 = Omega_(0, Lambda) + (Omega_(0,k))/a^2 + (Omega_(0,m))/a^3 + (Omega_(0,r))/a^4
$
Note, we sometimes use $1+z = a^(-1)$.

We would now like to measure the density parameters!

== The CMB
We know the CMB is a _perfect_ blackbody with
$
  T_("CMB", 0) tilde.eq 2.73
$
Within the Big Bang model we believe the CMB to originate when the Universe cooled sufficiently so $p$ and $e^-$ could form hydrogen
$
  p + e^- -> H + gamma tilde "recombination"
$
When recombination occured the temperature was $T_("CMB","rec") = 3700 "K"$ implying
$
  z_"rec" tilde 1300
$
During this epoch $gamma$ could not propagate freely due to Thomson scattering with free $e^-$. We have
$
  Gamma_T tilde n_e sigma_T c tilde X(z) n_B sigma_T c
$
where $n_e = X(z) n_B$. When recombination begins $X(z)$ will decrease.

After a long enough time $gamma$ will eventually propagate freely. We call this _decoupling_ and is defined by
$
  X(z_"dec") n_B sigma_T c =^! H(z_"dec")
$
implying
$
  z_"dec" tilde 1100
$
Note, the $gamma$ will of course still scatter with stuff in the Universe.

The CMB is important due observed temperature fluctuations $dd(T, d: delta)$ defined by
$
  dd(T(theta,phi), d: delta)/expval(T) = (T(theta,phi) - expval(T))/expval(T)
$
The largest anisotropy is a dipole anisotropy on the order of $10^(-3)$ caused by our peculiar motion. We can characterise the dipole by
$
  T^* = T (1 + v/c cos theta)
$
with $v tilde.eq 370 "km""s"^(-1)$. We usually remove the dipole after which we are left with fluctuations on the order of $10^(-5)$.

We would like to analyse the fluctuations. This is done by expanding in terms of spherical harmonics $Y_(l m)$
$
  dd(T(theta,phi), d: delta)/expval(T) = sum_(l=0)^oo sum_(m=-l)^l a_(l m) Y_(l m) (theta,phi)
$
with
$
  a_(l m) = integral_0^(2 pi) dd(phi) integral_0^pi dd(theta) sin theta dd(T, d: delta)/expval(T) Y_(l m)
$
Consider the correlation function $C(alpha)$
$
  C(alpha) &= expval(dd(T, d: delta)/expval(T) (hat(u)) dd(T, d: delta)/expval(T) (hat(u)'))_(hat(u) dot hat(u)' = cos alpha) \
  &= 1/(4 pi) sum_(l=0)^oo (2 l + 1) C_l P_l (cos alpha)
$
where $P_l$ are the Legendre polynomials and $C_l$ are the multipole moments given by
$
  C_l = 1/(2 l+1) sum_(m=-l)^l abs(a_(l m))^2
$
The multipole $l$ is related to an angular size $alpha$ by
$
  alpha = (180 degree)/l
$
We work with the _power spectrum_ normalised as
$
  D_l = (l(l+1))/(2 pi) C_l
$
and
$
  expval(dd(T, d: delta)^2) = D_l expval(T)^2
$
== $Omega_i$ from the CMB
Quantum fluctuations pre-inflation lead to density fluctuations in the early Universe after they are stretched to cosmic scales. These generate potential hills and wells acting on the photon-baryon plasma before recombination. Gravity leads to a pressure in this potential landscape which works to collapse the plasma. We likewise have a radiation pressure which will resist this collapse. The plasma then oscillates in a damped fashion. On potential hills we have rarefaction of density which acts to cool down the plasma, while compressions heats the plasma up.

These oscillations are imprinted on the power spectrum as peaks. The fundamental mode has
$
  k_1 = pi/lambda_"sound"
$
with $lambda_"sound"$ being the sound horizon. We can imagine the peaks as corresponding to standing waves in the plasma where the temperature anisotropy would be maximal.#footnote[The first peak corresponds to one compression and $l$ gives the scale of $lambda_"sound"$. The second peak corresponds to one compression and one rarefication, etc.]

We can immediately determine $Omega_k$ from the power spectrum as the peaks would shift due to $gamma$ propagation in a curved Universe being different. The CMB we observe is consistent with $Omega_k tilde.eq 0$. This is usually done by using the first peak.

The baryon density $Omega_b$ can be determined by considering the relative size of peaks. Making $Omega_b$ larger would lead to more matter being compressed during infall into potential wells. However, the rebound would be unchanged. This implies compressions or odd peaks are enhanced. When we add mass the oscillations also slow down leading to a shift toward higher $l$. The matter density $Omega_m$ also affects the peaks as a higher $Omega_m$ would weaken the driving force so peak amplitudes would decrease. The third peak is most sensistive to $Omega_m$ as the first two peaks depend heavily on $Omega_b$.

The damping we observe at large $l$ generally acts like a self-consistency check as it depends on everything.

Also, at small $l$ we have large uncertainties referred to as _cosmic variance_ since we only have a single sky we can measure.

Using modern data we can infer
$
  Omega_Lambda & tilde.eq 0.68 \
       Omega_c & tilde.eq 0.27 \
       Omega_b & tilde.eq 0.05
$
However, there are certain problems with the $Lambda"CDM"$ model presented above such as the _Hubble tension_.

== BBN
Alternative way to determine $Omega_b$ which we do not cover. See any book on cosmology.

#pagebreak()
= Dark Matter Candidates
== Requirements for DM
We evidently need at large amount of dark matter, so what is it?

We know:

1. DM is unlikely to consist of baryons. We are quite sure that $Omega_b tilde 0.05$. This would also be consistent with structure formation.

2. Otherwise:

  1. DM should only interact _weakly_ with EM radiation. This is what makes it dark.

  2. DM must have the correct density $Omega_c tilde 0.25$.

  3. DM must be very light or very heavy, else we would have found them already.

  4. DM must be stable against decay on cosmological scales.

== dMACHOs
There is still the possibility that DM is baryonic and concentrated in _dark massive compact halo objects_ or dMACHOs. The belief originally was that DM consisted of BHs, NS, etc. However, we would likely have observed these objects already through gravitational lensing.

== WIMPs
Assuming DM particles are heavy they should interact weakly. We refer to such particles as _weakly interacting massive particles_ or WIMPs. We will denote them by $chi$. One potential candidate is the _neutralino_ which is a proposed particle in supersymmetric extensions of the Standard Model. WIMPs should be neutral and if we assume there only is one type of WIMP then it should be its own antiparticle.#footnote[i.e. a _Majorana particle_.]

WIMPs are believed to be produced in the early Universe through particle collisions, in particular through particle-antiparticle collisions. At temperatures $k_B T >> m_chi c^2$ particle-antiparticle pairs had enough energy to produce WIMPs at a rate#footnote[We assume DM is cold. This is consistent with observed LSS.]
$
  Gamma_chi = expval(sigma v) n_chi
$
These production processes would have been in equilibrium with WIMP annihilation. When the Universe expanded the temperature decreased and the number of particles able to produce WIMPs would drop exponentially leading to an equilibrium density
$
  n_chi tilde ((m_chi T)/(2 pi))^(3\/2) exp[-(m_chi c^2)/(k_B T)]
$
Also, $n_chi$ would decrease due to the expansion of the Universe. When
$
  Gamma_chi tilde H
$
WIMP production decouples. After this point $n_chi$ decreases as $a^(-3)$. We can show the decoupling temperature $T_chi$ is given by
$
  k_B T_chi tilde (m_chi c^2)/x
$
with $20 lt.tilde x lt.tilde 50$ for $10 "GeV" lt.tilde m_chi c^2 lt.tilde 10 "TeV"$ and
$
  Omega_chi/0.2 tilde ((m_chi c^2\/k_B T)/20) ((3 times 10^(-26) "cm"^3 "s"^(-1))/expval(sigma v))
$
Then the _relic density_ decreases as $expval(sigma v)$ increases. This is the case since the WIMPs are able to stay in thermodynamic equilibrium for longer, which suppreses their density by the Boltzmann factor in $n_chi$.

Assuming $x tilde 20$ and $v tilde 1/3 c$ we require $sigma tilde 3"pb"$ for the proper $Omega_chi$. Taking the interaction to be on the order of the weak interaction $g_chi tilde g_W$ then
$
  sigma tilde g_W^4/m_chi^2
$
implying
$
  m_chi tilde 100 "GeV"
$
for $sigma tilde 3"pb"$. This is called the _WIMP miracle_.

== Detection of WIMPs
The detection candidates for many DM candidates are generally model dependent. We necessarily also carry a prejudice with us when considering candidates since we can only probe certain things. We always need some _signature_ unique to a given model!

One can imagine many different ways to try and observe DM using astroparticle observables.#footnote[e.g. $gamma$- and cosmic rays.] We can generally consider three different search strategies. We can attempt indirect detection by considering a process of the form
$
  underbracket(chi + chi, "dark sector") ->^"model dependent process" underbracket(f + overline(f), "Standard Model")
$
with $chi$ being a DM candidate and $f$ being some known particle in the Standard Model. We could then possibly observe $f$ and $overline(f)$ and try to check predictions from some model. Likewise, we could attempt to go the other way by a process of the form
$
  f + overline(f) -> chi + chi
$
using colliders. Lastly, we could possible conceive a detector being able to directly detect $chi$.

== Indirect searches for DM
Consider the case of WIMP DM annihilation. We could have annihilation of the form
$
  chi + chi ->^"model dependent" cases(W^-\/Z\/q, W^+\/Z\/overline(q)) ->"decays" cases(pi^0, pi^+, pi^-)->^"further decays" cases(gamma"-rays", nu, p\/overline(p)", " d\/overline(d)", " "anti-matter")
$
or
$
  chi + chi -> underbracket(l + overline(l), "lepton pairs") ->^("Bremsstrahlung", "Synchrotron", "Compton") gamma
$
or
$
  chi + chi -> gamma
$
The expected $gamma$ flux from WIMP annihilation can be written as
$
  dd(Phi_gamma)/dd(E_gamma, Omega) (E_gamma, Psi) = overbracket(1/(4 pi) integral_("l.o.s") rho_chi^2 (bold(r)) dd(cal(l) (Psi)), "amount of DM") underbracket([expval(sigma v)_"ann"/(2 S_chi m_chi^2) sum_f dv(N_gamma^f, E_gamma) B_f], "amount of" gamma)
$
where $Psi$ is the direction of observation. $dd(cal(l)(Psi))$ is the line of sight. $rho_chi (bold(r))$ is the DM density at the point $bold(r)$, which we square because we are talking about annihilation. $expval(sigma v)_"ann"$ is the velocity averaged annihilation cross section. $S_chi$ is the symmetry factor. $dd(N_gamma^f)\/dd(E_gamma)$ is the $gamma$-ray spectrum for each annihilation channel $f$ with branching ratio $B_f$. We can write
$
  dv(J, Omega) = 1/(4 pi) integral_("l.o.s") rho_chi^2 (bold(r)) dd(cal(l)(Psi)) tilde "astrophysics"
$
Then
$
  dv(Phi_gamma, E_gamma) (E_gamma, dd(Omega, d: Delta)) = underbracket(1/(4 pi) integral_(dd(Omega, d: Delta)) dd(Omega) integral_("l.o.s") rho_chi^2 (bold(r)) dd(cal(l)(Psi)), J tilde "astrophysics") overbracket([expval(sigma v)_"ann"/(2 S_chi m_chi^2) sum_f dv(N_gamma^f, E_gamma) B_f], "particle physics")
$
Assuming DM has uniform density $rho_chi$ is annihilating in a spherical galaxy of radius $r$ at distance $d$ we find
$
  J tilde.eq (4 pi r^3 rho_chi^2)/(3 d^2)
$
Then nice search targets have high $rho_chi$ are close (small $d$), large (big $V$) and have well understood backgrounds.

An example people use are dwarf spheroidal galaxies in the DM halo of the Milky Way. Assuming steady state and spherical symmetry one can approximate the DM density through the $2^"nd"$ Jeans equation which then provides an estimate of the $J$ factor.#footnote[using measured velocities of stars.] The dSphs also have high mass to light ratios and little background contamination, thereby making them excellent search targets.

The idea is then to pick our favourite dSph and compute the corresponding $J$ factor. We can model the emissions as
$
  dv(Phi_gamma, E_gamma) (cal(l),b) = underbracket(dv(Phi_(gamma,"diffuse"), E_gamma) + dv(Phi_(gamma,"PS"), E_gamma), "background") + underbracket(dv(Phi_(gamma,"DM"), E_gamma), "DM signal")
$
We can then compare with $gamma$-ray data using the maximum likelihood method and see if a model including DM is favorable. One uses the Poisson likelihood
$
  scr(L) (m_chi, expval(sigma v), theta | D_(i j)) = product_(i,j) (e^(-mu_(i j)) mu_(i j)^(N_(i j)))/(N_(i j) !)
$
where $N_(i j)$ are the number of observed $gamma$-ray counts and $mu_(i j) equiv mu_(i j) (m_chi, expval(sigma v), theta)$ are the expected number of counts. Here, $theta$ represents all nuisance parameters which are not of interest. We can compute $mu_(i j)$ by various methods. Then we can do a likelihood ratio test
$
  "TS" = -2 sum_(i,j) ln [(scr(L)(m_chi, expval(sigma v) =^! 0, hat(theta) | D_(i j)))/(scr(L)(m_chi, hat(expval(sigma v)), hat(theta) | D_(i j)))]
$
where $hat(square)$ represents parameters optimised in a fit to data with fixed $m_chi$. When $"TS"$ is large their is a preference for including DM.

When people have done the above they found no significant detection of a DM signal. However, the data we have can be used to constrain $expval(sigma v)$. This is simply done by another hypothesis test comparing a fixed $expval(sigma v)$ to the global best fit value. These limits obviously become stronger when more dSphs are included.

Aside from using dSphs many other astrophysical sources can obviously be used, however, this is very difficult since one needs to be able to confidently model the background. As an example there has been observed an excess signal from the galactic center. This excess could be due to DM annihilation, but it could also be from various other things. Many of these alternative explanations will probably be ruled out or confirmed eventually as we get better models and detectors.#footnote[e.g. emission from millisecond pulsars is favored.]

== Direct searches for DM
Consider the elastic scattering of WIMPs on nucleons
$
  chi + N -> chi + N
$
We use the lab frame wherein $N$ is stationary. Then
$
  bold(p)_chi = m_chi bold(v)";  " bold(p)_N = 0
$
After scattering we have
$
  bold(q)' = m_chi bold(v)'";  " abs(bold(q)) = sqrt(2 m_N E_(n r))
$
with scattering angles $theta'$ and $theta$.

The velocity of collisionless DM particles should follow the Maxwell-Boltzmann distribtution in a non-rotating halo
$
  f(bold(v)) = N abs(bold(v))^2 exp[-(3 abs(bold(v))^2)/(2 sigma^2)]
$
with
$
      sigma(R) & = sqrt(3/2) v_c (R) \
  v_c(R_dot.o) & = 220" km""s"^(-1)
$
The distribution is truncated at $v_"esc" = sqrt(2 phi.alt) tilde.eq 544 "km""s"^(-1)$ with $phi.alt$ being the Milky Way gravitational potential. Note, since the Earth goes around the Sun $v_c$ varies annually.

The WIMP scattering rate is given by
$
  dv(R, E_(n r)) = rho_chi/m_chi M/m_N integral_(v_"min")^(v_"esc") v f(v) dv(sigma, E_(n r)) dd(v) prop exp[-(4 mu E_(n r))/E_0]
$
With $M$ being the target detector mass. $v_"min"$ is given by
$
  v_"min" = sqrt((E_(n r) m_N)/2 1/mu^2)
$
$mu$ is the reduced mass and $E_0$ is the most probable interaction energy. The number of observed recoils in a time $T$ is given by
$
  N_"obs" = T integral_(E_"min")^(E_"max") dd(E_(n r)) epsilon(E_(n r)) dv(R, E_(n r))
$
with $epsilon$ being the efficiency of detecting recoil at $E_(n r)$ and
$
  E_"max" & = (2 mu v_"esc")/m_N \
  E_"min" & tilde "set by detector"
$
The WIMP scattering cross section is given by
$
  dv(sigma, E_(n r)) = m_N/(2 v^2 mu^2) [sigma_"SI" underbracket(F_"SI"^2 (E_(n r)), "form factor") + sigma_"SD" F_"SD"^2 (E_(n r))]
$
with $square_"SI"$ being the spin independent part arising from scalar or vector interactions and $square_"SD"$ being the spin dependent part arising from axial vector interactions
$
  scr(L)_S &tilde overline(chi) chi overline(q) q \
  scr(L)_V &tilde overline(chi) gamma^mu chi overline(q) gamma_mu q \
  scr(L)_A &tilde overline(chi) gamma^mu gamma_5 chi overline(q) gamma_mu gamma_5 q
$
We can write
$
  sigma_"SI" = (sigma_n mu^2 A^2)/mu_n^2
$
with $sigma_n$ being the WIMP nucleon cross section, $mu_n$ is the reduced mass of the WIMP nucleon system, and $A$ is the atomic number. The $sigma_"SD"$ is more complicated and depends on many things.

We can write the energy loss of WIMPs in detectors as
$
  (dv(E, x))_"tot" = underbracket((dv(E, x))_"elec", "WIMP electron interactions" #linebreak() "excitation or ionisation") + underbracket((dv(E, x))_"nucl", "WIMP nucleon scattering" #linebreak() "heat" -> "phonons")
$
This leads to a breadth of different detectors trying to detect side effects of a recoiling nucleus:

1. Anorganic crystals: $E_"kin"$ of recoiling nucleus

  1. in semiconducting crystals creates electron-hole pairs.

  2. in scintillating crystals excites and ionises atoms in the crystal lattice. This leads to scintillation light.

2. Cryogenic detectors: $E_"kin"$ of recoiling nucleus

  1. leads to heating of crystals and the creation of phonons.

3. Noble liquids in TPCs: $E_"kin"$ of recoiling nucleus

  1. leads to heat.

  2. excited states $X^*$ detected through scintillation.

  3. ionised states detected through scintillation (S1). $bold(E)$-field removes ionised $e^-$ leading to secondary scintillation signal (S2) in gaseous phase.

4. Bubble chambers: $E_"kin"$ of recoiling nucleus

  1. creates bubbles in a superheated liquid.

All of the above are heavily affected by backgrounds which must be reduced. Many detectors have managed to heavily constrain the scattering cross section.

== Complementarity
Comparing direct with indirect searches is usually not possible. We therefore use simplified models to compare with collider searches

#let d1 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a", label: $g_chi$),
    vertex("b", label: $g_q$),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $chi$),
    edge("i2", "a", type: "fermion", label: $chi$),
    edge("a", "b", type: "ghost", label: $M_"med"$),
    edge("b", "f1", type: "antifermion", label: $overline(f)$),
    edge("b", "f2", type: "antifermion", label: $f$),
  ),
)

#figure(
  scale(d1, 100%),
)

We consider EFT to compare direct and indirect searches. Assuming $E_"c.o.m" << M_"med"$ we have
#let d2 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a", shape: "blob", size: .5em, stroke: purple),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $chi$),
    edge("i2", "a", type: "fermion", label: $overline(f)$),
    edge("a", "f1", type: "antifermion", label: $chi$),
    edge("a", "f2", type: "antifermion", label: $f$),
  ),
)

#figure(
  scale(d2, 100%),
)

with the EFT scale $M_star$ given by
$
  M_star^2 tilde M_"med"/(g_chi g_q)
$
Different interactions suffer from different suppresions so one typically combines search strategies. As an example
$
  "scalar interactions" tilde m_q/M_star^3 (overline(chi) chi) (overline(q) q) tilde underbracket(v^2, "indirect" #linebreak() "suppression")"  and " overbracket(1, "direct" #linebreak() "suppression")
$
CTA upper limits on the cross section can be converted to lower limits on the EFT scale $M_star$.


