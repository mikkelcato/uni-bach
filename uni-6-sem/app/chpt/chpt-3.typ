#import "../../temp.typ": *
#show: chpt-note.with()

= Non-thermal Emission Processes
== Synchroton radiation
We have seen that accelerating charged particles will emit bremsstrahlung  as described by the Larmor formula. With _synchrotron radiation_ the idea is similar. The emitted spectrum is typically low-energy and usually due to $e^-$.

We consider some particle in a uniform magnetic field. Then the equation of motion is#footnote[Assuming $bold(E) = 0$.]
$
  dv(, t) (gamma m bold(v)) = Z e bold(v) times bold(B)
$
The LHS can be written as
$
  dv(, t) (gamma m bold(v)) &= gamma m bold(a) + underbracket(m bold(v) dv(gamma, t), -> 0)
$
where the second term vanishes since $bold(a) perp bold(v)$ in a magnetic field. We write $bold(v) = bold(v)_parallel + bold(v)_perp$ with respect to $bold(B)$ and define the _pitch angle_ $alpha$ by $abs(bold(v)) = v_perp sin alpha$ or
$
  tan alpha eq v_parallel/v_perp
$
Then $bold(v)_parallel$ does not change and we find ($a_parallel = 0$)
$
  gamma m bold(a) &= Z e v_perp abs(bold(B)) hat(bold(v)) times hat(bold(B)) \
  &= Z e abs(bold(v))abs(bold(B)) sin alpha hat(bold(v)) times hat(bold(B))
$
implying $bold(a) perp hat(bold(v))$ meaning we can equate
$
  a_perp = v_perp^2/r = (Z e abs(bold(v)) abs(bold(B)) sin alpha)/(gamma m)
$
we define the _gyro radius_
$
  r_g equiv (gamma m abs(bold(v)) sin alpha)/(Z e abs(bold(B)))
$
Then we find helical motion along $bold(B)$! We define the _cyclotron frequency_ by
$
  omega_g = v_perp/r_g = (Z e abs(bold(B)))/(gamma m)
$
We can also write (using $p = abs(bold(p)) = gamma m abs(bold(v))$)
$
  r_g & = ((p c)/(Z e)) (sin alpha)/(abs(bold(B)) c) \
      & = (R sin alpha)/(abs(bold(B)) c)
$
where we define the _rigidity_
$
  R equiv (p c)/(Z e)
$

Consider a slowly varying magnetic field i.e. $dd(B, d: Delta) \/B tilde "small"$ over $T = v_g^(-1)$. Then one can show
$
  p_perp^2/B = "const" => dd((B r_g^2), d: Delta) = 0
$
These are called _adiabatic invariants_.

We can determine the energy loss by
$
  - dv(E, t) &= ((Z e)^2 gamma^4)/(6 pi epsilon_0 c^3) (abs(a_perp)^2 + underbracket(gamma abs(a_parallel)^2, 0)) \
  &= (e^2 gamma^4)/(6 pi epsilon_0 c^3) (e^2 v^2 B^2 sin^2 alpha)/(gamma^2 m_e^2) \
  &= 2 c underbracket((e^4/(6 pi epsilon_0^2 c^4 m_e^2)), sigma_T) (v/c)^2 underbracket(B^2/(2 mu_0), U_B) gamma^2 sin^2 alpha \
  &= 2 sigma_T c U_B beta^2 gamma^2 sin^2 alpha
$
We assume $alpha$ is described by the distribution
$
  P(Omega_alpha) dd(alpha) = (2 pi sin alpha dd(alpha))/(4 pi) = 1/2 sin alpha dd(alpha)
$
with $0 <= alpha <= pi$. Then averaging we find
$
  - expval(dv(E, t)) = 4/3 sigma_T c U_B beta^2 gamma^2
$
which is quite nice.

We would now like the emitted spectrum. Within the instantaneous rest frame $S'$ of the particle then we will see the usual _dipole radiation_. We let the angle between $bold(v)$ and the observer be $phi.alt'$ and the angle between $bold(a)'$ and the observer $theta'$. We can then show that
$
  I_nu prop sin^2 theta' = cos^2 phi.alt'
$
which is the typical dipole radiation. By relativistic abberation we have
$
  sin phi.alt & = 1/gamma (sin phi.alt')/(1 + beta cos phi.alt') \
              & tilde.eq phi.alt \
              & tilde.eq^"highly relativistic" plus.minus 1/gamma
$
Meaning the radiation is strongly beamed forward in the lab frame $S$ with opening $phi.alt$.#footnote[We assume $bold(v)$ is toward us. However this does not matter since if the particle was moving away from us the emitted radiation would be weakened.] We measure a radiation pulse whenever this elongated beam sweeps past the observer. Which happens when $bold(v)$ is within an angle $plus.minus 1\/gamma$. The duration of this pulse can be found by considering two points $A$ and $B$ separated by $theta tilde 1\/gamma$. The oberver receives emission when the particle is between these points. We have
$
  t_A = R/c";  " t_B = 1/c (R-L) + L/v
$
implying
$
  dd(t, d: Delta) = L/v (1 - v/c) << L/v
$
We can write
$
  L/v = (r_g theta)/v tilde.eq 1/(gamma omega_g) equiv 1/omega_(g,"non-rel")
$
and
$
  1- beta tilde.eq 1/(2 gamma^2)
$
implying
$
  dd(t, d: Delta) tilde.eq 1/(2 gamma^2 omega_(g,"non-rel"))
$
The maximal Fourier component is then expected at the frequency
$
  nu tilde 1/dd(t, d: Delta) tilde gamma^2 nu_(g,"non-rel")
$
with
$
  nu_(g,"non-rel") = (e B)/(2 pi m)
$
we call this the _critical frequency_.

Computation of the full spectrum is more complicated and requires using the _Liénard-Wiechert_ potentials. We will not do this and simply state the result for a single $e^-$
$
  I_perp (omega) &= (sqrt(3) e^2 gamma sin alpha)/(8 pi epsilon_0 c) (F(x) + G(x)) \
  I_parallel (omega) &= (sqrt(3) e^2 gamma sin alpha)/(8 pi epsilon_0 c) (F(x) - G(x))
$
with
$
  F(x) = x integral_x^oo K_(5\/3) (z) dd(z)";  " G(x) = x K_(2\/3) (x)
$
where $x equiv omega\/omega_c = nu\/nu_c$ and
$
  nu_c equiv 3/2 gamma^2 nu_(g,"non-rel") sin alpha
$
which is complicated. These represent the spectra radiated by a single $e^-$ in the two polarisations during one period of the $e^-$ orbit. The radiated power is then found by dividing with the orbital time
$
  T_g equiv 1/nu_g = (2 pi gamma m_e)/(e B)
$
We find
$
  P_perp (omega) &= (sqrt(3) e^3 B sin alpha)/(16 pi^2 epsilon_0 c m_e) (F(x) + G(x)) \
  P_parallel (omega) &= (sqrt(3) e^3 B sin alpha)/(16 pi^2 epsilon_0 c m_e) (F(x) - G(x))
$
The total radiated power is then
$
  P(omega) = (sqrt(3) e^3 B sin alpha)/(8 pi^2 epsilon_0 c m_e) F(x)
$
above $x = 1$ i.e. $nu = nu_c$ the spectrum has a sharp cutoff.

Consider an electron density that follows a power-law in energy
$
  n equiv dd(N)/(dd(V, gamma)) = kappa gamma^(-p)
$
We want the _radiated emmissivity_
$
  j(nu) dd(nu) & = dd(E)/dd(V, t, nu) dd(nu) \
               & = (- dv(E, t)) n(gamma) dv(gamma, nu) dd(nu)
$
with $nu tilde gamma^2 nu_(g,"non-rel")$
$
  dv(gamma, nu) & = dv(, nu) (nu/nu_(g,"non-rel"))^(1\/2) \
                & = nu^(-1\/2)/(2 sqrt(nu_(g,"non-rel")))
$
Then using $nu_(g,"non-rel") tilde B$ and $- dd(E)\/dd(t) tilde gamma^2 B^2$ we find
$
  j(nu) tilde kappa nu^(-(p-1)/2) B^((p+1)/2)
$
Then the spectrum should follow a power-law with index
$
  a = (p-1)/2
$
To be correct we would need to compute
$
  j(omega) = integral_1^oo P(omega,gamma) n(gamma) dd(gamma)
$
For $n(gamma) tilde gamma^(-p)$ this guy can be computed analytically. The total _luminosity_ of some $e^-$ filled volume can then be found by
$
  L_nu = integral_V j(nu) dd(V)
$
The total _flux_ is then
$
  F_nu = L_nu/(4 pi d_L^2)
$
We define the _fractional polarisation_ by
$
  Pi (omega) &equiv (I_perp (omega) - I_parallel (omega))/(I_perp (omega)+I_parallel (omega)) = G(x)/F(x)
$
Then
$
  Pi = (integral_0^oo G(x) x^((p-3)/2) dd(x))/(integral_0^oo F(x) x^((p-3)/2) dd(x)) = (p+1)/(p+7/3)
$
Typically we observe $p tilde 2.5$ giving $Pi tilde 72 %$.

Synchrotron radiation comes alongside absorption due to $gamma$ interacting with charges in $bold(B)$ giving up some energy. This is called _synchrotron self absorption_. We define the absorption coefficient $kappa_nu$ and it is given by
$
  kappa_nu = (p+2)/(8 pi m_e nu^2) integral_0^oo dd(gamma) P(nu) n(gamma)/gamma tilde B^((p+2)/2) nu^(-(p+4)/2)
$
The flux from a homogeneous ball of radiating plasma with radius $r_b$ taking into account absorption is then
$
  F_nu = L_nu/(4 pi d_L^2) (3 u(2 r_b kappa_nu))/(2 r_b kappa_nu)
$
Assuming $2 r_b kappa_nu << 1$ i.e. optically thin we recover the previous while for $2 r_b kappa_nu >> 1$ i.e. optically thick we have
$
  F_nu tilde (nu^(-(p-1)/2))/(nu^(-(p+4)/2)) = nu^(5\/2)
$

== Inverse Compton scattering
We have _inverse Compton scattering_ whenever relativistic $e^-$ propagate through a $gamma$ gas. Whenever the $e^-$ loses energy due to this scattering we have inverse Compton scattering.

Before continuing recall _Doppler boosting_ which refers to the change in observed frequency of a moving source. We have
$
  nu = delta_D nu'";  " delta_D = gamma^(-1) (1 - beta cos theta)^(-1)
$
where $delta_D$ is the _Doppler factor_. We now consider an $e^-$ moving with $beta$ through an isotropic $gamma$ field in the lab frame $S$. Within the rest frame of the $e^-$ denoted by $S'$ all $gamma$ are incident within a narrow cone. A $gamma$ incident with $theta$ in $S$ will be incident with $theta'$ determined by
$
  tan theta' tilde.eq - 1/gamma cot theta/2
$
The $gamma$ energy in $S'$ is
$
  epsilon' = epsilon/delta_D
$
implying
$
  epsilon_"min"' tilde epsilon/(2 gamma)";  " epsilon'_"max" tilde 2 gamma epsilon
$
Two limiting cases are important:

1. $epsilon' << m_e c^2 tilde$ _Thomson regime_ where the $e^-$ loses small amounts of energy continuously.

2. $epsilon' >> m_e c^2 tilde$ _Klein-Nishina regime_ where the $e^-$ loses a large amount of energy in one scattering event.

After scattering we can compute the $gamma$ energy as
$
  epsilon.alt'_1 = epsilon.alt'/(1 + (epsilon.alt'\/m_e c^2) (1-cos theta'_1)) => epsilon.alt_1 = gamma epsilon.alt'_1 (1 + beta cos (pi - theta'_1))
$
Within the Thomson regime $epsilon'_1 tilde.eq epsilon'$ we have $ epsilon_(1,"max") tilde 2 gamma epsilon'_(1,"max") tilde 4 gamma^2 epsilon $
This process is very efficient! The energy loss for an $e^-$ accelerated by some incident radiation in the Thomson regime can be shown to be
$
  - dv(E, t) = sigma_T c u_"rad"
$
where $u_"rad"$ is the _total radiation density_.

Then in $S'$ we have
$
  - dv(E, t) = - dv(E', t') = sigma_T c u'_"rad"
$
where $u_"rad"$ is isotropic in $S$. To transform $u_"rad"$ we note that the total number of $gamma$ in some volume is invariant
$
  dd(N) (epsilon) = n(epsilon) dd(x, y, z) tilde "invariant"
$
The spacetime volume is also invariant
$
  dd(V_4) = dd(x, y, z, t) tilde "invariant"
$
which can be seen explicitly by computing#footnote[Follows immediately since $det g_(mu nu) = 1$ in Minkowski space.]
$
  dd(t, x) = matrixdet(pdv(t, t'), pdv(x, t'); pdv(t, x'), pdv(x, x')) dd(x', t')
$
Then
$
  dv(N(epsilon), V_4) = n(epsilon)/dd(t) tilde "invariant"
$
Since $dd(t)$ transforms the same as $epsilon$ by
$
  epsilon = delta_D epsilon'";  " t = delta_D t'
$
Then we obtain
$
  n(epsilon)/epsilon tilde "invariant"
$
Consider now $u'_"rad"$
$
  u'_"rad" & = epsilon' n' = delta_D^(-2) epsilon n \
           & = delta_D^(-2) u_"rad" = gamma^2 (1 - beta cos theta)^2 u_"rad"
$
By assumption $u_"rad"$ is isotropic so
$
  dv(u_"rad", Omega) = u_"rad"/(4 pi)
$
implying
$
  dv(u'_"rad", Omega) & = delta_D^(-2) dv(u_"rad", Omega) \
                      & = delta_D^(-2)/(4 pi) u_"rad"
$
We can integrate to find
$
  u'_"rad" &= u_"rad"/(4 pi) gamma^2 integral dd(Omega) (1- beta cos theta)^2 \
  &= u_"rad"/(4 pi) gamma^2 integral (1- beta cos theta)^2 sin theta dd(theta, phi.alt) \
  &= 4/3 u_"rad" (gamma^2 - 1/4)
$
The energy lost by the $e^-$ is the same as the energy gained by the $gamma$ giving
$
  dv(E, t) = 4/3 sigma_T c u_"rad" (gamma^2 - 1/4)
$
We also need to subtract the initial energy  of the $gamma$ being $sigma_T c u_"rad"$
$
  dv(epsilon_1, t) & = 4/3 sigma_T c u_"rad" (gamma^2 - 1/4 - 3/4) \
                   & = 4/3 sigma_T c u_"rad" beta^2 gamma^2 \
                   & = -dv(E, t)
$
This is the energy gained by $gamma$ in the Thomson regime. We notice this has the same structure as synchrotron radiation with the replacement $U_B -> u_"rad"$. This is the case because synchrotron radiation can be treated as scattering with virtual $gamma$ from $bold(B)$. We can also find the average $gamma$ energy after scattering by
$
  - dv(E, t) &= expval(epsilon_1) dv(n, t) \
  &= 4/3 gamma^2 beta^2 epsilon underbracket(sigma_T c u_"rad"/epsilon, "number of scattered" gamma)
$
implying
$
  expval(epsilon_1) = 4/3 gamma^2 beta^2 epsilon
$
Due to the similarity with synchrotron radiation we expect a power-law behaviour
$
  j(nu) tilde nu^(-(p-1)/2)
$

Additional notes:

1. The scattered radiation is polarised with $ Pi = (1-cos^2 theta)/(1+ cos^2 theta) $

2. Thomson scattering is very important for opacity. Assuming some medium has $e^-$ density $n_e$ then the $gamma$ density $n_gamma$ decreases as $ - dv(n_gamma, t) & = sigma_T c n_e n_y \
   -dv(n_gamma, x) & = sigma_T n_e n_gamma $
  implying
  $
    n_gamma = n_0 exp[- underbracket(integral sigma_T n_e dd(x), "optical depth" tau)]
  $
  The $gamma$ scatter and perform random walks with $ lambda = 1/(sigma_T n_e) $ This is the cause of the _cosmic microwave background_ where the opacity $tau$ drops below $1$.

  Going from a single $e^-$ to an $e^-$ cloud is complicated. We consider a $gamma$ gas with $gamma$ density $dd(n)\/dd(epsilon)$ and an $e^-$ distribution with $e^-$ density $dd(N_e)\/dd(V, gamma)$. We eventually find for $e^-$ with energy $gamma m c^2$ and $gamma$ with energy $epsilon$
  $
    dd(N)/(dd(t, epsilon_1, gamma, epsilon)) = (3 sigma_T c)/4 1/gamma^2 (1/epsilon dv(n, epsilon)) f(q = epsilon_1/(Gamma_epsilon (1-epsilon_1)); Gamma_epsilon = (4 epsilon gamma)/(m_e c^2))
  $
  with
  $
    f(epsilon_1,epsilon,gamma) = 2 q ln q + (1+ 2 q) (1- q) + 1/2 (Gamma_epsilon q)^2/(1 + Gamma_epsilon q) (1-q)
  $
  Then the total _IC emissivity_ is given by
  $
    j_(nu_1)^("IC") = h nu_1 integral dd(gamma) integral dd(epsilon) dd(N)/(dd(t, epsilon_1, gamma, epsilon)) dd(N_e)/dd(V, gamma)
  $

  To determine the energy loss in the Klein-Nishina regime we consider $Gamma_e >> 1$
  $
    - dot(gamma)^"IC" &= 1/(m_e c^2) integral dd(epsilon_1) integral dd(epsilon) (epsilon_1 - epsilon) dd(N)/(dd(t, epsilon_1, gamma, epsilon)) \
    &tilde.eq 3/8 sigma_T/(m_e c^3) integral dd(epsilon) (1/epsilon dv(n, epsilon)) (ln (4 epsilon gamma)/(m_e c^2) - 11/6)
  $
  A black body $gamma$ gas with temperature $T$ gives
  $
    - dot(gamma)^("IC", "KN") tilde T^2 ln gamma T
  $
  While in the Thomson regime we find
  $
    - dot(gamma)^("IC", "T") = 4/3 sigma_T/(m_e c) beta^2 gamma^2 U_"ph" tilde gamma^2 T^4
  $
  We see the cooling is slower in the Klein-Nishina limit!

  We can assume
  $
    dv(N_e, gamma) tilde gamma^(-p)
  $
  as before. Then we find
  $
     j_nu^"IC, T" & tilde nu^(-(p-1)/2) \
    j_nu^"IC, KN" & tilde nu^(-p)
  $
The simplest _synchrotron self Compton_ energy distribution clearly exhibits these two limits and the synchrotron spectrum with self absorption. Here the $e^-$ responsible for synchrotron radiation do inverse Compton scattering with the produced $gamma$.

== Hadronic processes
Aside from the above _leptonic_ emission models which produce high energy $gamma_"gamma"$ rays we need _hadronic_ models to explain observed $nu$. We consider binary collisions $p p -> pi X$ and _photohadronic_ interactions $N gamma -> dots$.

We start with $p p -> pi X$. Consider the distribution of produced _secondaries_#footnote[The $pi$'s in this case.] per time per kinetic energy $T_s$
$
  dot(n)_s = dd(N_s)/dd(t, V, T_s)
$
Assuming we have an isotropic $p$ distribution interacting with stationary target $p$ with density $n_p$ we have
$
  dot(n)_s (T_s) = 4 pi n_p integral_0^oo dd(T_p) J_p (T_p, Omega_p) dv(sigma_(p p -> pi), T_pi)
$
where $J_p$ is the cosmic-ray number intensity and $dd(sigma_(p p-> pi))\/dd(T_pi)$ is the inclusive cross-section for $pi$ production.#footnote[This is independent of $X$.] This depends on $n_p$ which is typically small.

We have multiple photohadronic processes of interest:
$
  N gamma & -> N pi tilde "photo-pion production" \
  N gamma & -> N e^+ e^- tilde "photo-pair production" \
  N gamma & -> N' N'' tilde "disintegration"
$
The cross-section depends on
$
  "invariant" = E_r = gamma_p E_gamma (1- beta_p mu) = E'_gamma
$
$E_gamma$ is the $gamma$ energy in the $p$ rest frame. We also have $gamma_p = E_p\/m_p c^2$ and $mu = cos theta$. We consider a general $gamma$ field with
$
  n_gamma = dd(N_gamma)/dd(E, V, Omega)
$
Then the interaction rate is
$
  dot(N)_"sc" = c integral.cont dd(Omega) integral_0^oo dd(E_gamma) (1-beta_p mu) n_gamma (E_gamma, Omega) sigma_(gamma p) (E_r)
$
The $p$ loses a fraction $K(E_r)$ of its energy on average in each interaction we refer to this quantity as the _inelasticity_. Then
$
  t_(gamma p)^(-1) equiv abs(dot(gamma)_p/gamma_p) = c integral_0^oo dd(E_gamma) underbracket(integral_0^(2 pi) dd(phi.alt) integral_(-1)^1 dd(mu), integral.cont dd(Omega)) (1- beta_p mu) n_gamma (E_gamma, Omega) sigma_(gamma p) (E_r) K (E_r)
$
We generally have $t_"sc" = dot(N)_"sc"^(-1) << t_(gamma p)$. We have $sigma_(p gamma) tilde 10^(-2) sigma_(p p)$ however typically $n_gamma >> n_p$ so
$
  lambda_(p gamma) << lambda_(p p)
$
implying $p gamma -> dots$ dominates.

The $pi$'s produced will decay after being produced through the following _secondary production channels_
$
  pi^+ & -> underbracket(mu^+, display(e^+ + overline(nu)_mu + nu_e)) + nu_mu \
  pi^0 & -> gamma + gamma
$
This leads to the observed $nu$!

We would like the spectrum of $pi^0$ decay. Within the rest frame of the $pi^0$ we have
$
  E'_gamma = m_pi/2 tilde.eq 67.5 "MeV"
$
Given $pi^0$ has spin zero the emitted $gamma$ have no preferred direction meaning their spectrum will be isotropic
$
  dv(N, Omega') = N/(4 pi)
$
We transform to the rest frame of $pi^0$ using $beta_pi$ and $gamma_pi = E_pi\/m_pi$.#footnote[Using $c = 1$.] Then
$
  E_gamma & = gamma_pi E'_gamma - beta_pi gamma_pi p'_gamma mu \
          & = (m_pi gamma_pi)/2 (1-beta_pi mu)
$
implying
$
  E_gamma^"min" = E_pi/2 (1 - beta_pi)";  " E_gamma^"max" = E_pi/2 (1+beta_pi)
$
We have
$
  dd(N)/dd(E'_gamma, Omega') = 2 delta(E'_gamma - m_pi/2)/(4 pi)
$
and
$
  "invariant" tilde 1/E_gamma dd(N)/dd(E_gamma, Omega)
$
Then $dd(N)\/dd(E_gamma, Omega)$ transforms as $E'_gamma = delta_D^(-1) E_gamma$ meaning
$
  dd(N)/dd(E_gamma, Omega) = delta(E'_gamma - m_pi/2)/(2 pi gamma_m (1- beta_pi mu))
$
We integrate over $dd(Omega) = - dd(phi.alt, mu)$ with $mu -> E_gamma$ giving
$
  dv(N, E_gamma) = 2/(beta_pi gamma_pi m_pi) Theta(E_gamma, E_pi/2 (1-beta_pi), E_pi/2 (1 + beta_pi))
$
The final $gamma$ spectrum is found numerically.
