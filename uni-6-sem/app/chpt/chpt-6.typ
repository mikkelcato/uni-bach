#import "chpt-temp.typ": *
#show: chpt-note.with()

= Galactic Sources of CRs
== SNR
The best candidates for energies up to the _knee_ are supernova remnants or SNRs. We consider stars with $M gt.tilde 8 M_dot.circle$. These will experience _core-collapse_ when fusion can no longer keep them stable.#footnote[These massive stars can not create energy through the fusion of iron thereby collapsing.] When collapse stops in-falling material will bounce back forming a shockwave. The now dense core quickly cools through $nu$-emission, with some $nu$ being absorbed by the shell driving the shockwave.#footnote[$99%$ of the released $E_"grav"$ is released through these $nu$.] After the dust has settled we are left with a dense _neutron core_ or neutron star surrounded by an $tilde$ spherical shell.#footnote[This is caused by inverse $beta$-decay.] Less massive stars will experience a similar core-collapse with less energy. These leave behind _white dwarfs_ since they do not have sufficient mass to form a neutron core. We observe these core-collapse supernovae through $nu$-emissions and GWs since these events are quite extreme.#footnote[Aside from the visible explosion.]

We now show these are candidates due to $1^"st"$ order Fermi acceleration. We want particles to gain energy upon crossing the shockwave, we need enough SNe with sufficient energy in their ejecta to account for the observed CR density, and the maximal energy should be high enough, per Hillas' criterion to accelerate CRs up to the knee. We also expect $alpha tilde 2$ which we showed is possible in the last section. We estimate $v_"ejecta"$. Assume a star is a sphere of constant density $rho$. Then
$
  M = 4/3 pi r^3 rho => dv(M, r) = 4 pi r^2 rho
$
implying
$
  dd(U_"grav") & = (G M dd(M))/r \
               & = (16 pi^2 G r^4 rho dd(r))/3
$
and
$
  dd(U, d: Delta) & = integral_R^(R_"NS") dv(U_"grav", r) dd(r) \
                  & tilde.eq (3 G M_"NS"^2)/(5 R_"NS")
$
Assuming $M_"NS" = 1.4 M_dot.circle$ and $R_"NS" tilde 10 "km"$ we find
$
  dd(U, d: Delta) tilde 10^53 "erg"
$
This is the amount of liberated energy. We assume $1%$ drive the shockwave while the rest escapes in the form of $nu$
$
  K tilde 10^51 "erg"
$
implying
$
  eta = v_s/c tilde 1/c sqrt((2 K)/M) tilde 0.02
$
Which is non-relativistic but still much larger than typical $v_s$. The shockwave expands and sweeps up the ISM with
$
  rho_"ISM" tilde 1.6 times 10^(-24) "g"/"cm"^3
$
This happens until $rho_"SN" tilde rho_"ISM"$ after which the shockwave slows down significantly. This is called the _Sedov phase_. This occurs when
$
  rho_"SN" = (3M_"SN")/(4 pi R_"SN"^3) = rho_"ISM"
$
With $M_"SN" tilde 10 M_dot.circle$ we find
$
  R_"SN" tilde 5 "pc"
$
and
$
  tau_"SN" tilde R_"SN"/v_s tilde cal(O)(10^3 "yr")
$
This is a long time which can be used to accelerate particles! The typical SNe is currently in the Sedov phase. An important expection is the _Crab nebula_ which is being driven continuously by a central _pulsar_. We observe $3$ core-collapse SNe per century so
$
  L = K f_"SN" tilde 10^42 "erg"/"s"
$
which we can compare with the total CR luminosity. We have
$
  rho_"CR" tilde 1 "eV"/"cm"^3
$
and with $V_"Galactic" tilde 10^66 "cm"^3$ and $tau_"esc" tilde 10^7 "yr"$ we find
$
  L_"CR" tilde (rho_"CR" V_"Galactic")/(tau_"esc") tilde 10^40 "erg"/"s"
$
Then $1%$ of the energy released in SN ejecta is enough to explain the observed CR density. The maximal energy is given by Hillas' criterion
$
  E_"max" & = e B c R_"SN" Z beta_u \
          & tilde^(B = 4 mu"G") 300 Z "TeV"
$
Then SNRs can accelerate CRs up to the knee. Then SNRs are nice candidates.

When observing SNRs we see a $pi^0$ bump at low energies. This is expected when we assume CRs are accelerated since this leads to hadronic interactions. However, the identification of PeVatron sources is weird. The detection of $gt.tilde 100 "TeV"$ $gamma$-rays would suggest the presence of freshly accelerated CR protons with PeV energies. These $gamma$-rays would be produced in the $pi^0$ decay produced in $2 p -> pi X$ while leptonic interactions would be suppressed at these energies. Whenever a $pi^0$ bump is observed SNRs are generally the preferred candidate. However, PWNe can also lead to $gt.tilde 100 "TeV"$ $gamma$-rays.

== White dwarfs
// https://indico.global/event/17334/registrations/4202/?token=e809a16a-0b6b-4a82-adc6-f98d2ccba9eb
We would like to understand compact objects such as NS since these are observed at high energies and drive _pulsar wind nebulae_. They also emit GW by existing and through mergers as we will see.#footnote[BH binaries can also form _micro-quasars_ etc.]

To understand these objects we start with the simplest being white dwarfs. These WD stop collapsing due to two quantum mechanical effects since their densities becomes immense. Those being Heisenberg's uncertainty principle and Fermi's exclusion principle. When the intermolecular spacing $dd(x, d: Delta)$ becomes small the spread in momenta $dd(p, d: Delta)$ becomes very small thereby providing _degeneracy pressure_ to withstand further collapse. We would like to estimate this pressure. We have
$
  dd(x, d: Delta) dd(p, d: Delta) tilde p dot d tilde hbar
$
with each $e^-$ occupying a cube of volume $d^3$ implying
$
  n tilde 1/d^3
$
and
$
  p tilde hbar/d tilde hbar n^(1\/3)
$
Assuming the particles are non-relativistic we have
$
  E = p^2/(2 m_e) tilde (hbar^2 n^(2\/3))/(2 m_e)
$
This is basically the _Fermi energy_. A proper treatment gives
$
  E_F = hbar^2/(2 m) (3 pi^2 n)^(2\/3)
$
Consider an ionized hydrogen gas with
$
  rho tilde m_p/d^3
$
and
$
  epsilon & = (E rho)/m_p \
          & = hbar^2/(2 m_e) (rho/m_p)^(5\/3)
$
We assume the gas is ideal and adiabatic
$
  p tilde rho^gamma
$
where
$
  gamma = (pdv(ln p, ln rho))_S
$
This follows from
$
  p V^gamma = "constant"
$
implying
$
  dd(p)/p = - gamma dd(V)/V
$
and by definition of $rho$
$
  dd(rho)/rho = - dd(V)/V
$
Also,
$
  p = (gamma-1) epsilon
$
which follows from $dd(U) = -p dd(V)$. Then
$
  p tilde epsilon tilde rho^(5\/3)
$
implying
$
  gamma = 5/3
$
and
$
  p = 2/3 epsilon tilde hbar^2/(3 m_e) (rho/m_p)^(5\/3)
$
A proper treatment gives
$
  p = (3 pi^2)^(2\/3)/5 hbar^2/m_e (rho/(mu_e m_u))^(5\/3) tilde rho^(5\/3)
$
which is valid for all compositions. We would like to know what happens to the equation of state when $e^-$ become relativistic, and when this occurs. We can write
$
  Delta p tilde^"relativistic" m_e c
$
implying
$
  rho_"crit" = m_p/d^3 tilde m_p ((m_e c)/hbar)^3 tilde 3 times 10^10 "kg"/"m"^3
$
A proper treatment gives#footnote[We require $hbar k_F = m_e c$]
$
  rho_"crit" = (m_u mu_e)/(3 pi^2) ((m_e c)/hbar)^3 tilde 10^9 mu_e "kg"/"m"^3
$
We can determine the equation of state as before. However, now we use
$
  E = sqrt(p^2 c^2 - (m c^2)^2) tilde p c tilde (hbar c)/d tilde hbar c n^(1\/3)
$
implying
$
  epsilon tilde E/d^3 tilde (hbar c)/d^4
$
with
$
  rho tilde m_p/d^3
$
we find
$
  epsilon tilde hbar c (rho/m_p)^(4\/3)
$
implying
$
  gamma = 4/3
$
and
$
  p = epsilon/3 tilde (hbar c)/3 (rho/m_p)^(4\/3)
$
A proper treatment gives
$
  p = (3 pi^2)^(1\/3)/4 hbar c (rho/(mu_e m_u))^(4\/3) tilde rho^(4\/3)
$
These two equations are the equations of state for degenerate $e^-$ gasses, including those within WDs and NS.

We can balance the pressure of an $e^-$ gas like those above with $U_"grav"$ to determine the maximal allowed mass of a WD.#footnote[or NS.] We apply the _virial theorem_ to get a rough estimate
$
  2 U = U_"grav"
$
with
$
  U & = epsilon V tilde hbar c underbracket(((rho V)/m_p)^(4\/3), (M\/m_p)^(4\/3)) underbracket(V^(-1\/3), tilde R^(-1)) \
  U_"grav" & = 1/2 (G M^2)/R
$
implying
$
  M_"Ch" tilde ((hbar c)/G)^(3\/2) 1/m_p^2 tilde 2 M_dot.circle
$
which we call the _Chandrasekhar mass_. We can define a _gravitational fine structure constant_
$
  alpha_G = (G m_p^2)/(hbar c) tilde 10^(-39)
$
and write
$
  M_"Ch" tilde m_p alpha_G^(-3\/2)
$
implying a WD consists of $tilde 10^60$ protons. A proper treatment gives
$
  M_"Ch" tilde 1.4 M_dot.circle
$
We can now estimate the size of a WD
$
  R_"WD" tilde (M_"Ch"/rho_"crit")^(1\/3) = \u{019B}_e/sqrt(alpha_G) tilde R_dot.circle/100
$
where $\u{019B}_e$ is the _Compton wavelenght of $e^-$_
$
  \u{019B}_e = hbar/(m_e c)
$

== NS and pulsars
We now consider NS. As core-collapse occurs we have stages of increasing $rho$. The following occur in order

1. $e^-$ gas becomes degenerate $p tilde rho^(5\/3)$.

2. $e^-$ gas becomes relativistic $p tilde rho^(4\/3)$ at $rho_"crit"$.

3. At $rho tilde 10^10 "kg"/"m"^3$ we have _neutronisation_. The $e^-$ have enough energy to perform inverse $beta$-decay: $p + e^- -> n + nu_e$ $ E_"tot" = gamma m_e c^2 >= (m_n - m_p) c^2 $ The neutrons cannot perform $beta$-decay as all $e^-$ states are occupied. This occurs for $ p_F >= sqrt(E_"tot"^2\/c^2 - m_e^2 c^2) $

4. At $rho tilde 10^14 "kg"/"m"^3$ we have the _neutron drip process_. The heavy nuclei become enriched with neutrons from inverse $beta$-decay making them unstable.

5. At $rho tilde 10^17 "kg"/"m"^3$ we have a degenerate $n$ gas. This occurs after all $p -> n$. All equations for NS are the same as those for WD "derived" above.

We note that the magnetic field of a star is conserved during core-collapse. The surface magnetic field $B$ scales as
$
  B tilde r^(-2)
$
implying NS can have very large magnetic fields. We also observe NS that rotate rapidly. These are called _pulsars_. We observer them as very stable, short periods of pulses and polarized radio emission. The idea is that radio emission is along a magnetic axis which is misaligned with the rotational axis. This would lead to pulses as we observe.

We can estimate the surface magnetic field using the pulse period $P$ and the _spin-down rate_ $dot(P)$. The rotation slows down with _braking index_ $n$
$
  dot(Omega) = - kappa Omega^n
$
where $Omega$ is the angular frequency and $n$ provides information about which processes act to slow the rotation. The important braking mechanism for our purposes is _magnetic braking_. The idea is that a magnetic dipole forms at an angle with respect to the rotational axis. This leads to a varying dipole moment and thereby radiation. We can use the _Larmor formula_ again
$
  - dv(E, t) = (mu_0 abs(dot.double(p)_m)^2)/(6 pi c^3)
$
a rotating dipole has
$
  p_m = p_(m 0) sin Omega t
$
Then
$
  expval(-dv(E, t)) = (mu_0 Omega^4 p_(m 0)^2)/(12 pi c^3)
$
The rotational energy is given by
$
  E = 1/2 I Omega^2
$
implying
$
  -dv(E, t) = - I dot(Omega) Omega
$
With the previous
$
  dot(Omega) tilde Omega^3
$
so $n = 3$. We can find
$
  n = (Omega dot.double(Omega))/dot(Omega)^2 = (nu dot.double(nu))/dot(nu)^2 = 2 - (dot.double(P) P)/dot(P)^2
$
Assuming $n$ is constant during the pulsars lifetime $tau$ we can find
$
  tau = (P)/((n-1) dot(P)) =^(n=3) P/(2 dot(P))
$
With
$
  B_s tilde (mu_0 p_(m 0))/(4 pi R^3)
$
we can also find
$
  B_s = ((3 mu_0 c^3 M)/(80 pi^3 R^4))^(1\/2) sqrt(P dot(P))
$
where we used
$
  I = 2/5 M R^2
$
We typically plot pulsars in a _$P$-$dot(P)$ diagram_ which can be treated as the _Hertzsprung-Russell diagram_ for pulsars.

The rotating magnetic field of any pulsar induces an electric field much stronger the gravitational field. This leads to charges being ripped of the surface. These charges will rearragne such that no Lorentz force acts on them
$
  0 = bold(E) + (bold(Omega) times bold(r)) times bold(B)
$
The distribution of charges is given by
$
  rho = (div bold(E))/(4 pi) tilde - (bold(Omega) dot bold(B))/(2 pi c) = (B_s)/(2 pi c) (R/r)^3 (3 cos^2 theta - 1)
$
The charges cannot co-rotate with field lines beyond the _light cylinder_ where $R_"LC" Omega = c$. Any field lines that would close beyond $R_"LC"$ become open and the charges can escape. The last closed field line has angle
$
  sin theta_"PC" = sqrt((Omega R)/c)
$
This defines the _polar cap region_ on the surface from which particles can escape to infinity. There is a large potential difference between the pole and the edge of the polar cap region given by
$
  - Delta Phi tilde 1/2 B_0 R theta_"PC" tilde 10^13"V"
$
When observing pulsars we see radiation across the entire spectrum. The produced emission is believed to originate due to charged particles being accelerated along field lines i.e. through synchrotron radiation. This is also believed to cause cascades like in our atmosphere. There is likely also some IC emission.

== PWNe
When a pulsar is embedded in a nebula from some SN explosion then the spin-down luminosity from said pulsar can drive emission from the nebula. The $e^-$ from the pulsar can be accelerated to $"PeV"$ energies by the SNR shockwave. These are referred to as _pulsar wind nebulae_.#footnote[An example of such an object is the _Crab nebula_.] Also, see _$"TeV"$ halos_.
