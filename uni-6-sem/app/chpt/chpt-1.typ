#import "../../temp.typ": *
#show: chpt-note.with()

= Basics of Energy Loss
== Ionisation
When high energy particles pass through matter they can lose energy by ionising or exciting atoms and molecules in said matter. This loss of energy influences the propagation of these particles in the Universe and leads to a heating of the interstellar medium. We can use this to detect particles and measure their properties.

We start by considering the non-relativistic case where a charged particle with mass $M$, charge $Z e$ and velocity $v$ scatters of an electron at rest with some impact parameter $b$ and scattering angle $theta$. Assuming particles are hard spheres the maximal velocity of the electron after scattering is given by
$
  v_(e,"max") = (2 M)/(M+m_e) v tilde.eq 2 v
$
The fractional energy loss of the particle is
$
  (1/2 m_e (2v)^2)/(1/2 M v^2) = (4 m_e)/M <<^(M>>m_e) 1
$
This means the energy loss per collision is very small so we can assume the particles movement is unchanged due to the collision. The scattering is mediated by
$
  abs(bold(F)) = (Z e^2)/(4 pi epsilon_0) 1/r^2
$
Consider the momentum impulse given to the electron
$
  P equiv integral_(-oo)^oo bold(F) dd(t) =_"by symmetry" integral_(-oo)^oo F_perp
$
We use
$
  F_perp = (Z e^2)/(4 pi epsilon_0) 1/r^2 sin theta
$
and $v dd(t) = dd(x)$. With $x tan theta = b$ and $r sin theta = b$ we eventually find
$
  p = (Z e^2)/(2 pi epsilon_0 b v)
$
Then the kinetic energy transferred to the electron is
$
  p^2/(2 m_e) = (Z^2 e^4)/(8 pi^2 epsilon_0^2 b^2 v^2 m_e)
$
We want the average energy loss per unit length our particle will experience moving through a material with electron density $n_e$. We consider a cylindrical geometry with the impact parameter being the radial coordinate $b$ and the height being $dd(x)$. The volume element of interest is then $dd(V) = 2 pi b dd(b) dd(x)$ and we find
$
  underbracket(- dv(E, x), "energy loss") &= n_e integral_(b_"min")^(b_"max") p^2/(2 m_e) 2 pi b dd(b) \
  &= (Z^2 e^4 n_e)/(4 pi epsilon_0^2 v^2 m_e) integral_(b_"min")^(b_"max") dd(b)/b \
  &= (Z^2 e^4 n_e)/(4 pi epsilon_0^2 v^2 m_e) ln (b_"max"/b_"min")
$
To estimate $b_"max"$ we relate the collision time to the orbital period of the electron. As in this case the electron would only feel a slowly varying field and no ionisation takes place. Then
$
  tau tilde.eq (2 b_"max")/v tilde.eq (2 pi)/omega_0
$
implying
$
  b_"max" tilde.eq (pi v)/omega_0
$
For $b_"min"$ we can use
$
  p_"max"^2/(2 m_e) tilde.eq 1/2 m_e (2 v)^2
$
However, this breaks for high velocities. To do better we use
$
  dd(p, d: Delta) = m v_"max" = 2v m
$
and
$
  dd(p, d: Delta) dd(x, d: Delta) gt.tilde hbar
$
Then we find
$
  b_"min" tilde dd(x, d: Delta) >= hbar/(2 m_e v)
$
Using these guys we finally obtain
$
  - dv(E, x) = (Z^2 e^4 n_e)/(4 pi epsilon_0^2 v^2 m_e) ln ((2 pi m_e v^2)/(hbar omega_0))
$
Assuming the atom is in the ground state then
$
  I equiv (hbar omega_0)/2
$
would be the ionisation energy. In reality we need the average $overline(I)$ which is hard to compute. Given this we write
$
  - dv(E, x) = (Z^2 e^4 n_e)/(4 pi epsilon_0^2 v^2 m_e) ln ((m_e v^2)/overline(I))
$
This is the non-relativistic result.

We denote the rest frame of $M$ by $S'$ and the rest frame of the electron by $S$. We assume $bold(B) = 0$ so the only thing that happens is that $E_perp$ is enhanced by the Lorentz factor $gamma$. At the same time the collision time decreases by $gamma^(-1)$. Then
$
  tau tilde (2 b)/v tilde (2 b)/(gamma beta c)
$
implying
$
  b_"max" tilde (pi gamma beta c)/omega_0
$
And
$
  b_"min" tilde hbar/dd(p, d: Delta) = hbar/(gamma m v) = (hbar c)/(gamma m beta)
$
This is the only thing that changes since everything else cancels. Then
$
  -dv(E, x) tilde ln((2 gamma^2 beta^2 m_e c^2)/overline(I))
$
A more correct treatment gives what is called the _Bethe formula_
$
  expval(-dv(E, x)) = K underbracket(z^2, "charge") overbracket(Z/A, "of absorber") 1/beta^2 [1/2 ln((2 m_e c^2 beta^2 gamma^2 W_"max")/overline(I)^2)- underbracket(beta^2 - dd((beta gamma), d: delta)/2, "density correction" #linebreak() "due to polarization")]
$
with
$
  K &= 4 pi N_A r_e^2 m_e c^2 \
  W_"max" &= (2 m_e c^2 beta^2 gamma^2)/(1+(2 gamma m_e)\/M + (m_e\/M)^2) ->^(2 gamma m_e << M) 2 m_e c^2 beta^2 gamma^2
$
Typically this is shown against $beta gamma = p\/M c$.

== Bremsstrahlung
Ionisation losses as described above underestimate the energy loss of relativstic electrons. This additional energy loss is attributed to the acceleration of electrons in the electric field of an atom. This type of scattering amounts to flipping the picture described in the previous section. This acceleration leads to the emission of radiation called _bremsstrahlung_ and the energy loss due to this radiation is given by _Larmor's formula_
$
  - dv(E, t) = abs(dot.double(bold(p)))^2/(6 pi epsilon_0 c^3) = (q^2 abs(bold(a)_0)^2)/(6 pi epsilon_0 c^3)
$
where $bold(p) = q bold(d)$ is the electric dipole moment and $bold(a)_0$ is the proper acceleration. We need to rewrite this in the lab frame. Before seeing how this is done consider
$
  integral_(-oo)^oo dv(E, t) dd(t) &= q^2/(6 pi epsilon_0 c^3) integral_(-oo)^oo abs(dot(bold(v))(t))^2 dd(t) \
  &=^"Parseval" q^2/(6 pi epsilon_0 c^3) integral_(-oo)^oo abs(dot(bold(v))(omega))^2 dd(omega) \
  &= integral_0^oo q^2/(3 pi epsilon_0 c^3) abs(dot(bold(v))(omega))^2 dd(omega) \
  &= integral_0^oo I(omega) dd(omega) \
  &= "total emitted radiation"
$
Where we have defined the _emitted spectrum_
$
  I(omega) equiv q^2/(3 pi epsilon_0 c^3) abs(dot(bold(v))(omega))^2
$
This is the quantity we would like to find. We proceed by considering the four-acceleration
$
  a_mu &= gamma vecrow(c dot(gamma), pdv(, t) (gamma bold(v))) \
  &= vecrow(gamma^4 c ((bold(v) dot bold(a))/c^2), gamma^2 bold(a) + gamma^4 bold(v) ((bold(v) dot bold(a))/c^2))
$
with $bold(a)$ and $bold(v)$ measured in the lab frame $S$. In the rest frame $S'$ we just have
$
  a'_mu = vecrow(0, bold(a)_0)
$
We want to relate $bold(a)_0$ with $bold(a)$ and $bold(v)$. Using $a_mu a^mu = a'^mu a'_mu$ we find
$
  abs(bold(a)_0)^2 = gamma^4 (abs(bold(a))^2 + gamma^2 ((bold(v) dot bold(a))/c)^2)
$
Then
$
  -dv(E, t) = (q^2 gamma^4)/(6 pi epsilon_0 c^3) [abs(bold(a))^2 + gamma^2 ((bold(v) dot bold(a))/c)^2]
$
where we have used
$
  dv(E, t) = underbracket(dv(E, E'), gamma) dv(E', t') underbracket(dv(t', t), gamma^(-1)) = dv(E', t') tilde "Lorentz invariant"
$
Writing $bold(a) = bold(a)_parallel + bold(a)_perp$ we obtain
$
  [dots] &= abs(bold(a)_parallel)^2 + abs(bold(a)_perp)^2 + gamma^2 (v a_parallel)^2/c^2 \
  &= abs(bold(a)_perp)^2 + gamma^2 abs(bold(a)_parallel)^2
$
So we are left with
$
  - dv(E, t) = (q^2 gamma^4)/(6 pi epsilon_0 c^3) (abs(bold(a)_perp)^2 + gamma^2 abs(bold(a)_parallel)^2)
$
To determine the acceleration we use $bold(F) = q bold(E)$ to find
$
  a_parallel & = dot(bold(v))_x = - (Z e E_x)/m_e \
      a_perp & = dot(bold(v))_z = - (Z e E_z)/m_e
$
Determining $E_x$ and $E_z$ is very annoying due to relativitistic business. We eventually find
$
  a_parallel &= dot(bold(v))_x = (gamma Z e^2 v t)/(4 pi epsilon_0 m_e [b^2 + (gamma v t)^2 ]^(3\/2)) \
  a_perp &= dot(bold(v))_z = (a_parallel b)/(v t)
$
We need the Fourier transform of these
$
  dot(bold(v))_x = (Z e^2)/(4 pi epsilon_0 m_e) 1/sqrt(2 pi) integral_(-oo)^oo (gamma v t)/([b^2 + (gamma v t)^2]^(3\/2)) e^(i omega t) dd(t)
$
Letting $b xi equiv gamma v t$ we find
$
  dot(bold(v))_x = (Z e^2)/(4 pi epsilon_0 m_e) 1/sqrt(2 pi) 1/(b v gamma) 2 i y K_0 (y)
$
where $K_0$ is the modified Bessel function and we have defined
$
  y equiv (omega b)/(gamma v)
$
Similarly we find
$
  dot(bold(v))_z = (Z e^2)/(4 pi epsilon_0 m_e) 1/sqrt(2 pi) 1/(b v) 2 y K_1 (y)
$
The emitted spectrum is then
$
  I(omega) &= e^2/(3 pi epsilon_0 c^3) (abs(a_parallel (omega))^2 + abs(a_perp (omega))^2) \
  &= e^2/(3 pi epsilon_0 c^3) 1/(2 pi) ((Z e^2)/(4 pi epsilon_0 m_e) 1/(b v))^2 ((omega b)/(gamma v))^2 \
  &#h(1em) times 2^2 [1/gamma^2 K_0^2 ((omega b)/(gamma v)) + K_1^2 ((omega b)/(gamma v))]
$
For $y << 1$ or large $v$ we have
$
  K_0 (y) & -> ln y \
  K_1 (y) & -> 1/y
$
Meaning the $perp$ contribution becomes constant and completely dominates the $parallel$ contribution! We find
$
  I (omega) &tilde.eq^(y << 1) (2^2 e^2)/(3 pi epsilon_0 c^3) 1/(2 pi) ((Z e^2)/(4 pi epsilon_0 m_e) 1/(b v))^2 \
  &tilde.eq I_0/(2 pi) 1/(b^2 v^2)
$
where we have defined
$
  I_0 equiv (Z^2 e^6)/(12 pi^3 epsilon_0^3 c^3 m_e^2)
$
The last thing we need is to integrate over $b$. We have
$
  I(omega') &= integral_(b'_"min")^(b'_"max") 2 pi b' underbracket(gamma n, "lab frame") v I_0/(2 pi b'^2 v^2) dd(b') \
  &= (I_0 n gamma)/v ln (b'_"max"/b'_"min")
$
To estimate $b_"min"$ we consider
$
  hbar omega =^"maximal energy" gamma m_e v^2
$
with $tau tilde 2 b_"min"\/v$ and $omega' tilde.eq 2 pi\/tau$ we find
$
  b_"min" = hbar/(m_e v)
$
To estimate $b_"max"$ we consider when the charge of the nucleus gets screened by the electron cloud of the atom. By Thomas-Fermi the potential of the nucleus is
$
  V (r) tilde.eq (Z e^2)/(4 pi epsilon_0 r) e^(-r\/a)
$
with
$
  a equiv (1.4 a_0)/z^(1\/3)
$
The screening distance is then $b_"max" tilde a$. With everything we find
$
  I(omega) & = dd(E)/(dd(omega, t)) \
           & = dd(E)/(dd(omega', t)) dv(omega', omega) \
           & = I(omega') gamma^(-1) \
           & = (I_0 n)/v ln ((1.4 a_0 m_e v)/(z^(1\/3) hbar))
$
Then the total loss is
$
  - dv(E, t) = integral_0^(E\/hbar) I(omega) dd(omega)
$
where we integrate until the electron has given up all its energy in one collision. With $v tilde.eq c$ we find
$
  -dv(E, t) = (I_0 n)/(hbar c) E underbracket(ln (192/z^(1\/3)), "mainly screening")
$
A more correct treatment by Bethe and Heitler gave
$
  -dv(E, t) = (z(z+1.3) e^6 n)/(16 pi^3 epsilon_0^3 m_e^2 c^4 hbar) E [ln((183)/z^(1\/3))+1/8]
$
We see the energy loss is exponential. To quantify this we define the radiation length $X_0$ by
$
  -dv(E, x) = - dv(E, t) 1/c = E/X_0
$
The above is quite good for $z gt.tilde 4$ where
$
  1/X_0 tilde.eq 4 (hbar/(m_e c))^2 z(z+1.3) alpha^3 n ln(183/z^(1\/3))
$
with $alpha$ being the fine structure constant. We have found
$
  expval(-dv(E, x)) &tilde.eq E/X_0 \
  &eq^"total loss" underbracket(a(E), "ionisation") + underbracket(b(E) E, "bremsstrahlung" #linebreak() + "others")
$
Sometimes the radiation length is written using $xi_0 = rho X_0$.

Consider the number of radiated photons
$
  dd(N)/(dd(t, omega)) = I(omega)/(hbar omega)
$
since $I(omega)$ is constant we get either few very high energy photons or many low energy photons.

== Photon energy losses
We have
$
  N = N_0 e^(-x\/lambda)
$
with $lambda$ being the mean free path length
$
  1/lambda = n sigma
$
with $sigma$ being the cross-section which depends on the given process we are considering e.g. the photoelectric effect, Rayleigh scattering, Compton scattering, etc. The process we care about and which dominates at large energies ($>1 "MeV"$) is pair-production
$
  gamma + gamma -> e^+ + e^-
$
The cross-section for this process is
$
  sigma & = 7/9 A/N_A 1/(rho X_0)
$
so
$
  lambda = 9/7 X_0
$

== Electromagnetic showers
When high energy $gamma$ or $e^(plus.minus)$ enter an absorber they will initiate an _electromagnetic cascade_. The $e^-$ will produce $gamma$ through bremstrahlung which will then pair-produce $e^(plus.minus)$ again producing $gamma$ through bremstrahlung. The $gamma$ just pair-produce $e^(plus.minus)$ directly. The difference in showers initiated by the two is the point of first interaction since $lambda_gamma > lambda_e$.

As a simple model we assume the incoming particle has energy $E_0 >> E_c$ with $E_c$ being the _critical energy_. At energies above $E_c$ the energy loss is dominated by bremsstrahlung and below it is dominated by ionisation. We assume $e^(plus.minus)$ travels $X_0$ before losing half its energy to a $gamma$ and $gamma$ travels $X_0$ before pair-producing $e^(plus.minus)$ with each carrying half the energy. We define
$
  t = X/X_0";  " epsilon = E/E_c
$
An example of a shower is:

#diagram(
  node((0, 0), name: <0>),
  edge(<0>, <1>, "wave", label: $gamma$),
  node((0, 1), name: <1>),
  node((1.5, 2), name: <2A>),
  node((-1.5, 2), name: <2B>),
  edge(<1>, <2A>, "-", label: $e^-$),
  edge(<1>, <2B>, "-", label: $e^+$),

  node((2.5, 3), name: <3D>),
  node((0.5, 3), name: <3C>),
  node((-0.5, 3), name: <3B>),
  node((-2.5, 3), name: <3A>),

  edge(<2A>, <3D>, label: $e^-$),
  edge(<2A>, <3C>, "wave", label: $gamma$),
  edge(<2B>, <3B>, "wave", label: $gamma$),
  edge(<2B>, <3A>, label: $e^+$),
)


After $t$ radiation lengths we have $N(t) &= 2^t$ and roughly equal amounts of $gamma$ and $e^(plus.minus)$. The average energy is then
$
  E(t) = E_0/2^t
$
The cascading stops when $E = E_c$ meaning
$
  E_0/(2^(t_"max")) = E_c => t_"max" = 1/(ln 2) ln E_0/E_c
$
implying
$
  N_"max" = E_0/E_c equiv y
$
A more sophisticated model gives
$
  dv(E, t) = (E_0 b (b t)^(a-1) e^(-b t))/(Gamma(a))
$
and
$
  t_"max" = (a-1)/b = ln E_0/E_c - c_j
$
with $c_e = 1$ and $c_gamma = 1/2$. Another quantity sometimes used is the _shower age_
$
  s = (3 t)/(t+2t_"max")
$
so $s_"first" = 0$, $s(t_"max") = 1$, and $s_"death" = 3$. The transversal development of the shower scales $tilde$ with the _Moliére radius_
$
  R_M tilde.eq (21 "MeV")/E_c X_0
$
Aside from these electromagnetic showers then _hadronic showers_ are also possible. Here most particles end up being $pi$ with $tilde$ a third $pi^0$ decaying into $gamma$.

== Cherenkov radiation
_Cherenkov radiation_ is caused by the charged particles propagating in a medium with $beta > n^(-1)$. When this happens a coherent waveform is created in similar fashion to a _Mach cone_ for sound. We call this cone the _Cherenkov cone_ and it is characterised by the _Cherenkov angle_
$
  cos theta_c equiv 1/(n beta)
$
The energy loss due to Cherenkov radiation is incredibly small and therefore negligible. However it is very important for particle detection as we will see later. The emitted spectrum for Cherenkov radiation follows
$
  dd(N)/dd(x, lambda) tilde.eq (2 pi alpha z^2)/lambda^2 sin theta_c
$
importantly the emitted photons are optically visible.
