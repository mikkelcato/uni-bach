#import "chpt-temp.typ": *
#show: chpt-note.with()

= Fermi Acceleration
We deal with plasmas in space and these typically have high conductivity. This is problematic since any electric fields are short-circuited and therefore unable to efficiently accelerate particles. However, we do observer very energetic charged particles. The mechanism we will now describe accounts for this and is called _Fermi acceleration_. The idea is that particles are able to scatter of magnetic irregularities _frozen_ into clouds in the ISM.

== Magnetic mirrors
The above process is sometimes referred to as _magnetic mirrors_. We assume a charged particle travels in a non-uniform magnetic field with $v << c$. Then the particle will not only follow helical motion around the field lines. Because of non-uniformity the particle feels a time-varying magnetic field in its instantaneous rest frame. This induces an electric field
$
  curl bold(E) = - pdv(bold(B), t)
$
The induced $bold(E)$ will act to change $v_perp$ and also the energy as
$
  Delta E_perp & = 1/2 m Delta v_perp^2 \
               & = integral.cont q bold(E) dot dd(bold(l)) \
               & = q integral.cont curl bold(E) dd(bold(S)) \
               & = -q integral.cont pdv(bold(B), t) dd(bold(S))
$
We assume the magnetic field changes by a small amount during one period of one Larmor radius $r_g$
$
  T_g = (2 pi)/omega_g
$
with frequency
$
  omega_(g,"non-rel") = (q B)/m
$
We approximate
$
  -pdv(bold(B), t) tilde.eq (Delta B)/T_g
$
Then
$
  Delta E_perp & tilde.eq q Delta B omega_g/(2 pi) pi r_g^2 \
               & = 1/2 m v_perp^2 (Delta B)/B \
               & = E_perp (Delta B)/B
$
or
$
  (Delta E_perp)/E_perp = (Delta B)/B
$
implying
$
  E_perp/B tilde "constant"
$
With $v_perp = v sin alpha$ we find
$
  (sin^2 alpha)/B tilde "constant"
$
Then if $B$ changes the pitch angle $alpha$ must change accordingly thereby reflecting the particle.

== $2^"nd"$ order Fermi acceleration
We consider particles accelerated by collisions with gas clouds acting as magnetic mirrors. We denote the particles mass and velocity by $m$ and $bold(v)$ respectively and the clouds mass and velocity by $M$ and $bold(u)$ respectively. We will assume $M >> m$ and $u << v$. The velocity $v'$ after scattering is simply
$
  bold(v)' & = ((m-M) bold(v) + 2 M bold(u))/(m+M) \
        v' & tilde.eq -v plus.minus 2 u
$
With $-$ corresponding to the _head-on_ case. The energy after scattering is
$
  E' & = 1/2 m (-v plus.minus 2 u)^2
$
implying
$
  Delta E_(1,2) & = 1/2 m (minus.plus u v + 4 u^2) \
                & = minus.plus 4 u/v E +4 (u/v)^2 E \
                & tilde.eq plus.minus 4 u/v E
$
We expect the following rates
$
  f_(1,2) = (v plus.minus u)/cal(l)
$
with $cal(l)$ being the particle mean free path. Then
$
  (Delta E)/(Delta t) & = f_1 E_1 + f_2 E_2 \
                      & = (8 u^2)/(cal(l) v)E \
                      & equiv tau''_F E
$
Where we have defined the time constant $tau''_F$ which is second order in $u$. We can compute
$
  tau''_F tilde 10^7 "years"
$
So the process is slow and competes with ionisation losses meaning it is difficult to accelerate particles with low energies to sufficiently high energies. This is related to the _injection problem_.

== $1^"st"$ order Fermi acceleration
We can obtain higher efficiency by considering particles moving between mutually approaching clouds or if particles cross a shock front. We consider two frames: $S$ the observer reference frame and $S'$ the rest frame of a magnetic cloud traveling with velocity $u$ along the $x$-axis. We define the _bulk_ Lorentz factor as
$
  Gamma = 1/sqrt(1-(u\/c)^2)
$
Then
$
    E' & = Gamma (E + u p_x) \
  p'_x & = Gamma (p_x + u/c^2 E)
$
We assume the collision is elastic so
$
  E'^* = E'";  " p'^*_x = - p'_x
$
with $*$ indicating after scattering. Then
$
  E^* & = Gamma (E'^* - u p'^*_x) \
      & = Gamma (E' + u p'_x) \
      & = Gamma^2 [E + u p_x + u p_x + u^2/c^2 E] \
      & = Gamma^2 E [1 + (2 u p_x)/E + u^2/c^2]
$
We use
$
  v_x/c^2 = p_x/E
$
with $v_x = v cos theta$ and since $u << c$ we have
$
  Gamma^2 tilde.eq 1 + (u/c)^2
$
implying
$
  Delta E & = E^* - E \
          & tilde.eq 2E [(u v cos theta)/c^2 + u^2/c^2] \
          & tilde.eq^(v tilde c) (2 E u cos theta)/c
$
Then the particle gains energy if the collision is head-on
$
  -pi/2 < theta < pi/2
$
We average $cos theta$ over this range using $dd(Omega) = sin theta dd(phi.alt, theta)$ and a probability distribution $Delta E tilde cos theta$
$
  expval(cos theta) &= (integral_0^(2 pi) dd(phi.alt) integral_0^(pi\/2) dd(theta) cos^2 theta sin theta)/(integral_0^(2 pi) dd(phi.alt) integral_(-pi\/2)^(pi\/2) dd(theta) cos theta sin theta)
$
The $phi.alt$ integration cancels and we substitute $mu = cos theta$
$
  expval(cos theta) & = (integral_0^1 mu^2 dd(mu))/(integral_0^1 mu dd(mu)) \
                    & = 2/3
$
Averaging over all angles we find $expval(cos theta) = 0$ recovering second order Fermi acceleration. However, for $cos theta> 0$ we find
$
  expval(Delta E) & = (4 u)/(3 c) expval(E) equiv eta expval(E) \
      expval(E^*) & = (1 + (4 u)/(3 c)) expval(E) equiv scr(B) expval(E)
$
We could obtain the same result by only considering head-on collisions in the previous section
$
  (Delta E)/(Delta t) & = f_1 E_1 \
                      & = 4 u/v E (u+ v)/cal(l) \
                      & tilde.eq (4 u)/cal(l) E \
                      & equiv E/tau'_F
$
With
$
  tau'_F = cal(l)/(4 u) = 2 tau''_F u/c << tau''_F
$
We could imagine scenarios where particles only experience head-on collisions such as due to diffusive shocks. We imagine a shock wave propagating through a medium. This shock wave is characterised by a discontinuous change in characteristics of the medium such as density, velocity and temperature. We say the shock is supersonic when the shock velocity $v_s$ is larger than the speed of sound in the medium. The shock velocity is by definition much larger than the thermal velocity of particles $v_s >> v_1$. Within the downstream rest frame $S_"down"$ scattering ensures the particle distribution is isotropic. However, the upstream gas as a uniform velocity $v_2$ in $S_"down"$. Then a particle crossing the shock wave will scatter head-on with the upstream gas since
$
  v_2 = 3/4 v_s
$
The same thing happens within the upstream rest frame $S_"up"$ since the downstream gas appears to have a uniform velocity $v_1$ in this frame. Then a particle crossing the shock wave will scatter head-on with the downstream gas since
$
  v_1 = - 3/4 v_s
$
These numbers are for a monoatomic gas where $gamma = 5/3$ since
$
  v_1/v_2 = (gamma+1)/(gamma-1)
$
The shock velocity is $v_s$ in the lab frame $S_"obs"$. Then in the rest frame of the shock $S_"shock"$ we have
$
  v_1 = -v_s";  " v_2 = -v_s/4
$
Then the numbers above follow.

With
$
  u = 3/4 v_s
$
we have
$
  eta = v_s/c
$
and
$
  scr(B) = 1 + v_s/c
$

== Emitted spectrum
We have shown how Fermi acceleration could accelerate particles to high energies. We still need to show that these particles produce the power-law particle distribution we observe. We will assume a simple model.

We consider a particle gaining an amount of energy $scr(B)$. The particle could escape the acceleration region. We denote the probability this happens by $P_"esc"$ so the probability that the particle keeps accelerating is $ P = 1 - P_"esc" $ Then after $k$ cycles we have
$
  N = N_0 P^k
$
particles remaining with energy larger than $E = E_0 scr(B)^k$. We can find
$
  N(>=E) = N_0 (E/E_0)^(ln P\/ln scr(B))
$
This quantity is related to the spectrum $N(E)$ by
$
  N(>= E) = integral_E^oo N(E') dd(E')
$
implying
$
  N(E) tilde (E/E_0)^(ln P \/ ln scr(B) - 1)
$
We know $scr(B)$ so we need to compute $P_"esc"$. We consider an isotropic distribution with a particle flux
$
  Phi (E) = dv(N, E)
$
hitting the _planar_ shock front. The flux emitted from one hemi-sphere of a source is
$
  F(E) &= integral dd(N, 2)/dd(E, Omega) cos theta dd(Omega) \
  &= Phi(E) integral_0^(2 pi) dd(phi.alt) integral_0^(pi\/2) cos theta sin theta dd(theta) \
  &= pi Phi(E)
$
Then
$
  n & = (4 pi)/c Phi(>= E) \
    & = 4/c F(>= E)
$
implying
$
  F(>= E) = (c n)/4
$
This is the average flux of particles through the shock front. Within $S_"shock"$ particles can escape the acceleration region with velocity $ v_2 = v_s/4 $ Then the escape flux is
$
  F_"esc" = n v_2 = (n v_s)/4
$
implying
$
  P_"esc" = F_"esc"/F = v_s/c
$
Then
$
  P = 1-v_s/c
$
and
$
  ln P tilde.eq -v_s/c";  " ln scr(B) tilde.eq v_s/c
$
implying
$
  N(E) tilde (E/E_0)^(-2)
$
which is reasonable!

== Hillas' criterion
The rate of energy gain is roughly
$
  dv(E, t) tilde.eq eta E 1/T_"cycle"
$
with $T_"cycle"$ being the time between scatterings. We approximate
$
  lambda_"cycle" tilde r_g = p_perp/(Z e B) = E/(Z e B c)
$
Then in $S_"shock"$
$
  T_"cycle" & = r_g/u \
            & = E/(Z e B u c)
$
and
$
  dv(E, t) & tilde.eq (Delta E)/T_"cycle" \
           & = eta Z e B u c \
           & tilde.eq Z e B u^2
$
The maximal is
$
  E_"max" & = dv(E, t) T_"shock" \
          & = Z e B u L
$
where $T_"shock" = L\/u$ with $L$ being the size of the accelerating region. We write this as
$
  (E_"max")/(Z beta_u) = e B c L
$
with $beta_u = u\/c$. This is often written in the form
$
  E_"max"/(Z beta_u) & tilde.eq 10^21 "eV" (B/"G") (L/"pc")
$
which is called _Hillas' criterion_. We see either $B$ must be large or $L$ must be large.
