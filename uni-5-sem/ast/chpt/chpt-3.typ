//**** init-ting
#import "@preview/physica:0.9.8": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()




#pagebreak()
= The equation of state
At any given point in a star we assume the gas is in thermodynamic equilibrium, this assumption allows us to just consider average properties of the gas. These are described by local state variables, e.g. temperature $T$ and density $rho$, and their relations. The relations between various state variables define the _equation of state_ for the gas. Due to the high temperatures of stars we can assume that our gas is fully ionized. Under this assumption we can fully ignore internal degrees of freedom. To a first approximation we can also neglect interactions between particles. This is what we call an ideal gas.

== The ideal gas
We just consider an ideal gas consisting of a single type of non-relativistic particles. In this case the ideal gas law is usually written as
$
  P V = N k_B T
$
in the case of stars we write
$
  P = n k_B T
$
with $n$ being the number density of particles. Defining the mass density $rho$ and atomic mass $mu$ of the particles we have
$
  n = rho/(mu m_u)
$
with $m_u$ being the atomic mass unit, we can write
$
  P = (rho k_B T)/(mu m_u)
$
this is the form we will use. By the equipartition theorem the internal energy per particle is $3k_B T\/2$ giving the total internal energy per volume
$
  u = 3/2 n k_B T = 3/2 (rho k_B T)/(mu m_u) = 3/2 P
$
note this assumes the Maxwell-Boltzmann distribution
$
  f(v) = 4 pi (m/(2 pi k_B T))^(3\/2) exp(- (m v^2)/(2 k_B T)) v^2
$
which holds for non-relativistic _classical_ particles.

=== Adiabatic processes
Consider the first law of thermodynamics
$
  dd(Q) = dd(U) + P dd(V)
$
let a system have unit mass, then $V = rho^(-1)$ and $U = u rho^(-1)$ so
$
  U = 3/2 (k_B T)/(mu m_u)
$
For a process with $dd(V)=0$ we have
$
  dd(Q) = 3/2 k_B/(mu m_u) dd(T) equiv c_V dd(T)
$
We can use the ideal gas law to write
$
  P dd(V) + V dd(P) = k_B/(mu m_u) dd(T)
$
so
$
  dd(Q) = dd(U) - V dd(P) + k_B/(mu m_u) dd(T) = 5/2 k_B/(mu m_u) dd(T) - V dd(P)
$
For a process with $dd(P)=0$ we have
$
  dd(Q) = 5/2 k_B/(mu m_u) dd(T) equiv c_P dd(T)
$

We care about adiabatic processes, which have $dd(Q)=0$. So
$
  c_V dd(T) = - P dd(V)
$
for an ideal gas
$
  c_V ( dd(P)/P + dd(V)/V) = - k_B/(mu m_u) dd(V)/V = (c_V-c_P) dd(V)/V
$
so
$
  dd(P)/P = - c_P/c_V dd(V)/V = -gamma dd(V)/V = gamma dd(rho)/rho
$
where we define the adiabatic index $gamma equiv c_P\/c_V$. We can also write this as
$
  (pdv(ln P, ln rho))_s = gamma
$
with $s$ implying $dd(Q)=0$. For an ideal gas we can also find
$
    (pdv(ln P, ln T))_s & = gamma/(gamma-1) \
  (pdv(ln T, ln rho))_s & = gamma -1
$

Under more general assumptions these are obviously not valid. But, they motivate defining adiabatic exponents $Gamma_i$ by the following
$
              Gamma_1 & equiv (pdv(ln P, ln rho))_s \
  Gamma_2/(Gamma_2-1) & equiv (pdv(ln P, ln T))_s \
          Gamma_3 - 1 & equiv (pdv(ln T, ln rho))_s
$
by the chain rule
$
  Gamma_2/(Gamma_2-1) = Gamma_1/(Gamma_3 -1)
$
and for an ideal gas we have $Gamma_i = gamma = 5\/3$.

We will come back to these later.

=== Radiation pressure
Due to the high temperature in a star photons contribute significantly to the pressure and energy of the gas. One can show that the pressure due to photons is
$
  P_R = 1/3 a T^4
$
and the energy density is
$
  u_R = a T^4
$
with $a = 4 sigma_"SB"\/c$ being the radiation energy constant, note $P_R = u_R\/3$ as expected for relativistic particles. If the particles behave like an ideal gas and radiation the total pressure and internal energy become
$
  P & = (rho k_B T)/(mu m_u) + 1/3 a T^4 \
  u & = 3/2 (rho k_B T)/(mu m_u) + a T^4
$

== Degenerate matter
At low temperatures and high densities we need to incorporate quantum mechanical effects to describe our gas. This is because of Pauli's exclusion principle.

To do this we work with the Fermi-Dirac distribution
$
  f(p) = 1/(exp[(E\/k_B T - psi)+1])
$
which gives the number density
$
  n(p) dd(p) = 2/h^3 4 pi p^2 dd(p) f(p)
$
note that $psi$ is related to the chemical potential and the factor $h^(-3)$ gives the density of states.

In the _complete degeneracy_ limit with $T->0$ where $f(p) -> Theta$ with $Theta$ being the Heaviside step-function. One can compute
$
  n = (8 pi p_F^3)/(3 h^3)
$
in the non-relativistic case using $E = p^2\/(2 m)$ one can also compute
$
  u = 3/5 n E_F
$
note in this case $E_F = p_F^2\/(2m)$ so
$
  P = 2/5 n E_F
$
in the relativistic case we instead use $E = p c$ giving
$
  u & = 3/4 n E_F \
  P & = 1/4 n E_F
$

#pagebreak()
= Hydrostatic equilibrium
Since most stars do not change on small time scales this indicates that the forces acting on stars are in perfect balance.

To see what this means consider a mass shell of some spherically symmetric star of area $dd(A)$ extending from $r$ to $r + dd(r)$. We want an equation of motion for this shell. The force due to gravity will be due to all the mass contained within $r$, we denote this by $m equiv m(r)$. Then the gravitational acceleration is $-G m\/r^2$. The volume of the shell is $dd(r, A)$, so if the density of the shell is $rho$ then is has mass $rho dd(r, A)$. Therefore
$
  dd(F_g) = - rho (G m)/r^2 dd(r, A)
$
The other force acting on the shell is due to pressure, and by definition we have
$
  dd(F_P) tilde.eq P(r) dd(A) - P(r+dd(r)) dd(A) tilde.eq^"Taylor" - dv(P, r) dd(r, A)
$
So the equation of motion reads
$
  rho dd(r, A) dv(r, t, 2) = dd(F_g)+dd(F_P)
$
implying
$
  rho dv(r, t, 2) = - rho (G m)/r^2 - dv(P, r)
$
Assuming the star is in equilibrium then $dot.double(r) =^! 0$ meaning
$
  dv(P, r) = - (G m rho)/r^2
$
this is the equation of hydrostatic equilibrium---and is the first equation of stellar structure.

By the definition of $m$ we immediately find
$
  dd(m) = 4 pi r^2 rho dd(r) => dv(m, r) = 4 pi r^2 rho
$
this is the equation of mass continuity---and is the second equation of stellar structure.

These equations can not be solved as they are without further information. But, if we know $rho = rho(r)$ or $rho = rho(P)$ then they can be solved. These simple models are not necessarily good models.

== The linear model
We assume $rho$ is linear in $r$, with $rho (r=0) = rho_c$ and $rho (r=R) = 0$. This can be written as
$
  rho = rho_c (1- r/R)
$
then by using mass continuity we can find
$
  m = (4 pi rho_c r^3)/3 (1- (3 r)/(4 R)) = M (4 x^3 - 3 x^4 )
$
where we let the total mass of the star be $m(r=R) = M$ and define $x equiv r\/R$. It follows that
$
  rho_c = (3 M)/(pi R^3)
$
Then by hydrostatic equilibrium one can obtain
$
  P = 5/(4 pi) (G M^2)/R^4 (1 - 24/5 x^2 + 28/5 x^3 - 9/5 x^4)
$
where we assume $P(r=R)=0$. Assuming the gas is ideal it is also possible to determine the temperature
$
  T = 5/12 (G mu m_u M)/(k_B R) (1 + x- 19/5 x^2 + 9/5 x^3)
$

== The isothermal atmosphere
We assume that the temperature $T$ is constant in the atmosphere. We also assume that the atmosphere is small compared to the radius of the star. Then $g = G M\/R^2$ can be taken as constant. Assuming the gas is ideal we find
$
  dv(P, r) = - g rho = - P/H
$
with the pressure scale height $H$ being defined as
$
  H equiv (k_B T)/(g mu m_u)
$
We can integrate the above to find
$
  P = P_0 e^(- h\/H)
$
where $h = r-r_0$ with $r_0$ being some arbitrary reference level in the atmosphere, and $P_0 = P(h=0)$. The by the ideal gas law
$
  rho = rho_0 e^(-h\/H)
$
with $rho_0 = rho(h=0)$.

So both density and pressure decrease exponentially in an isothermal atmosphere.

== Polytropic models
We assume a relation of the form
$
  P = K rho^gamma
$
with $K$ and $gamma$ being constants.

Consider
$
  dv(, r) (r^2/rho dv(P, r)) = - G dv(m, r) = - 4 pi G rho r^2
$
using $P = K rho^gamma$
$
  dv(, r) (r^2 K gamma rho^(gamma-2) dv(rho, r)) = - 4 pi G rho r^2
$
Now we define the polytropic index $n$ by
$
  n equiv 1/(gamma-1)
$
and introduce the dimensionless $theta$ by
$
  rho = rho_c theta^n
$
then the equation becomes
$
  ((n+1)K rho_c^(1\/n-1))/(4 pi G) 1/r^2 dv(, r) (r^2 dv(theta, r)) = - theta^n
$
Introducing $r = alpha xi$ with
$
  alpha^2 equiv ((n+1) K rho_c^(1\/n-1))/(4 pi G)
$
we can write
$
  1/xi^2 dv(, xi) (xi^2 dv(theta, xi)) = - theta^n
$
which is known as the Lane-Emden equation, with $theta = theta_n (xi)$ being the Lane-Emden function. This must satisfy the boundary condition $theta_n = 1$ for $xi = 0$, and the surface is defined by $xi = xi_1$ where $theta = 0$.

Given $theta_n$ we can find
$
  R = [((n+1) K rho_c^(1\/n-1))/(4 pi G)]^(1\/2) xi_1
$
and from the mass continuity equation one can find
$
  m(xi) &= integral_0^(a xi) 4 pi r^2 rho dd(r) \
  &= 4 pi alpha^3 rho_c integral_0^xi xi^2 theta_n^n dd(xi) \
  &= - 4 pi alpha^3 rho_c integral_0^xi dv(, xi) (xi^2 dv(theta_n, xi)) dd(xi) \
  &= - 4 pi alpha^3 rho_c xi^2 dv(theta_n, xi)
$
with the total mass being $M = m(xi=xi_1)$.

By assumption it follows that
$
  P_c = K rho_c^((n+1)\/n)
$
and
$
  P = P_c theta_n^(n+1)
$
similarly to above assuming the gas is ideal one can find an expression for the temperature.

In general the Lane-Emden equation must be solved numerically, but we will not cover this here.

#pagebreak()
= Radiative energy transport
Above we had $u_R = a T^4$ and since $T$ decreases with increasing $r$ photons moving away from the centre of a star carry more energy than photons moving toward the centre. This leads to a net transport of energy toward the surface.

We consider the energy transport in a time interval $dd(t)$ through an area $dd(A)$, orthogonal to the direction to the centre of the star, at $r = r_0$. We describe the motion of photons by $theta$ which we define to be the angle between the outward normal of $dd(A)$ and the direction of motion. Assuming this motion is isotropic the fraction of photons with directions between $theta$ and $theta + dd(theta)$ is
$
  ("band with" dd(A))/"unit sphere"=(2 pi sin theta dd(theta))/(4 pi) = 1/2 sin theta dd(theta)
$
The motion of photons through the gas is determined by their mean free path $lambda$. In the interior of a star it is fair to assume $lambda << r$. The average photon going through $dd(A)$ with a direction between $theta$ and $dd(theta)$ will come from the distance $r' = r_0 - lambda cos theta$, so they correspond to the energy density $u_R (r')$. Their contribution to the energy transport through $dd(A)$ is therefore
$
  underbrace(1/2 sin theta dd(theta) u_R (r_0 - lambda cos theta), "contribution") underbrace(cos theta dd(A), "area") underbrace(c dd(t), "path length")
$
we integrate over all directions to find the energy
$
  dd(E) &= 1/2 c integral_0^pi u_R (r_0-lambda cos theta) cos theta sin theta dd(theta, A, t) \
  &tilde.eq^"Taylor" 1/2 c integral_0^pi u_R (r_0) cos theta sin theta dd(theta, A, t) - 1/2 lambda c integral_0^pi dv(u_R, r) cos^2 theta sin theta dd(theta, A, t) \
  &= - (lambda c)/3 dv(u_R, r)dd(A, t)
$
we write this as
$
  dd(E) = F_R dd(A, t)";  " F_R equiv - (lambda c)/3 dv(u_R, r)
$
Typically one uses the opacity $kappa$ defined by
$
  lambda = 1/(kappa rho)
$
with $u_R$ written in terms of $T$ to obtain
$
  F_R = - (4 a c T^3)/(3 kappa rho) dv(T, r)
$
Assuming that energy transport only occurs through radiation we can write
$
  L(r) = 4 pi r^2 F_R
$
giving
$
  dv(T, r) = - (3 kappa rho L(r))/(16 pi a c r^2 T^3)
$
which is the equation describing radiative energy transport---an alternative derivation uses $P_R$.

At the surface of a star the density is very low meaning the free path is large. Therefore the above is invalid. However, we can approximate the energy radiated by only taking the contribution from outgoing photons, those with $theta < pi\/2$. This gives
$
  dd(E) tilde.eq 1/2 c integral_0^(pi\/2) u_R cos theta sin theta dd(theta, A, t) = 1/4 c u_R dd(A, t)
$
which is the familiar Stefan-Boltzmann law
$
  F_R tilde.eq (a c)/4 T^4 = sigma T^4
$
Using this we can define an effective temperature
$
  F_R = sigma T_"eff"^4
$
Then the surface luminosity $L_s = L(r)$ is
$
  L_s = 4 pi R^2 sigma T_"eff"^4
$
with $R$ being the surface radius.

== The energy equation
The equation of radiative energy transport must be supplemented by an equation for the luminosity $L$ in terms of $r$. By definition we can write
$
  dv(L, r) = 4 pi r^2 rho epsilon.alt
$
where $epsilon.alt$ is the rate of energy production per unit mass.

This does not include gravitational energy, and is therefore not valid during early stellar evolution. Instead consider the first law of thermodynamics
$
  dd(Q) = dd(U) + P dd(V)
$
for a system containing unit mass, meaning $V = rho^(-1)$ and $U = u rho^(-1)$, the heat $dd(Q)$ has two contributions. The first is given by $epsilon.alt$, the heat liberated by fusion. The second is the heat deposited from or extracted by the energy flowing through the layer. Combining these we obtain
$
  underbrace(epsilon.alt dd(t), "fusion") + underbrace((L(r) - L(r+dd(r)))/(4 pi r^2 rho dd(r)) dd(t), "energy flow") = dd((u/rho)) + P dd((1/rho))
$
or rewriting
$
  dv(L, r) = 4 pi r^2 [rho epsilon.alt - rho dv(, t) (u/rho) + P/rho dv(rho, t)]
$
this is the energy equation. The second and third term are typically negligible during much of a stars life.

#pagebreak()
= Energy transport by convection
Assuming all energy transport happens through radiation is a very bad assumption, since such stars are generally unstable.

A typical instability is caused by having a layer of high density on top of a layer of low density. This type of instability happens if the temperature decreases too rapidly with distance from the centre. The pressure decrease is determined by hydrostatic equilibrium, and is essentially given. The only way to compensate for the decrease in temperature is, by the ideal gas law, that the density decreases slowly or increases. Due to the instability hotter, light elements of fluid rise and cooler, heavier elements sink. When this motion becomes sufficiently strong, the elements are dissolved and the gas is mixed. As a result rising elements deposit excess heat, leading to a net transport of energy. This is called convection, and the instability is usually called convective instability.

== Instability condition
We want to determine the condition for instability. Consider an element of gas which is moved a distance $Delta r$ outwards. Denote the pressure and density outside (inside) the element by $P_i, rho_i$ $(P_i^*, rho_i^*)$, with $i = 1$ before moving and $i = 2$ after moving.

Initially the element is identical to its surroundings $P_1^* = P_1$ and $rho_1^* = rho_1$. The motion is determined by buoyancy
$
  f_"buoy" = - g (rho_2^* - rho_2) equiv -g Delta rho
$
if $f_"buoy" > 0$ the motion is accelerated, hence it is unstable. So we want to determine $Delta rho$. To do this we make to assumptions: the motion is so slow that there is pressure balance between the element and its surroundings, the motion is so fast that there is no heat loss to the surroundings. By the first we have $P_2^* = P_2$ the second says that the motion is adiabatic meaning
$
  dd(rho^*)/rho^* = 1/Gamma_1 dd(P^*)/P^* = 1/Gamma_1 dd(P)/P
$
