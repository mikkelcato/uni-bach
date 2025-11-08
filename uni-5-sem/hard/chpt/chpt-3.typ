//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures
#import "@preview/mannot:0.3.0": *
#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
#text(size: 35pt, strong("Electrons"))
= Free electron models
Now we'll begin to discuss electrons and their properties---important for thermal/electrical conductivity, heat capacity etc. We proceed in the following way:
$
  "FEM" ->^"folding" "empty lattice" ->^"weak potential" "NFE" ->^"strong potential" "TB"
$
we'll see that the addition of a periodic crystal potential leads to interesting band structure, such as band gaps.

The idea of free electron models is that valence electrons of the constituent atoms of metals become conduction electrons and move freely through the volume of the metal. This started as a classical model with great success, but changed with the introduction of quantum mechanics which explained many failures with the classical model. When we talk about a free electron Fermi gas, we mean a gas of free electrons subject to the exclusion principle.
== Sommerfeld theory
=== In one dimension
Letting $V = 0$ in some box of length $L$ we have
$
  H psi_n = - hbar^2/(2 m) dv(psi_n, x, 2) = epsilon.alt_n psi_n
$
by orbital we mean a solution of the wave equation for a system with only one electron. In this way we distinguish between a system of $N$ interacting electrons, and some approximate system of $N$ eletroncs in $N$ different orbitals---this is exact assuming no interactions. Applying the usual $psi_n (0) = 0$ and $psi_n (L) = 0$. Then
$
  psi_n = A sin (2 pi x)/lambda_n " with " (n lambda_n)/2 = L
$
with energy
$
  epsilon.alt_n = hbar^2/(2 m) ((n pi)/L)^2
$
by the exclusion principle each orbital can only be occupied by at most one electron. In a linear solid the quantum numbers of interest for a conduction electorn are $n$ and $m_s$, with $n$ being any positive integer, and $m_s = plus.minus 1/2$---the spin. So a pair of orbitals labeled by $n$ can be occupied by at most two electrons. If multiple orbitals have the same energy we say they are degenerate.

Define $n_F$ to be the topmost filled energy level, if we start filling from $n = 1$ and continue filling until all $N$ electrons occupy some orbital. Let $N$ be even for convenience then $2 n_F = N$. We define the Fermi energy $epsilon.alt_F$ to be the energy of $n_F$ (the ground state of the $N$ electron system),
$
  epsilon.alt_F = hbar^2/(2 m) ((n_F pi)/L)^2 = hbar^2/(2m) ((N pi)/(2 L))^2
$

=== With temperature
The ground state is the state which the $N$ electron system takes at absolute zero. If we increase $T$ then we introduce kinetic energy and some previously empty energy levels can now be occupied by excited electrons. The Fermi-Dirac distribution gives the probability that an orbital with energy $epsilon.alt$ is occupied
$
  f(epsilon.alt) = 1/(exp[(epsilon.alt-mu)\/k_B T]+1)
$
this distribution is valid for all fermions---note the Bose-Einstein distribution valid for all bosons would have a minus sign.

The chemical potential $mu$ is defined such that the total number of electrons to come out to $N$. At absoulute zero $mu = epsilon.alt_F$, since in the limit $T -> 0$ the distribution changes discontinuously from $1 -> 0$ (filled $->$ empty) at $epsilon.alt = epsilon.alt_F = mu$. In the limit $epsilon.alt - mu >> k_B T$ we rederive the Boltzmann distribution.

=== In three dimensions
The above easily generalizes to three dimensions giving
$
  psi_bold(k) (bold(r)) = exp(i bold(k) dot bold(r))
$
given that all components of $bold(k)$ satisfy
$
  k_i = 0, plus.minus (2 pi)/L, dots, (2 n pi)/L, dots
$
these ensure periodicity for a cube of length $L$, where $n$ is any integer. The components of $bold(k)$ act as our quantum numbers, alongside $m_s$. The energies are then
$
  epsilon.alt_bold(k) = hbar^2/(2 m) k^2 = hbar^2/(2m) (k_x^2 +k_y^2+ k_z^2 )
$
note that
$
  bold(p) psi_bold(k) (bold(r)) = - i hbar nabla psi_bold(k) (bold(r)) = hbar bold(k) psi_bold(k) (bold(r))
$
so $psi_bold(k)$ is an eigenfunction of $bold(p)$ with eigenvalue $hbar bold(k)$. In the ground state of a system of $N$ electrons, the occupied orbitals correspond to points inside a sphere in $bold(k)$ space. The energy at the surface of this sphere is the Fermi energy, the wavevectors at the Fermi surface have magnitude $k_F$ defined by
$
  epsilon.alt_F = hbar^2/(2 m) k_F^2
$
by the restrictions on the triplet $k_i$ we see there is one distinct triplet for each volume element $(2 pi\/L)^3$ of $bold(k)$ space. So in the sphere with volume $4 pi k_F^3\/3$ the total number of orbitals $N$ is
$
  2 (4 pi k_F^3\/3)/(2pi\/L)^3 = V/(3 pi^2) k_F^3 = N
$
so
$
  k_F = ((3 pi^2 N)/V)^(1\/3) => epsilon.alt_F = hbar^2/(2 m) ((3 pi^2 N)/V)^(2\/3)
$
depending only on $N\/V$, we can further define the electron velocity on the Fermi surface
$
  v_F = (hbar k_F)/m = hbar/m ((3 pi^2 N)/V)^(1\/3)
$
we could also define the Fermi temperature $T_F = epsilon.alt_F\/k_B$, which can be treated like the relevant temperature scale. From the expression for $epsilon.alt_F$ we can also find the density of states (strictly the density of orbitals)
$
  D(epsilon.alt) equiv dv(N, epsilon.alt) = V/(2 pi^2) ((2 m)/hbar^2)^(3\/2) epsilon.alt^(1\/2) = (3 N)/(2 epsilon.alt)
$

== Heat capacity
Classically the predicted heat capacity should be $3/2 N k_B$ by the equipartition theorem, since we have $N$ free particles---the observed electronic contribution is much, much smaller. This problem is solved by the Fermi-Dirac distribution, when we heat a specimen not every electron gains energy $tilde k_B T$, but only those within $k_B T$ of the Fermi level can be excited thermally.

We derive an expression valid for $k_B T << epsilon.alt_F$. The increase in energy $dd(U, d: Delta) equiv U(T)-U(0)$ for a system of $N$ electrons heated from $0 -> T$ is
$
  dd(U, d: Delta) = integral_0^oo dd(epsilon.alt) epsilon.alt D(epsilon.alt) f(epsilon.alt) - integral_0^epsilon.alt_F dd(epsilon.alt) epsilon.alt D(epsilon.alt)
$
with $f(epsilon.alt)$ being the Fermi-Dirac distribution. Consider
$
  N = integral_0^oo dd(epsilon.alt) D(epsilon.alt) f(epsilon.alt) = integral_0^epsilon.alt_F dd(epsilon.alt) D(epsilon.alt)
$
and multiply by $epsilon.alt_F$ then
$
  (integral_0^epsilon.alt_F + integral_(epsilon.alt_F)^oo) dd(epsilon.alt) epsilon.alt_F f(epsilon.alt) D(epsilon.alt) = integral_0^epsilon.alt_F dd(epsilon.alt) epsilon.alt_F D(epsilon.alt)
$
this gives
$
  dd(U, d: Delta) = integral_(epsilon.alt_F)^oo dd(epsilon.alt) (epsilon.alt-epsilon.alt_F) f(epsilon.alt) D(epsilon.alt) + integral_0^epsilon.alt_F dd(epsilon.alt) (epsilon.alt_F-epsilon.alt) [1-f(epsilon.alt)] D(epsilon.alt)
$
this gives the compact result
$
  C_"el" = dv(U, T) = integral_0^oo dd(epsilon.alt) (epsilon.alt-epsilon.alt_F) dv(f, T) D(epsilon.alt)
$
is should be clear that $(epsilon.alt-epsilon.alt_F) dd(f)\/dd(T)$ peaks sharply around $epsilon.alt_F$ for this reason we approximate
$
  C_"el" tilde.equiv D(epsilon.alt_F) integral_0^oo dd(epsilon.alt) (epsilon.alt - epsilon.alt_F) dv(f, T)
$
also when $k_B T << epsilon.alt_F$ we can ignore the $T$ dependence of $mu$, so we can replace $mu -> epsilon.alt_F$. Then with $x equiv (epsilon.alt-epsilon.alt_F)\/tau$ with $tau equiv k_B T$ we can write
$
  C_"el" = k_B^2 T D(epsilon.alt_F) integral_(-epsilon.alt_F\/tau)^oo dd(x) x^2 e^x/(e^x+1)^2
$
extending the lower limit to $- oo$, since the integrand is essentially zero for lower values we find
$
  C_"el" = 1/3 pi^2 D(epsilon.alt_F) k_B^2 T =^"f.e.g" (pi^2 N k_B)/(2) ( T/T_F )
$

Now for temperatures below both the Debye temperature $theta$ and the Fermi temperature $T_F$ then we can write the heat capacity as
$
  C = C_"lat" + C_"el" tilde.equiv gamma T + A T^3
$
it is convenient in experiment to use
$
  C/T tilde.equiv gamma + A T^2
$
and then treat the variable as $T^2$.

== Drude model and conductivity
By the Lorentz force we have
$
  bold(F) = m dv(bold(v), t) = hbar dv(bold(k), t) = - e (bold(E) + 1/c bold(v) times bold(B))
$
if there are no collisions then the Fermi sphere moves in $bold(k)$ space at a uniform rate due to some applied constant electric field. Taking $bold(B)=0$ we can integrate this to obtain
$
  bold(k) (t) - bold(k) (0) = - e bold(E) t\/hbar
$
if $bold(F) = - e bold(E)$ is applied at $t = 0$ then at some later time the center of our Fermi sphere will be displace by $ dd(bold(k), d: delta) = - e bold(E) t\/hbar $
If the collision time is given by $tau$ (we assume stopping is instanteneous as is typical in the Drude model)---typically dependent on impurities and phonons---then the displacement of the Fermi sphere in steady state is given by this with $t = tau$. The velocity is $bold(v) = dd(bold(k), d: delta)\/m = - e bold(E) tau\/m$. If in a constant $bold(E)$ we have $n$ electrons with charge $q = - e$ then the current density is
$
  bold(j) = n q bold(v) = n e^2 tau bold(E)\/m = sigma bold(E)
$
which is Ohm's law, and we've defined the electrical conductivity
$
  sigma = (n e^2 tau)/m
$
the resistivity is just the inverse $rho = sigma^(-1)$, experimentally we have $rho = rho_L + rho_i$ since this guy depends on thermal phonons $rho_L$ (through e.g. Umklapp-scattering) and impurities $rho_i$. Usually the first is independent of the number of defects, while the second is independent of temperature---Matthiessen's rule.

=== Thermal conductivity
We previously found $K = 1/3 C v cal(l)$. For a Fermi gas we use $epsilon.alt_F = 1/2 m v_F^2$, and using $cal(l) = v_F tau$ we find
$
  K_"el" = pi^2/3 (n k_B^2 T)/(m v_F^2) v_F cal(l) = (pi^2 n k_B^2 T tau)/(3 m)
$

==== Wiedemann-Franz law
We find
$
  K/sigma = (pi^2 k_B^2 T n tau\/3 m)/(n e tau^2 \/m) = pi^2/3 (k_B/e)^2 T
$
where we use the $sigma$ and $K$ previously found. We define the Lorenz number
$
  L = K/(sigma T)
$
notably this doesn't depend on $n$ or $m$. So the ratio of thermal conductivity to electrical conductivity is $prop T$.

== Magnetic fields
We can write the equation of motion for $dd(bold(k), d: delta)$ acted on by a force $bold(F)$ and a friction represented by $1\/tau$
$
  hbar (dv(, t) + 1/tau) dd(bold(k), d: delta) = bold(F)
$
if $bold(F)$ is given by the Lorentz force with $m bold(v) = hbar dd(bold(k), d: delta)$ then
$
  m (dv(, t)+1/tau) bold(v) = -e (bold(E) + 1/c bold(v) times bold(B))
$
if we let $bold(B) = B hat(bold(z))$ and assume steady-state with a static $bold(E)$ then the drift velocities are
$
  v_x & = - (e tau)/m E_x - omega_c tau v_y \
  v_y & = - (e tau)/m E_y + omega_c tau v_x \
  v_z & = - (e tau)/m E_z
$
where we've defined the cyclotron frequency
$
  omega_c = (e B)/(m c)
$

/*
=== Hall effect
The Hall field is the electric field developed across two faces of a conductor in the direction $bold(j) times bold(B)$ when a current $bold(j)$ flows accros a magnetic field $bold(B)$. If we consider some rod-shaped object with a longitudinal electric field $E_x$ and a transverse magnetic field. We assume $dd(v_y, d: delta) = 0$ for this we require
$
  E_y = - omega_c tau E_x = - (e B tau)/(m c) E_x
$
so there is a transverse electric field. We define the Hall coefficient
$
  R_H = E_y/(j_x B)
$
then using Ohm's law $j_x = n e^2 tau E_x \/m$ then
$
  R_H = - 1/(n e c)
$
with CGS.
*/

#pagebreak()
= Nearly free $e^-$ model
The previous free electron model is very naive, meaning we miss a lot of detail including band structure---so we'd like a more detailed model incorporating the crystal lattice.

Before we found
$
  epsilon.alt_bold(k) = hbar^2/(2 m) (k_x^2 + k_y^2 + k_z^2)
$
where for periodicity over a cube of length $L$ we have
$
  k_i = 0; plus.minus (2pi)/L dots
$
with the free electron wavefunctions having the form
$
  psi_bold(k) (bold(r)) = exp(i bold(k) dot bold(r))
$
these represent free waves with $bold(p) = hbar bold(k)$.

The nearly free electron model builds upon this by supposing electrons are weakly perturbed by the periodic potential from the ion cores---we'll see that this causes band gaps.
/*
Bragg reflection of electron waves then occurs, and causes the band gaps. As an example consider a linear lattice with lattice constant $a$, then the Bragg condition becomes
$
  k = plus.minus 1/2 G = plus.minus (n pi)/a
$
the first reflections (and energy gap) occur at $k = plus.minus pi\/a$. The wavefunctions at these $k$ are not travelling waves $exp(plus.minus i pi x \/a)$ of free electrons. Instead they consist of equal parts left and right travelling wave---it is a standing wave. These are
$
  psi(plus.minus) = exp((i pi x)/a) plus.minus exp((-i pi x)/a) = cases(2 cos pi x\/a, 2 i sin pi x\/a)
$
but why does this lead to an energy gap? The two wavefunctions gather electrons at different regions, leading to them having different potential energies---this difference is the energy gap. We can see this by
$
  rho(plus.minus) = abs(psi(plus.minus))^2 prop cases(cos^2 pi x\/a, sin^2 pi x\/a)
$
so $psi(+)$ gathers electrons preferably at the ion cores where the potential is smallest, oppositely for $psi(-)$. So we find that energetically $rho(+) < "free wave" < rho(-)$, we let the energy gap be $rho(-)-rho(+) = E_g$. Below $E_g$ the wavefunction is $psi(+)$ and above $E_g$ it is $psi(-)$.
*/

== Bloch's theorem
Bloch proved that solutions to the Schrödinger equation for a periodic potential takes the form
$
  psi_bold(k) (bold(r)) = u_bold(k) (bold(r)) exp(i bold(k) dot bold(r))
$
with $u_bold(k) (bold(r)) = u_bold(k) (bold(r)+bold(T))$. So any solution $psi_bold(k)$ takes the form of a free wave with some periodic modulation $u_bold(k)$---these are sometimes called Bloch functions. A proof of Bloch's theorem is given below.

/*
=== An example---Kronig-Penney
We consider a simple periodic potential, the square-well array. We know
$
  - hbar^2/(2 m) dv(psi, x, 2) +U(x) psi = epsilon.alt psi
$
in $0 < x < a$ where $U = 0$ then
$
  psi = A exp(i K x) + B exp(-i K x)
$
with
$
  epsilon.alt = (hbar^2 K^2)/(2 m)
$
in $-b < x < 0$ (within the barrier) we have
$
  psi = C exp(Q x) + D exp(-Q x)
$
where
$
  U_0 - epsilon.alt = (hbar^2 Q^2)/(2 m)
$
this must be related to the solution in $a < x < a+b$ by Bloch's theorem
$
  psi(a < x < a +b) = psi(-b < x < 0) exp(i k (a+b)))
$
we pick ${A,B,C,D}$ such that $psi$ and its derivative are continuous at $x=0$ and $x = a$, at the first
$
     A + B & = C + D \
  i K(A-B) & = Q(C-D)
$
and by the second (also using Bloch's theorem)
$
  A exp(i K a) + B exp(-i K a) &= (C exp(-Q b) + D exp(Q b)) exp(i k(a+b)) \
  i K ( A exp(i K a) - B exp(-i K a)) &= Q (C exp(-Q b) - D exp(Q b)) exp(i k (a+b))
$
this gives
$
  [(Q^2-K^2)/(2 Q K)] sinh Q b sin K a + cosh Q b cos K a = cos k(a+b)
$
this is simpler if we write the potential as a periodic $delta$-function obtained by the limit $b -> 0$ and $U_0 -> oo$ keeping $P = Q^2 b a\/2$ finite. In this case $Q >> K$ and $Q b << 1$ giving
$
  P/(K a) sin K a + cos K a = cos k a
$
this only has solutions for certain $K a$ and these determine the possible energies and corresponding gaps---they must be between $plus.minus 1$ due to the $cos k a$ on the RHS.
*/

=== Central equation
We now treat the wave equation for a general potential at general $bold(k)$. Let $U(bold(r))$ be the potential energy of an electron in the lattice, then we know $U(bold(r)+bold(T))=U(bold(r))$. So we can write the potential as a Fourier series
$
  U(bold(r)) = sum_bold(G) U_bold(G) e^(i bold(G) dot bold(r))
$
the $U_bold(G)$ tend to decrease rapidly.

We now assume $U(bold(r))$ is symmetric about $bold(r)=bold(0)$ and $U_bold(0) = bold(0)$. The Schrödinger equation then becomes
$
  (bold(p)^2/(2 m) + sum_bold(G) U_bold(G) e^(i bold(G)dot bold(r))) psi(bold(r)) = epsilon.alt psi(bold(r))
$
this describes an electron in the potential of ion cores, and in the average potential of other electrons. Note that $U(bold(r)) in RR$ since $U_(-bold(G))=U_bold(G)^*$ and any Bravais lattice has inversion symmetry $U_(-bold(G))=U_bold(G)$.

We can also expand $psi(bold(r))$ as
$
  psi(bold(r)) = sum_bold(k) C(bold(k)) e^(i bold(k)dot bold(r))
$
where $bold(k)$ ensures this satisfies boundary conditions---$k_i = 2 pi n_i\/L$. Not all $bold(k)$ of the given form enter the expansion, since if any $bold(k)$ is in the expansion, then all $bold(k) + bold(G)$ enter aswell. So we label the Bloch functions by the $bold(k)$ in the first Brillouin zone. Note that $hbar bold(k)$ is called the crystal momentum, since $bold(k)$ represents a phase factor picked up by $psi(bold(r))$ after a translation $bold(r)->bold(r)+bold(T)$, and since this is what enters in conservation laws.

Substituting both expansions into the Schrödinger equation gives
$
  sum_bold(k) (hbar^2 k^2)/(2 m) C(bold(k)) e^(i bold(k) dot bold(r)) + sum_(bold(G),bold(k)) U_bold(G) C (bold(k)) e^(i(bold(k)+bold(G))dot bold(r)) = epsilon.alt sum_bold(k) C(bold(k)) e^(i bold(k)dot bold(r) )
$
by matching exponential factors one can see this implies:
$
  markrect((epsilon.alt_k^0 - epsilon.alt) C(bold(k)) + sum_bold(G) U_bold(G) C (bold(k)-bold(G)) = 0"  with  " epsilon.alt_k^0 = (hbar^2 k^2)/(2 m), outset: #.4em)
$
which is called the central equation!---notice $epsilon.alt_k^0$ is the free electron energy, as can also be seen in the limit $U_bold(G) = 0$. There is one central equation per Fourier component of $psi_bold(k)$.

#proof[ of Bloch's theorem][

  We can write
  $
    psi_bold(k) (bold(r)) & = sum_bold(G) C (bold(k)-bold(G)) e^(i(bold(k)-bold(G)) dot bold(r)) \
    & = e^(i bold(k) dot bold(r)) (sum_bold(G) C (bold(k)-bold(G)) e^(-i bold(G) dot bold(r))) \
    &= e^(i bold(k) dot bold(r)) u_bold(k) (bold(r))
  $
  where
  $
    u_bold(k) (bold(r)) equiv sum_bold(G) C (bold(k)-bold(G)) e^(-i bold(G) dot bold(r))
  $
  this satisfies $u_bold(k) (bold(r) + bold(T)) = u_bold(k) (bold(r))$
  $
    u_bold(k) (bold(r) + bold(T)) & = sum_bold(G) C (bold(k)-bold(G)) e^(-i bold(G) dot (bold(r)+bold(T))) \
    & =e^(-i bold(G)dot bold(T)) [sum_bold(G) C (bold(k)-bold(G)) e^(-i bold(G) dot bold(r))] \
    & = underbrace(e^(-i bold(G) dot bold(T)), equiv 1) u_bold(k) (bold(r))
  $
]

== Kronig-Penney potential
We consider a (one-dimensional) potential of the form
$
  U(x) = 2 sum_(G > 0) U_G cos G x = U_0 a sum_s delta(x - s a)
$
for a linear lattice with spacing $a$ and length $L$ we find
$
  U_G & = 1/L integral_0^L dd(x) U(x) cos G x \
      & = U_0/L a sum_s integral_0^1 dd(x) delta(x-s a) cos G x \
      & = U_0/L a sum_s underbrace(cos G s a, cos(2 pi n s)) \
      & = U_0
$
then the central equation becomes
$
  (epsilon.alt_k^0 - epsilon.alt) C_k + U_0 sum_n C(k - (2 pi n)/a) = 0
$
define
$
  f(k) equiv sum_n C(k - (2 pi n)/a)
$
then
$
  C_k = - ((2 m U_0\/hbar^2) f(k))/(k^2 - (2 m epsilon.alt\/hbar^2))
$
using
$
  f(k) = f(k - (2 pi n)/a)
$
we can then write
$
  C(k - (2 pi n)/a) = - ((2 m U_0)/hbar^2) f(k) [(k-(2 pi n)/a)^2 - (2 m epsilon.alt)/hbar^2]^(-1)
$
and summing over all $n$ we obtain
$
  hbar^2/(2 m U_0) = - sum_n [(k-(2pi n)/a)^2 - ((2 m epsilon.alt)/hbar^2)]^(-1)
$
using
$
  "ctn" x = sum_n 1/(n pi + x)
$
defining $(k_epsilon.alt^0)^2 equiv 2 m epsilon.alt\/hbar^2$ we can write the result as
$
  -underbracket((m U_0 a^2)/(2 hbar^2), equiv P) 1/(k_epsilon.alt^0 a) sin k_epsilon.alt^0 a + cos k_epsilon.alt^0 a = cos k a
$
this shows band gaps since $abs(cos k a)<= 1$ so only certain $K tilde epsilon.alt$ are allowed---see Figure ...

== For nearly free $e^-$
When applying the central equation to nearly free electrons we assume that the $U_bold(G)$ are small.

We start by considering the simplest case, namely that $U_bold(G) eq.not bold(0)$ when $bold(G)$ is the shortest reciprocal lattice vector---so it only couples $bold(k)$ differing by $plus.minus bold(G)$. We get two equations:
$
  bold(k)\:"  "& (epsilon.alt_bold(k)^0-epsilon.alt) C(bold(k)) + U_(bold(G)) C(bold(k)-bold(G)) = 0 \
  bold(k)-bold(G)\:"  "& (epsilon.alt_(bold(k)-bold(G))^0-epsilon.alt) C(bold(k)-bold(G))+U_(bold(G)) C(bold(k)) = 0
$

this system of equations is easily solved to give
$
  epsilon.alt_plus.minus = (epsilon.alt_bold(k)^0 + epsilon.alt_(bold(k)-bold(G))^0)/2 plus.minus [(epsilon.alt_bold(k)^0-epsilon.alt_(bold(k)-bold(G))^0)^2/4 + U_bold(G)^2]^(1\/2)
$
we get two solutions!---meaning we get a band gap. In particular if $bold(k)$ lies on the edge of the first Brillouin zone with $bold(k) = 1/2 bold(G)$ then $epsilon.alt_bold(k)^0 = epsilon.alt_(bold(k)-bold(G))^0$ giving
$
  markrect(epsilon.alt_plus.minus = epsilon.alt_bold(k) plus.minus abs(U_bold(G)) => E_g = 2U_bold(G), outset: #.3em)
$
Consider now a $bold(k)$ close to the edge of the first Brillouin zone, with $bold(K) = bold(k) - 1/2 bold(G) tilde.eq 0$. This eventually gives
$
  markrect(epsilon.alt_plus.minus &tilde.eq^("Taylor") epsilon.alt_plus.minus^0 + (hbar^2 bold(K)^2)/(2 m) (1 plus.minus (2 epsilon.alt_(1/2 bold(G))^0)/U_bold(G)) " with " epsilon.alt_plus.minus^0 = underbrace(epsilon.alt_(1/2 bold(G))^0, "free energy" #linebreak() "at edge") plus.minus U_bold(G), outset: #.3em)
$
see Figure ...

#pagebreak()
= Tight-binding model
We've already discussed multiple ways to calculate energy bands. The last one we treat in this course is the tight-binding method.

We start with two seperated neutral atoms, take hydrogen, and consider what happens when they are brought together to form a crystal. Say the atoms have wavefunctions $psi_A$ and $psi_B$. As they are brought together we consider the linear combinations $psi_A plus.minus psi_B$. An electron in the $+$ state will have lower energy than an electron in the $-$ state. Since in the $+$ state it will sometimes be in the middle of the atoms, where there is a stronger binding energy, this never happens in the $-$ state. If we instead had $N$ neutral atoms we'd find $N$ orbitals for each orbital of the isolated atom. The Coulomb interaction between atom cores and the electron leads to energy splitting and energy bands.

This approximation based on the free atom wavefunctions is called LCAO or tight-binding, and works pretty well for inner orbitals. We assume the ground state of an electron moving in $U(bold(r))$ of a free atom is $phi(bold(r))$. Taking the influence of one atom on another to be small we approximate the wavefunction of the crystal to be a sum of these
$
  psi_bold(k) (bold(r)) = sum_j C_(bold(k) j) phi(bold(r)-bold(r)_j)
$
where we sum over the entire lattice. For a crystal of $N$ atoms this can be written in Bloch form if $C_(bold(k) j) = N^(-1\/2) e^(i bold(k) dot bold(r)_j)$
$
  psi_bold(k) (bold(r)) = 1/sqrt(N) sum_j e^(i bold(k) dot bold(r)_j) phi(bold(r)-bold(r)_j)
$
#proof[
  Consider a translation $bold(T)$ connecting any two lattice points
  $
    psi_bold(k) (bold(r)+bold(T)) &= N^(-1\/2) sum_j e^(i bold(k) dot bold(r)_j) phi(bold(r)+bold(T)-bold(r)_j) \
    &= e^(i bold(k) dot bold(T)) N^(-1\/2) sum_j e^(i bold(k) dot (bold(r)_j-bold(T))) phi(bold(r)-(bold(r)_j-bold(T))) \
    &= e^(i bold(k) dot bold(T)) psi_bold(k) (bold(r))
  $
  by the Bloch condition.
]

We can then find the first-order energy by calculating $braket(bold(k), H, bold(k))$ (for $s$-orbitals)
$
  braket(bold(k), H, bold(k)) &= 1/N sum_(j,m) e^(i bold(k) dot overbrace((r_j-r_m), -bold(rho)_m)) braket(phi_m, H, phi_j) \
  &=^(sum_j = N) sum_m e^(- i bold(k) dot bold(rho)_m) braket(phi(bold(r)-bold(rho)_m-r_j), H, phi(bold(r)-bold(r)_j)) \
  &=^"translational invariance" sum_m e^(- i bold(k) dot bold(rho)_m) braket(phi(bold(r)-bold(rho)_m), H, phi(bold(r)))
$
this is quite nice. To evaluate the matrix elements we assume only nearest neighbors affect eachother, we define
$
  braket(phi(bold(r)), H, phi(bold(r))) = -alpha";  " braket(phi(bold(r)-bold(rho)), H, phi(bold(r))) = - gamma
$
with $bold(rho)$ being nearest neighbor distance, then
$
  braket(bold(k), H, bold(k)) &= underbracket(braket(phi(bold(r)), H, phi(bold(r))), -alpha) + sum_"nearest neighbors" e^(-i bold(k) dot bold(rho)_m) underbracket(braket(phi(bold(r)-bold(rho)_m), H, phi(bold(r))), -gamma) \
  &= - alpha - gamma sum_m e^(-i bold(k) dot bold(rho)_m) = epsilon.alt_bold(k)
$
for a simple cubic structure
$
  bold(rho)_m = (plus.minus a,0,0)";  " (0,plus.minus a,0)";  " (0,0,plus.minus a)
$
so
$
  epsilon_bold(k) = -alpha -2gamma (cos k_x a + cos k_y a + cos k_z a)
$
so the energies are confined to a band of width $12 gamma$.

#alt[

  What we've done above is essentially used the variational principle
  $
    braket(bold(k), H, bold(k)) >= E_"gs"
  $
  and taken the interaction with nearest neighbors as a pertubative term
  $
    H = underbracket(H_0, "free-electron sum") + V
  $
  so we could've written (still assuming only neighbor terms matter)
  $
    braket(bold(k), H, bold(k)) = underbracket(braket(phi, H_0, phi), epsilon.alt_0)+underbracket(braket(phi, V, phi), -beta) + sum_"nearest" e^(-i bold(k) dot bold(rho)_m) [underbracket(braket(phi_bold(rho), H_0, phi), epsilon.alt_0 braket(phi_bold(rho), phi))+underbracket(braket(phi_bold(rho), V, phi), -t)]
  $
  since we take the wavefunctions to be highly localized then $braket(phi_bold(rho), phi) tilde 0$ so we get
  $
    braket(bold(k), H, bold(k)) = underbracket(epsilon.alt_0 - beta, -alpha) - t sum_"nearest" e^(-i bold(k) dot bold(rho)_m) = epsilon.alt_bold(k)
  $
  which is basically the same as above. This more clearly shows that in the limit $V -> 0$ we just get back
  $epsilon.alt_bold(k) = epsilon.alt_0$.
]

