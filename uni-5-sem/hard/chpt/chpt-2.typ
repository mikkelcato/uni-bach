//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Phonons
Phonons are essentially quantized elastic waves in a crystal, as we'll see later. These are really important for the heat capacity and heat conduction in a crystal which is why we care.

We start discussing general vibrations in crystals, then we quantize these, and discuss phonon properties.

== Vibrations of Crystals
=== Monatomic basis
We start by considering a crystal with just one atom in the primitive  cell---we want to find the frequency of an elastic wave in this crystal. This is most simply done for waves with propagation direction $[100], [110] "or" [111]$ in a cubic crystal. In this case the problem becomes one-dimensional with planes of atoms displacing in phase either parallel (longitudinal---one mode) or perpendicular (transverse---two modes) to the wave vector, with the displacement of a plane $s$ being described by some $u_s$.

We assume the response is a linear function of the force, considering just the nearest interactions we can write
$
  F_s = C (u_(s+1)-u_s) + C(u_(s-1)-u_s)
$
where $C$ is the force constant, and we treat this as the force per atom. Then the equation of motion for an atom in the plane is
$
  M dv(u_s, t, 2) = C (u_(s+1) + u_(s-1) - 2 u_s)
$
we assume a time-dependence of the form $exp(-i omega t)$ so
$
  - M omega^2 u_s = C(u_(s+1)+u_(s-1)-2u_s)
$
this has solutions
$
  u_(s plus.minus 1) = u exp(i s K a) exp(plus.minus i K a)
$
with $a$ being the plane spacing. Substituting this immediately gives the dispersion relation
$
  omega^2 = (2C)/M (1-cos K a)
$
In our setup we essentially have a one-dimensional lattice so the lattice vectors are given by $n a$ with $n in ZZ$. This means the reciprocal lattice vectors are just
$
  G = (2 pi)/a m
$
so the distance between any two planes is $2 pi\/a$ making the midpoints $plus.minus pi/a$ so the boundary of the first Brillouin zone has $K = plus.minus pi\/a$. We see that at the boundary we have
$
  dv(omega^2, K) = 0
$
we can also write
$
  omega^2 = (4 C)/M sin^2 (K a)/2
$
we basically only care about the first Brillouin zone since
$
  u_(s+1)/u = exp(i K a) => - pi < K a <= pi => -pi/a < K <= pi/a
$
so values of $K$ outside $K_max = plus.minus pi\/a$ reproduce the lattice---since if we have some $K$ outside the first Brillouin zone we should be able to find some $K' = K - 2pi n\/a$. At the boundaries we have
$
  u_s = u exp(plus.minus i s pi) = u (-1)^s
$
this is just a standing wave, since alternating atoms oscillate with opposite phases, since alternating atoms oscillate with opposite phases. $K_"max"$ satisfies the Bragg condition with $lambda = 2 a$.

With the dispersion relation we can derive many things, take e.g. the group velocity
$
  v_g = nabla_bold(K) omega (bold(K))
$
which in this case becomes
$
  v_g = ((C a^2)/M)^(1\/2) cos (K a)/2
$
as we'd expect for a standing wave then this is zero at the boundary of the first Brillouin zone. In the long wavelength limit $K a << 1$ we can Taylor expand $cos K a = 1 - 1/2 (K a)^2$ giving the relation
$
  omega^2 = (C/M) K^2 a^2 => omega = sqrt((C)/M) K a
$
so the velocity is independent of the frequency---so $v = omega\/K$ which is why we sometimes call the limit the continuum limit.

=== Diatomic basis
In the case of two atoms in the primitive cell each polarization mode develops two branches---the acoustical and optical, so we have longitudinal acoustical (LA) and optical (LO) modes etc. With $p$ atoms we get $3 p$ branches, $3$ acoustical branches (3 spatial directions, everything moves together) and $3 p - 3$ optical branches.

We now consider a cubic crystal with interleaved planes consisting of atoms with $M_1$ and $M_2$ respectively with planar spacing $a$ (between identical planes)---we should be able to find a symmetry direction such that a single plane contains just one type of atom. Assuming only nearest-neighbor interactions we can find
$
  M_1 dv(u_s, t, 2) & = C (v_s + v_(s-1) - 2 u_s) \
  M_2 dv(v_s, t, 2) & = C (u_(s+1) + u_s - 2 v_s)
$
with $u_s$ and $v_s$ being different amplitudes. We try
$
  u_s & = u exp(i s K a) exp(- i omega t) \
  v_s & = v exp(i s K a) exp(- i omega t)
$
on substitution
$
  - omega^2 M_1 u & = C v [1 + exp(-i K a)]-2 C u \
  - omega^2 M_2 v & = C u [exp(i K a) + 1] - 2 C v
$
this has a solution when
$
  matrixdet(2 C - M_1 omega^2, - C [1+exp(- i K a)]; -C[1+exp(i K a)], 2 C - M_2 omega^2) = 0
$
or
$
  M_1 M_2 omega^4 - 2 C (M_1 + M_2) omega^2 + 2 C^2 (1-cos K a) = 0
$
for $K a << 1$
$
  M_1 M_2 omega^4 - 2 C (M_1 + M_2) omega^2 + C^2 K^2 a^2 = 0
$
giving the roots
$
  omega^2 &= (2 C(M_1+M_2))/(2 M_1 M_2) plus.minus 1/(2 M_1 M_2) (4 C^2 (M_1+M_2)^2 - 4 C^2 M_1 M_2 K^2 a^2)^(1\/2) \
  &= (C(M_1 + M_2))/(M_1 M_2) plus.minus (C(M_1+M_2))/(M_1 M_2) [1 - (M_1 M_2)/(M_1+M_2)^2 K^2 a^2]^(1\/2) \
  &= (C (M_1+M_2))/(M_1 M_2) plus.minus (C(M_1+M_2))/(M_1 M_2) (1 - (M_1M_2)/(M_1+M_2)^2 (K^2 a^2)/2) \
  &= C {1/M_1 + 1/M_2 plus.minus [1/M_1+1/M_2 - 1/(M_1+M_2) (K^2 a^2)/2]}
$
giving
$
  omega^2 & = 2 C (1/M_1+1/M_2) - C/(M_1+M_2) (K^2 a^2)/2 " (optical)" \
  omega^2 & = C/(M_1+M_2) (K^2 a^2)/2 " (acoustic)"
$
at $K_"max"$ we find
$
  omega^2 = (2 C)/M_1 "and" omega^2 = (2 C)/M_2
$
by back substitution of the optical branch at $K = 0$ we can find
$
  u/v = - M_2/M_1
$
so the atoms vibrate against each other, but their center of mass is fixed. Using the other branch in the limit $K = 0$ we find $u = v$, so the atoms move together as in a typical acoustical vibration. Wavelike solutions don't exist between $sqrt(2C\/M_1)$ and $sqrt(2 C\/M_2)$---there is a frequency gap at the boundary.

== Quantization of Elastic Waves
The energy of these vibrations is quantized. The quantum of this energy is called a phonon. The energy of a mode is given by
$
  epsilon = (n + 1/2) hbar omega
$
similarly to the simple harmonic oscillator.

=== Quantization of a Ring
To see this we model phonons as vibrations of a linear lattice of particles connected by springs. Let $N$ particles of mass $M$ be connected by springs with spring constant $C$ and length $a$---we let these form a ring. Now consider the transverse displacement of the particles out of the plane of this ring. We denote the displacement of particle $s$ by $q_s$ and its momentum by $p_s$. The Hamiltonian is
$
  H = sum_(s=1)^N {1/(2 M) p_s^2 + 1/2 C (q_(s+1)-q_s)^2}
$
this is just a sum of simple harmonic oscillators. To solve this Hamiltonian we Fourier transform to $P_k$ and $Q_k$ which we call phonon coordinates. We let
$
  q_s = N^(-1\/2) sum_k Q_k exp(i k s a) <-> Q_k = N^(-1\/2) sum_s q_s exp(-i k s a)
$
note that we have the boundary condition $q_s = q_(s+N)$ giving
$
  k = (2 pi n)/(N a) " with " n = 0, plus.minus 1, dots, 1/2 N
$
so $-pi\/a < k <= pi\/a$. The momentum $P_k$ is canonically conjugate to $Q_k$ so
$
  p_s = N^(-1\/2) sum_k P_k exp(-i k s a) <-> P_k = N^(-1\/2) sum_s p_s exp(i k s a)
$
these satisfy the canonical commutation relations
$
  [Q_k, P_k'] & = N^(-1) [sum_r q_r exp(- i k r a), sum_s p_s exp(i k' s a)] \
              & = N^(-1) sum_(r,s) [q_r,p_s] exp[-i(k r- k' s) a] \
              & = N^(-1) i hbar sum_r exp[-i(k-k') r a] \
              & = i hbar delta_(k,k')
$
transforming the Hamiltonian and using
$
  sum_s p_s^2 & = N^(-1) sum_(s,k,k') P_k P_k' exp[-i(k+k')s a] \
              & = sum_(k,k') P_k P_k' delta_(-k,k') \
              & = sum_k P_k P_(-k)
$
and
$
  sum_s (q_(s+1)-q_s)^2 & = 2 sum_k Q_k Q_(-k) (1- cos k a)
$
gives
$
  H = sum_k {1/(2 M) P_k P_(-k) + C Q_k Q_(-k) (1 - cos k a)}
$
defining
$
  omega_k equiv sqrt((2 C)/M (1 - cos k a))
$
we can write
$
  H = sum_k {1/(2 M) P_k P_(-k) + 1/2 M omega_k^2 Q_k Q_(-k)}
$
in the Heisenberg picture we can find
$
  i hbar dot(Q)_k = [Q_k, H] = (i hbar P_(-k))/M => i hbar dot.double(Q)_k = [dot(Q)_k, H] = i hbar omega_k^2 Q_k
$
so
$
  dot.double(Q)_k + omega_k^2 Q_k = 0
$
this is just a simple harmonic oscillator with energy
$
  epsilon_k = (n_k + 1/2) hbar omega_k
$
then the energy of the entire ring is
$
  U = sum_k (n_k + 1/2) hbar omega_k
$
=== Ladder operators
It is useful to write
$
  H = sum_k hbar omega_k (a^dagger_k a_k + 1/2)
$
with
$
  a^dagger ket(n) = sqrt(n+1) ket(n+1)",  " a ket(n) = sqrt(n) ket(n-1)
$
giving $a^dagger a ket(n) = n ket(n)$---when the phonon mode $k$ is in the eigenstate $ket(n_k)$ we say that there are $n_k$ phonons in the mode. The eigenvalues of the Hamiltonian then immediately gives the previous $U$.

Note that since $a a^dagger ket(n) = (n+1) ket(n)$ then $[a, a^dagger] = 1$. We can write
$
  Q_k & = sqrt(hbar/(2 M omega_k)) (a_k + a^dagger_(-k)) \
  P_k & = i sqrt((hbar M omega_k)/2) (a_k^dagger - a_(-k))
$
giving
$
  q_s = sum_k sqrt(hbar/(2 N M omega_k)) [a_k exp(i k s) + a_k^dagger exp(- i k s)]
$
note that $Q^dagger_(-k) = Q_k$ and $P_k^dagger = P_(-k)$. Inerting the relations for $Q_k$ and $P_k$ we can show
$
  [a_k, a_k^dagger] = delta_(k,k')
$
note $omega_k = omega_(-k)$.

=== Back to vibrations
Consider a standing wave move with amplitude
$
  u = u_0 cos K x cos omega t
$
with $u$ being the displacement of a volume element from the equilibrium position at $x$. The kinetic energy density is $1/2 rho(dd(u, d: partial)\/dd(t, d: partial))^2$ the volume integral of this guy is
$
  1/4 rho V omega^2 u_0^2 sin^2 omega t -->^"average" 1/8 rho V omega^2 u_0^2 = 1/2 (n+1/2) hbar omega
$
since for a harmonic oscillator the energy is split equally between kinetic and potential energy (when averaged). This gives
$
  u_0^2 = (4 (n+1/2) hbar)/(rho V omega)
$
this nicely relates $n$ to $u_0$.

A phonon with some $K$ would interact with stuff as if it had momentum $hbar K$---even though they don't carry physical momentum. We saw previously that scattering of a photon in a crystal is governed by
$
  bold(k)' = bold(k) + bold(G)
$
in this case the crystal would recoil with momentum $-hbar bold(G)$. If the scattering is inelastic, and creates a phonon, then
$
  bold(k)' + bold(K) = bold(k) + bold(G)
$
or if a phonon is absorbed then
$
  bold(k)' = bold(k) + bold(K) + bold(G)
$

== Heat Capacity
We define heat capacity as
$
  C_V equiv (pdv(U, T))_V
$
the phonon contribution is denoted by $C_"lat"$. The total phonon energy at some $tau equiv k_B T$ in a crystal can be written as a sum
$
  U_"lat" = sum_(K,p) U_(K,p) = sum_(K,p) expval(n_(K,p)) hbar omega_(K,p)
$
with $expval(n_(K,p))$ having the obvious meaning---it's form is given by
$
  expval(n) = 1/(exp(hbar omega\/tau) -1)
$

#proof[
  By the Boltzmann distribution
  $
    N_(n+1)/N_n = exp(- hbar omega\/tau)
  $
  and
  $
    P_n = N_n/(sum_(s=0)^oo N_s) = exp(- n hbar omega\/tau)/(sum_(s=0)^oo exp(-s hbar omega\/tau))
  $
  so
  $
    expval(n) = sum_s s P_s = (sum_s s exp(- s hbar omega \/ tau))/(sum_s exp(-s hbar omega \/ tau)) = 1/(exp(hbar omega \/tau) -1)
  $
]

Then
$
  U = sum_(K,p) (hbar omega_(K,p))/(exp(hbar omega_(K,p)\/tau)-1)
$
we replace the summation over $K$ by an integral
$
  U = sum_p integral dd(omega) D_p (omega) (hbar omega)/(exp(hbar omega\/tau)-1)
$
where $D_p (omega)$ is the density of states---it tells us the number of modes per unit frequeny range. Then
$
  C_"lat" = k_B sum_p integral dd(omega) D_p (omega) (x^2 exp x)/(exp x - 1)^2
$
with $x = hbar omega \/ tau$.

=== Finding $D_p (omega)$.
It should be obvious that this is the central problem.

==== One dimension
Consider a one-dimensional line of length $L$ carrying $N+1$ particles with separation $a$. We take the particles $s = 0$ and $s = N$ at the ends to be fixed. Each mode of polarization $p$ is a standing wave with
$
  u_s = u_0 exp(- i omega_(K,p) t) sin s K a
$
by the boundary condition
$
  K = pi/L, (2 pi)/L, dots, ((N-1)pi)/L
$
each of these is associated with a standing wave, importantly there is one mode for each interval $Delta K = pi\/L$ so the number of modes is $N = K\/Delta K = (L\/pi) K$. And we obtain
$
  D_1 (omega) dd(omega) = dv(N, omega) dd(omega) = L/pi dv(K, omega) dd(omega) = L/pi dd(omega)/(dd(omega)\/dd(K))
$
This could also be done using periodic boundary conditions.

==== Two dimensions
For a square lattice with periodic boundary conditions there is one allowed value of $K$ per area $(2 pi\/L)^2$ so within a circle of $pi K^2$ we have
$
  N(K) = L^2/(2pi)^2 pi K^2 = (L^2 K^2)/(4 pi) -> dv(N, K) = (L^2 K)/(2 pi)
$

==== Three dimensions
We apply periodic boundaries over $N^3$ primitive cells within a cube of side $L$. So
$
  exp[i(K_x x + K_y y + K_z z)] = exp{i [K_x (x+L) + K_y (y+L) + K_z (z+L)]}
$
giving
$
  K_i = 0, plus.minus (2 pi)/L, dots, (N pi)/L
$
so we have one $bold(K)$ value per $(2 pi \/L)^3$, so
$
  N = (4/3 pi K^3)/(Delta K)^3 = V/(6 pi^2) K^3
$
giving for each polarization
$
  D(omega) = dv(N, omega) = (V K^2)/(2 pi^2) dv(K, omega)
$

=== Debye model
In the Debye approximation we write
$
  omega = v K
$
then
$
  D(omega) = (V omega^2)/(2 pi^2 v^3)
$
with $N$ primitive cells we have $N$ acoustic phonon modes. We define the cutoff frequency by
$
  N = V/(6 pi^2) K^3 = V/(6 pi^2 v^3) omega_D^3 => omega_D^3 = (6 pi^2 v^3 N)/V
$
this has a corresponding $K_D = omega_D\/v$. In the Debye model we don't allow modes with $K > K_D$, since $K <= K_D$ exhausts the degrees of freedom of a monatomic lattice.

Then
$
  U = integral dd(omega) D(omega) expval(n(omega)) hbar omega = integral_0^omega_D dd(omega) ((V omega^2)/(2 pi^2 v^3))((hbar omega)/(exp(hbar omega\/tau)-1))
$
for each $p$. Assuming the velocity doesn't care about $p$ we multiply by a factor $3$ to obtain
$
  U = (3 V k_B^4 T^4)/(2 pi^2 v^3 hbar^3) integral_0^x_D dd(x) x^3/(exp x- 1)
$
with $x equiv hbar omega \/ tau$ and $x_D equiv hbar omega_D \/ k_B T equiv theta\/T$. With
$
  theta = (hbar v)/k_B ((6 pi^2 N)/V)^(1\/3)
$
being the Debye temperature. Then we can write the energy as
$
  U = 9 N k_B T (T/theta)^3 integral_0^x_D dd(x) x^3/(exp x - 1)
$
and the heat capacity can also easily be found
$
  C_V = 9 N k_B (T/theta)^3 integral_0^x_D dd(x) (x^4 exp x)/(exp x -1)^2
$
at low temperatures we can approximate the energy by letting $x_D -> oo$ giving
$
  U = (3 pi^4 N k_B T^4)/(5 theta^3)
$
since the integral just becomes a Bose-Einstein integral---and
$
  C_V = (12 pi^4)/5 N k_B (T/theta)^3
$
which is the Debye approximation for the heat capacity---which needs really low temperature else it kind of sucks.

=== Einstein model
We consider $N$ oscillators with the same $omega_0$ in one dimension. The Einstein density of states is then $D(omega) = N delta(omega-omega_0)$---the energy is then
$
  U = N expval(n) hbar omega = (N hbar omega)/(exp(hbar omega\/tau)-1)
$
with $omega_0 -> omega$ for niceness. The heat capacity becomes
$
  C_V = N k_B ((hbar omega)/tau)^2 exp(hbar omega\/tau)/[exp(hbar omega\/tau)-1]^2
$
in three dimensions we just replace $N -> 3N$ (three modes per oscillator). Then the high temperatur limit is $3 N k_B$ which is the Dulong-Petit law---at low temperature this model kind of sucks.

=== General $D(omega)$
In general we can write
$
  D(omega) dd(omega) & = (L/(2 pi))^3 integral_"shell" dd(K, 3) \
                     & = (L/(2 pi))^3 integral dd(S_omega) dd(omega)/v_g \
                     & => D(omega) = V/(2pi)^3 integral dd(S_omega)/v_g
$
where the integral is taken over the area of the surface $omega =$ constant in $bold(K)$ space.

=== Anharmonic effects
In the harmonic approximation phonons never die and lattice waves don't interact---many things are generally much simpler. We consider collisions next which are anharmonic by definition.

== Thermal Conductivity
The thermal conductivity coefficient $K$ is defined by (steady-state flow of heat down a long rod)
$
  j_U = - K dv(T, x)
$
where $j_U$ is the flux of thermal energy---that it depends on $nabla T$ implies this is a somewhat random process. From kinetic theory we know that
$
  K = 1/3 C v cal(l)
$
with $C$ being heat capacity per unit volume, $v$ being the average particle velocity and $cal(l)$ being the mean free path.

#proof[
  The flux in the $x$-direction is $n expval(abs(v_x))$ with $n$ being the concentration of particles. Given a particle has heat capacity $c$ then when moving from $T + Delta T -> T$ it gives up energy $c Delta T$---and
  $
    Delta T = dv(T, x) cal(l)_x = dv(T, x) v_x tau
  $
  then it follows that the flux of energy is
  $
    j_U = - n expval(v_x^2) c tau dv(T, x) = - 1/3 n expval(v^2) c tau dv(T, x)
  $
  for phonons $v$ is constant so
  $
    j_U = - 1/3 C v cal(l) dv(T, x) " with " cal(l) equiv v tau, C equiv n c
  $
  so $K = 1/3 C v cal(l)$.
]

=== Normal process
The $cal(l)$ for phonons is determined by scattering, due to geometry and other phonons. Anharmonic coupling between phonons predicts $cal(l) prop 1\/T$ at high temperatures (if the forces were harmonic, then there'd be no mechanism for collisions between phonons)---this can be understood by assuming that the number of excited phonons is $prop T$, then the collision frequency should be the inverse. To actually define conductivity there must be some mechanism which can bring a distribution of phonons into thermal equilibrium (locally), else we'd be unable to say e.g. one end of a rod has $T_1$ while the other has $T_2$. Phonon collisions with geometry will not establish this, since $omega_f = omega_i$ in such collisions. Even three-phonon collisions (an N-process)
$
  bold(K)_1 + bold(K)_2 = bold(K)_3
$
is not enough---since the total momentum is unchanged. An equilibrium distribution with $T$ can move down the crystal with some drift velocity which is not changed by N-processes. The phonon momentum
$
  bold(J) = sum_bold(K) n_bold(K) hbar bold(K)
$
is conserved. If $bold(J) eq.not 0$ then these collisions are unable to establish equilibrium since they leave $bold(J)$ unchanged---there is no thermal resistance, since the distribtuion of phonons will propagate without changing $bold(J)$.

\* splitting

=== Umklapp process
These processes have the form
$
  bold(K)_1 + bold(K)_2 = bold(K)_3 + bold(G)
$
where $bold(G)$ is some reciprocal lattice vector---here energy is conserved. We have already seen wave-interaction in crystals where the total wavevector is changed by some reciprocal lattice vector---in particularly for phonons we know everything meaningful happens in the first Brillouin zone, so any longer $bold(K)$ must be brought back by addition of some $bold(G)$. If $bold(G) = 0$ the process becomes normal.

At high temperatures $T > theta$, all phonon modes are excited since $k_B T > hbar omega_"max"$. Then most processes will be U-processes with high momentum change---by the previous argument it is expected that $cal(l) prop 1\/T$.

The energy of phonons $bold(K)_1, bold(K)_2$ suitable for a U-process is of order $1/2 k_B theta$ since each should have wavevectors of order $1/2 G$. If both phonons have low energy, so small $K$, then they can't get outside the first Brillouin zone. The number of suitable phonons with energy $1/2 k_B theta$ is given by $tilde exp(- theta\/2 T)$ from the Boltzmann factor. The needed $cal(l)$ is the mean free path for umklapp collisions between phonons.

=== Geometry
These effects include distribution of mass, edge effects, imperfections etc. At low $T$ the mean free path $cal(l)$ becomes $tilde$ size of sample. Then $cal(l)$ obviously becomes limited, and the conductivity becomes $prop "size of sample"$---at the same time the umklapp process becomes less effective at limiting conductivity. We obtain
$
  K tilde C v D
$
for $cal(l) tilde D$, and at low temperatures $C tilde T^3$. In a real material with defects we instead have
$
  cal(l) tilde 1/root(3, N_D)
$
with $N_D$ being the density of defects.

#pagebreak()
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

== Electrical conductivity
=== Drude model
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

#pagebreak()
= Energy Bands
== Nearly free $e^-$ model
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

We extend this to the nearly free electron model, wherein band electrons are treated as being pertubed weakly by the periodic potential of the ion cores. Bragg reflection of electron waves then occurs, and causes the band gaps. As an example consider a linear lattice with lattice constant $a$, then the Bragg condition becomes
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

=== Zone schemes
All $bold(k)$ live in reciprocal space as we've seen and any $bold(k)' = bold(k) + bold(G)$ can be translated back to the first Brillouin zone, so we can write
$
  epsilon.alt(bold(k)') = hbar^2/(2m) (bold(k)+bold(G))^2
$
so all energies out side the first Brillouin zone can be represented inside the it---this is the reduced zone scheme, and bands get folded within the first Brillouin zone, we basically still calculate free electron parabolas. Usually this is done using a direction e.g. from $Gamma:(0,0,0) -> X:(pi/a,0,0)$, so we'd only vary $k_x$.

This is the main idea in the empty lattice approximation, in terms of later stuff this corresponds to setting $V=0$.

== Bloch's theorem
Bloch proved that solutions to the Schrödinger equation for a periodic potential takes the form
$
  psi_bold(k) (bold(r)) = u_bold(k) (bold(r)) exp(i bold(k) dot bold(r))
$
with $u_bold(k) (bold(r)) = u_bold(k) (bold(r)+bold(T))$. So any solution $psi_bold(k)$ takes the form of a free wave with some periodic modulation $u_bold(k)$---these are sometimes called Bloch functions. We'll prove a restatement of this later.

=== Kronig-Penney
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

== Central equation
We now treat the wave equation for a general potential, at general $k$. Let $U(x)$ be the energy of an electron in a linear lattice with lattice constant $a$, we know $U(x)=U(x+a)$. We can then expand it as
$
  U(x) = sum_G U_G exp(i G x)
$
the $U_G$ tend to decrease rapidly. We rewrite this to ensure $U(x) in RR$
$
  U(x) = sum_(G>0) U_G (exp(i G x) + exp(-i G x)) = 2 sum_(G>0) U_G cos G x
$
we assume $U(x)$ is symmetric about $x=0$ and $U_0 = 0$. The Schrödinger equation becomes
$
  (p^2/(2 m) + sum_G U_G exp(i G x)) psi(x) = epsilon.alt psi(x)
$
this described one electron in the potential of ion cores, and in the average potential of other electrons. We can similarly expand $psi$ as
$
  psi = sum_k C_k exp(i k x)
$
where $k$ ensure this satisfies boundary conditions, these have the form $2 pi n \/L$. Note that we claim nothing about the periodicity of $psi(x)$, the properties of $psi(x)$ with respect to translation are fully determined by Bloch's theorem. Further not all $k$ of the given form enter the expansion, since if any $k$ is in the expansion then all $k + G$ enter aswell---for this reason when labeling the Bloch functions we pick the $k$ within the first Brillouin zone.

Substituting both expansions into the wave equation gives
$
  sum_k hbar^2/(2 m) k^2 C_k exp(i k x) + sum_(G,k) U_G C_k exp(i(k+G) x) = epsilon.alt sum_k C_k exp(i k x)
$
by matching we find
$
  (lambda_k - epsilon.alt) C_k + sum_G U_G C_(k-G) = 0
$
where
$
  lambda_k = (hbar^2 k^2)/(2 m)
$
this is called the central equation.

Now we could write
$
  psi_k (x) & = sum_G C_(k-G) exp(i(k-G)x) \
            & = (sum_G C_(k-G) exp(-i G x)) exp(i k x) = exp(i k x) u_k (x)
$
where
$
  u_k (x) equiv sum_G C_(k-G) exp(-i G x)
$
this satisfies $u_k (x + T) = u_k (x)$
$
  u_k (x + T) & = sum C_(k-G) exp(-i G (x+T)) \
              & =exp(-i G T) [sum C_(k-G) exp(-i G x)] \
              & = exp(-i G T) u_k (x)
$
and $exp(-i G T) = 1$ by definition. So we have proven Bloch's theorem.

==== Properties of $bold(k)$
We'll shortly state some properties of $bold(k)$ which acts as a label for the Bloch functions.

Under a translation we have
$
  psi_bold(k) (bold(r)+bold(T)) = exp(i bold(k) dot bold(T)) psi_bold(k) (bold(r))
$
so $exp(i bold(k) dot bold(T))$ is the phase factor we multiply our Bloch functions by when we do a translation. If the lattice potential vanishes $U_G = 0$ then
$
  (lambda_bold(k) - epsilon.alt) C_bold(k) = 0
$
so all $C_(bold(k)-bold(G))$ vanish except for $C_bold(k)$, meaning $u_bold(k)$ is constant or $psi_bold(k) (bold(r)) = exp(i bold(k) dot bold(r))$. We recover the free electron.

This $bold(k)$ also enters in all conservation laws, which is why $hbar bold(k)$ is called the crystal momentum in the first place. Take the example of an electron absorbing a phonon with momentum $bold(q)$, then
$
  bold(k)+ bold(q) = bold(k)' + bold(G)
$

=== Examples
==== Kronig-Penney
We use
$
  U(x) = 2 sum_(G > 0) U_G cos G x = A a sum_s delta(x - s a)
$
we sum over all integers $0 < s < a^(-1)$. We have periodicity over a ring of unit length ($a^(-1)$ atoms) so
$
  U_G = integral_0^1 dd(x) U(x) cos G x = A a sum_s integral_0^1 dd(x) delta(x-s a) cos G x = A a sum_s cos G s a = A
$
then the central equation becomes
$
  (lambda_k - epsilon.alt) C_k + A sum_n C(k - (2 pi n)/a) = 0
$
define
$
  f(k) = sum_n C(k - (2 pi n)/a)
$
then
$
  C_k = - ((2 m A\/hbar^2) f(k))/(k^2 - (2 m epsilon.alt\/hbar^2))
$
we sum over all $n$ so
$
  f(k) = f(k - (2 pi n)/a)
$
meaning we can write
$
  C(k - (2 pi n)/a) = - ((2 m A)/hbar^2) f(k) [(k-(2 pi n)/a)^2 - (2 m epsilon.alt)/hbar^2]^(-1)
$
and summing over all $n$ we obtain
$
  hbar^2/(2 m A) = - sum_n [(k-(2pi n)/a)^2 - ((2 m epsilon.alt)/hbar^2)]^(-1)
$
using
$
  "ctn" x = sum_n 1/(n pi + x)
$
and $K^2 = 2 m epsilon.alt\/hbar^2$ we finally get
$
  (m A a^2)/(2 hbar^2) 1/(K a) sin K a + cos K a = cos k a
$
agreeing with the previous result.

==== Near zone boundary
We take $U_G$ to be small relative to the kinetic energy of an electron at the zone boundary. We consider a wavevector at the boundary $1/2 G$ or $pi\/a$ then
$
  k^2 = (1/2 G)^2;"  "(k-G)^2 = (1/2 G-G)^2 = (1/2 G)^2
$
so at the boundary the kinetic energy of the two waves $k=plus.minus 1/2 G$ are equal.

If $C(1/2 G)$ is important then so is $C(- 1/2 G)$, we take these to be the important ones---in NFE we assume the crystal potential is weak meaning other contributions are negligible. Then we find
$
  (lambda - epsilon.alt) C (1/2 G) + U C(-1/2 G) & = 0 \
   (lambda - epsilon.alt) C(-1/2 G) + U C(1/2 G) & = 0
$
giving
$
  (lambda-epsilon.alt)^2 = U^2;"  "epsilon.alt = lambda plus.minus U = hbar^2/(2 m) (1/2 G)^2 plus.minus U
$
so the energy has two roots, and we've found an energy gap of $E_g = 2 U$ at the boundary. We find
$
  C(-1/2 G)/C(1/2 G) = (epsilon.alt - lambda)/U = plus.minus 1
$
so
$
  psi(x) = exp((i G x)/2) plus.minus exp((-i G x)/2)
$
these are identical to what we started with---one solution gives the wavefunction below $E_g$ and vice versa depending on the sign of $U$.

Now we solve for these in terms of $k$ (to see behavior around the edge)
$
  psi(x) = C(k) exp(i k x) + C(k-G) exp(i(k-G)x)
$
by the central equation
$
    (lambda_k - epsilon.alt) C(k) + U C(k-G) & = 0 \
  (lambda_(k-G)-epsilon.alt) C(k-G) + U C(k) & = 0
$
giving
$
  epsilon.alt = 1/2 (lambda_(k-G)+lambda_k) plus.minus [(lambda_(k-G)-lambda_k)^2/4 + U^2]^(1\/2)
$
each root describes a band. We expand this in terms of $tilde(K) equiv k - 1/2 G$ (assumed to be small) giving
$
  epsilon.alt_tilde(K) &= hbar^2/(2 m) (G^2/4 + tilde(K)^2) plus.minus [4 lambda ((hbar^2 tilde(K)^2)/(2 m)) + U^2]^(1\/2) \
  &tilde.eq hbar^2/(2 m) (G^2/4 + tilde(K)^2) plus.minus U [1+2(lambda/U^2) (hbar^2 tilde(K)^2)/(2 m)]
$
for $hbar^2 G tilde(K)\/2 m << abs(U)$. We can write this as
$
  epsilon.alt_tilde(K) (plus.minus) = epsilon.alt(plus.minus) + (hbar^2 tilde(K)^2)/(2 m) ( 1 plus.minus (2 lambda)/U)
$
with $epsilon.alt (plus.minus)$ being the energies at the boundary.

== Bands
We consider a linear crystal constructed of an even number $N$ of primitive cells with lattice constant $a$. By periodicity
$
  k = 0; plus.minus (2 pi)/L; dots (N pi)/L
$
we cut this off at $N pi\/L = pi\/a$ being the zone boundary. Each primitive cell contributes one independent value of $k$ to each energy band. Taking into account spin there are $2 N$ independent orbitals in each band---which carries over to three dimensions. If there is a single atom of valence, one per primitive cell, then the band can be half filled. If each atom contributes two valence electrons, then the band can be exactly filled. Similarly if there are two atoms of valence, one in each primitive cell.

=== Metals and insulators
If one or more bands are exactly filled, then the crystal will be an insulator, provided there is some $E_g$ to the next higher (empty) band. This can only occur if the number of valence electrons in the primitive cell is even---and even then band overlap can lead to two partly filled bands, making the crystal a metal (or semimetal). Alkali metals and noble metals have one valence electron per primitive cell, and are thus always metals.
