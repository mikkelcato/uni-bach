//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()


= Binding
We want to know why crystals are held together. There are four principal types of binding: the van der Waals interaction, the Coulomb interaction between charged particles, covalent binding due to hybridization and metallic binding. We will go through each in turn.

== van der Waals interaction
The van der Waals interaction is present everywhere and form the simplest crystals. Since it is present everywhere this is the principle interaction responsible for the binding of inert crystals.

The van der Waals interaction is caused by small distortions in electron distributions leading to induced dipoles. The interaction is given by
$
  dd(U, d: Delta) = - A/R^6 -> "very weak"
$
with $A tilde hbar omega_0 alpha^2$. This is proven below.

#proof[

  We consider two harmonic oscillators separated by a distance $R$. Each has two charges $plus.minus e$ with separations $x_1$ and $x_2$. These oscillate along $x$ with momenta $p_1$ and $p_2$. We obtain the Hamiltonian
  $
    H_0 = p_1^2/(2 m) + 1/2 m omega_0^2 x_1^2 + p_2^2/(2 m) + 1/2 m omega_0^2 x_2^2
  $
  Let $H_1$ be the Hamiltonian for the Coulomb interaction
  $
    H_1 = e^2/R + e^2/(R+x_1-x_2) - e^2/(R+x_1) - e^2/(R-x_2)
  $
  taking the dipole to be perfect ($x_1 << R$ and $x_2 << R$) we find
  $
    H_1 tilde.eq - (2 e^2 x_1 x_2)/R^3
  $
  We use the substitution
  $
    x_s equiv 1/sqrt(2) (x_1 + x_2)",  " x_a equiv 1/sqrt(2) (x_1 - x_2)
  $
  to write
  $
    H = [p_s^2/(2 m) + 1/2 (m omega_0^2-(2 e^2)/R_3) x_s^2] + [p_a^2/(2 m) + 1/2 (m omega_0^2+ (2 e^2)/R^3) x_a^2]
  $
  We find two coupled modes $s$ and $a$ with
  $
    omega_(a,s) = [(m omega_0^2 plus.minus (2 e^2)\/R^3)/m]^(1\/2) = omega_0 [1 plus.minus 1/2 ((2 e^2)/(m omega_0^2 R^3)) - 1/8 ((2 e^2)/(m omega_0^2 R^3))^2 + dots]
  $
  the zero-point energy is $ (hbar (omega_s + omega_a))/2 $ which is lower than the uncoupled value $hbar omega_0$ by
  $
    dd(U, d: Delta) = - hbar omega_0 1/8 ((2 e^2)/(m omega_0^2 R^3))^2 equiv - A/R^6
  $
  as we wanted.

]

=== Pauli repulsion
When atoms come together their electron distribution begin to overlap. Due to the Pauli exclusion principle this leads to a repulsive interaction since electrons are forced to occupy high energy states.

There are multiple ways to model this interaction, one way is taking it to be of the form $B R^(-12)$. Combined with the van der Waals interaction we find
$
  U(R) = 4 epsilon [(sigma/R)^12 - (sigma/R)^6]
$
which is called the Lennard-Jones potential. Later we use an interaction of the form $lambda e^(-R/rho)$ with $rho$ being the interaction range.

The total energy of a crystal with $N$ atoms can be approximated by summing over the Lennard-Jones potential
$
  U_"tot" = underbracket(N/2, "all pairs") #h(1em) underbracket(4 epsilon [sum_j (sigma/(p_(i j)R))^12 - sum_j (sigma/(p_(i j)R))^6], "energy due to one atom")
$
with $p_(i j) R$ being the distance between atom $i$ and any other atom $j$ in terms of the nearest neighbor distance $R$. For the face-centered cubic we find the equilibrium distance to be
$
  dv(U_"tot", R) = 0 => R_0/sigma = 1.09
$
then the cohesive energy is approximately given by
$
  U_"tot" (R_0) = -(2.15)(4 N epsilon)
$


== Coulomb interaction
The principle interaction between ions with $plus.minus q$ is the Coulomb interaction $tilde plus.minus q^2\/r$ and is thus responsible for binding of ionic crystals.

The van der Waals interaction is always present but for ionic crystals it is very weak. The main contribution to the binding energy is instead the Madelung energy. The energy associated to a reference ion $i$ is
$
  U_i = sum_j U_(i j)
$
where we sum over all $j eq.not i$. We write $U_(i j)$ as a combination of the Pauli repulsion and the Coulomb interaction
$
  U_(i j) = lambda e^(- r_(i j)/rho) plus.minus q^2/r_(i j)
$
letting $r_(i j) equiv p_(i j) R$ and only including the Pauli repulsion for nearest neighbors we obtain
$
  U_(i j) = cases(
    lambda e^(-R/rho)- q^2/R #h(15pt) & "nearest",
    plus.minus q^2/(p_(i j) R) & "else"
  )
$
where we assume the nearest-neighbors have opposite charge. The total energy for $2N$ ions is then
$
  U_"tot" = N U_i = N (z lambda e^(-R\/rho) - (alpha q^2)/R)
$
where $z$ is the number of nearest neighbors and $alpha$ is the Madelung constant
$
  alpha equiv sum_j (plus.minus 1)/p_(i j)
$
The equilibrium separation is given by
$
  N dv(U_i, R) = - (N z lambda)/rho e^(-R\/rho) + (N alpha q^2)/R^2 = 0 => R_0^2 e^(-R_0\/rho) = (rho alpha q^2)/(z lambda)
$
so we obtain
$
  U_"tot" = underbracket(- (N alpha q^2)/R_0, "Madelung energy") (1- rho/R_0)
$
and we see that for a crystal to be stable $alpha$ must be positive. Taking the reference ion $i$ to be negative then $plus$ applies to all positive ions. An equivalent definition of $alpha$ is
$
  alpha/R = sum_j (plus.minus 1)/r_(i j)
$
with $r_(i j)$ being the distance from the $j$th ion to the reference ion.

As an example consider an infinite line of alternating ions with $R$ being the distance between adjacent ions. Then
$
  alpha/R = 2 [1/R - 1/(2 R) + 1/(3 R) - 1/(4 R) + dots] => alpha = 2 [1-1/2+1/3-dots] = 2 ln 2
$
this is obviously way harder for three-dimensional structures.

== The exchange interaction
What the exchange interaction means is not immediately obvious. So we illustrate it with an example.

Consider the hydrogen molecule. When separated we have two protons each with their own electron. As they are brought close their orbitals can _hybridize_. When hybridized the electrons are localized in the region between the two atoms and have anti-parallel spin (it is the singlet state). This is the symmetric (bonding) state. We could also imagine that each atom tries to hold on to their own electron. This is the anti-symmetric (anti-bonding) state. By considering the Hamiltonian one can compute
$
  J = expval(H)_"anti" - expval(H)_"sym" tilde 10 "eV"
$
this is the exchange splitting and it tells us that it is favorable to be in the symmetric state.

Due to symmetry these covalent bonds are highly directional. Since the exchange splitting is high they are also very strong.


== Other interactions
Metallic bonding and hydrogen bonding.

#pagebreak()
= Elasticity
== Hooke's law
We approximate a crystal has a continuous medium. Then the relationship between stress and strain is given by Hooke's law
$
  bold(sigma) =^"eng" C bold(e) =^"math" C bold(epsilon)
$
with $C$ being the stiffness tensor, or equivalently
$
  bold(epsilon) = S bold(sigma)
$
with $S$ being the compliance tensor. The stress $bold(sigma)$ describes the response to some applied strain $bold(epsilon)$. Both $C$ and $S$ have $36$ components.

We consider displacing the crystal by
$
  bold(R) (bold(r)) = u(bold(r)) hat(x) + v (bold(r)) hat(y) + w (bold(r)) hat(z)
$
then the strain components are defined by
$
  epsilon_(x x) = pdv(u, x) "etc."
$
we will need the related quantities
$
  e_(i i) = epsilon_(i i)";  " e_(i j) = epsilon_(j i) + epsilon_(i j)
$

=== Elastic energy
To determine the various components we consider the elastic energy
$
  U = 1/2 sum_(lambda=1)^6 sum_(mu=1)^6 tilde(C)_(lambda mu) e_lambda e_mu
$
with $lambda, mu = {x x, y y, z z, y z, z x, x y}$. The stress components are defined by
$
  sigma_(x x) = pdv(U, e_(x x)) equiv pdv(U, e_1) = tilde(C)_(1 1) e_1 + 1/2 sum_(beta = 2)^6 (tilde(C)_(1 beta) + tilde(C)_(beta 1)) e_beta
$
we see
$
  C_(alpha beta) = 1/2 (tilde(C)_(alpha beta) + tilde(C)_(beta alpha)) = C_(beta alpha)
$
so $C$ and $S$ are symmetric tensors.

Given our crystal has some symmetry the situation is even simpler. As an example for a cubic crystal we claim
$
  U = 1/2 C_11 (e_(x x)^2 + e_(y y)^2 + e_(z z)^2) + 1/2 C_44 (e_(y z)^2 + e_(z x)^2 + e_(x y)^2) + C_12 (e_(y y) e_(z z) + e_(z z) e_(x x) + e_(x x) e_(y y))
$
In a cubic crystal we have four three-fold symmetric rotational axes, these act as
$
  x -> y -> z -> x & "  " -x->z->y->-x \
    x-> z-> -y-> x & "  " -x->y->z->-x
$
all terms are invariant under these in the claimed $U$. So we have just three stiffness constants in a cubic crystal. The specific surviving $C_(alpha beta)$ are determined by taking derivatives and comparing with $bold(sigma) = C bold(e)$.

=== Bulk modulus
Consider the case $e_(x x)=e_(y y)=e_(z z) = 1/3 delta$ then
$
  U = 1/6 (C_11 + 2 C_12) delta^2
$
and we define the bulk modulus $B$ by
$
  U = 1/2 B delta^2
$
which is equivalent to
$
  B = - V dv(p, V)
$
so for a cubic crystal
$
  B = 1/3 (C_11+2C_12)
$
the compressibility $K$ is defined as $K = B^(-1)$.

== The wave equation
Elastic waves in a crystal evolve according to
$
  rho pdv(u, t, 2) = pdv(sigma_(x x), x)+pdv(sigma_(x y), y)+pdv(sigma_(x z), z)
$
we can write this as
$
  rho pdv(u, t, 2) &= C_11 pdv(e_(x x), x) + C_12 (pdv(e_(y y), x)+pdv(e_(z z), x)) + C_44 (pdv(e_(x y), y)+pdv(e_(z x), z)) \
  &= C_11 pdv(u, x, 2) + C_44 (pdv(u, y, 2)+pdv(u, z, 2)) + (C_12+C_44) (pdv(v, x, y)+pdv(w, x, z))
$
by symmetry we have similar equations for $v$ and $w$.

=== Simple solutions
Consider the ansatz
$
  u = u_0 e^(i(K x-omega t)) "  longitudinal wave"
$
with $u$ being the $x$-component of the displacement and $bold(K) = K hat(x)$. We obtain the dispersion
$
  omega^2 rho = C_11 K^2
$
so the velocity of a longitudinal wave in the $[1 0 0]$ direction is
$
  v_s = nu lambda = omega/K = sqrt(C_11/rho)
$
Consider instead the ansatz
$
  v = v_0 exp[i(K x- omega t)] "  transverse wave"
$
with $v$ being the $y$-component of the displacement. We obtain
$
  omega^2 rho = C_44 K^2 => v_s = sqrt(C_44/rho)
$
which is the velocity for a transverse wave in the $[100]$ direction. Similarly can be found for $w$.

For more complicated $bold(K)$ one finds a system of equations. As an example take $bold(K) = K_x hat(x) + K_y hat(y)$ so the direction $[110]$. Consider the ansatz
$
  w = w_0 e^(i(K_x x + K_y y - omega t)) "  transverse"
$
we find
$
  omega^2 rho = C_44 (K_x^2 + K_y^2) = C_44 K^2
$
which is similar to before. Consider the ansatz
$
  u = u_0 e^(i(K_x x+ K_y y- omega t))";  " v = v_0 e^(i(K_x x + K_y y - omega t)) "  longitudinal/shear"
$
here we need both since $bold(K)$ has multiple directions. We find a system of equations
$
  omega^2 rho u & = (C_11 K_x^2 + C_44 K_y^2) u + (C_12 + C_44) K_x K_y v \
  omega^2 rho v & = (C_11 K_y^2 + C_44 K_x^2) v + (C_12 + C_44) K_x K_y u
$
in the longitudinal case $u = v$ with $K_x = K_y = K\/sqrt(2)$ we find
$
  omega^2 rho = 1/2 (C_11 + C_12 + 2C_44) K^2
$
in the shear case $u = - v$ we find
$
  omega^2 rho = 1/2 (C_11 - C_12) K^2
$


#pagebreak()
= Phonons
== Monoatomic basis
We consider a crystal with one atom in the primitive cell. We want the dispersion relation for an elastic wave in such a crystal. We take the wave to propagate along a direction with a certain amount of symmetry. This makes the problem one-dimensional with planes of atoms displacing in phase. This displacement can be longitudinal (one mode) or transverse (two modes).


We now use the harmonic approximation and assume the atoms in our crystal are connected by similar springs. Denoting the displacement of a plane $s$ by $u_s$ we can then write
$
  F_s = C (u_(s+1)-u_s) + C(u_(s-1)-u_s)
$
where $C$ is the spring constant, and we only take into account nearest-neighbor interactions. The equation of motion for an atom in the plane is then
$
  M dv(u_s, t, 2) = C (u_(s+1) + u_(s-1) - 2 u_s)
$
We use the ansatz
$
  u_s = u_0 e^(i(K a n - omega t))
$
with $a$ being the plane spacing. Using this we obtain the dispersion
$
  omega^2 = (2C)/M (1-cos K a) = (4 C)/M sin^2 (K a)/2
$
this is quite nice!

The lattice vectors are given by $n a$ with $n in ZZ$. Meaning
$
  G = (2 pi)/a m
$
so the distance between any two planes is $2 pi a^(-1)$. The boundary of the first Brillouin zone is then $K = plus.minus pi a^(-1)$. We only care about the first Brillouin zone since we can bring any $K'$ inside it by $K' -> K - 2 pi n a^(-1)$.

Consider the group velocity
$
  v_g & = nabla_bold(K) omega (bold(K)) \
      & = ((C a^2)/M)^(1\/2) cos (K a)/2
$
at $K = plus.minus pi a^(-1)$ we get standing waves! Taking the long wavelength limit $K a << 1$ (the continuum limit) we find
$
  omega = sqrt((C a^2)/M) K equiv v_s K
$
and we recover a linear dispersion relation.

Consider a finite system with $N$ atoms. By periodicity
$
  u_s =^! u_(s+N)
$
meaning
$
  e^(i K N a) = 1 => K = (2 pi n)/(N a)
$
with $n = 0, plus.minus 1 dots, N\/2$ since $K_"max" = pi a^(-1)$. We find that the allowed values of $K$ are discrete. So the elastic wave is only defined at discrete points!

For a monoatomic basis with $N$ atoms we have $3 N$ degrees of freedom.


== Diatomic basis
We consider a crystal with two atoms in the primitive basis. We again only consider directions where entire planes displacing in phase.

Since we have two atoms each polarization mode (one longitudinal and two transverse) develop two branches. These are called the acoustical and optical. For the acoustical branches the two atoms move together, so they act as if we just had one atom giving three branches. For the optical branch they move out of phase, so for two atoms this gives three optical branches. This generalizes to $p$ atoms, here we still have three acoustical branches, but now we have $3 p - 3$ optical branches. So in total we have $3 p$ branches. The optical branches are more energetic since they oscillate with larger frequency. The degrees of freedom are always concerved so if we have two atoms in the primitive cell we always get six branches, though they are typically degenerate.


We take our crystal to have planes of atoms with alternating masses $M_1$ and $M_2$. We should be able to find a direction with enough symmetry such that a single plane contains just one type of atom. Assuming only nearest-neighbors matter we find
$
  M_1 dv(u_s, t, 2) & = C (v_s + v_(s-1) - 2 u_s) \
  M_2 dv(v_s, t, 2) & = C (u_(s+1) + u_s - 2 v_s)
$
with $u_s$ and $v_s$ being different amplitudes. We use the ansatz
$
  u_s = u e^(i s K a) e^(- i omega t)";  " v_s & = v e^(i s K a) e^(- i omega t)
$
giving the set of equations
$
  - omega^2 M_1 u & = C v (1 + e^(-i K a))-2 C u \
  - omega^2 M_2 v & = C u (e^(i K a) + 1) - 2 C v
$
the solution is given by
$
  matrixdet(2 C - M_1 omega^2, - C (1+e^(- i K a)); -C(1+e^(i K a)), 2 C - M_2 omega^2) = 0
$
giving
$
  omega^2 = C/(M_1 M_2) [M_1+ M_2 plus.minus sqrt((M_1+M_2)^2 -2M_1 M_2 (1- cos K a))]
$
Taking the long wavelength limit we find
$
  omega^2 & tilde.eq 2 C (1/M_1+1/M_2) " (optical)" \
  omega^2 & tilde.eq C/(M_1+M_2) (K^2 a^2)/2 " (acoustical)"
$
at $K_"max"$ we find
$
  omega^2 = (2 C)/M_1 "and" omega^2 = (2 C)/M_2
$
By back-substitution one can show these are actually optical and acoustical modes.

== Quantization
The energy of elastic waves in a crystal is quantized, this is what we call a phonon!

The energy of a mode is given by
$
  epsilon.alt_k = hbar omega_k (n_k + 1/2)
$
where
$
  omega_k^2 equiv (2 C)/M (1 - cos k a)
$
and $n_k$ counts the number of phonons with wavevector $k$.

#proof[

  Let $N$ particles of mass $M$ be connected by springs with spring constant $C$ and length $a$. The Hamiltonian is
  $
    H = sum_(s=1)^N {p_s^2/(2 M) + 1/2 C (q_(s+1)-q_s)^2}
  $
  this is just a sum of simple harmonic oscillators.

  We now Fourier transform to the phonon coordinates $P_k$ and $Q_k$ by
  $
    q_s &= N^(-1\/2) sum_k Q_k e^(i k s a) <-> Q_k = N^(-1\/2) sum_s q_s e^(-i k s a) \
    p_s &= N^(-1\/2) sum_k P_k e^(-i k s a) <-> P_k = N^(-1\/2) sum_s p_s e^(i k s a)
  $
  these satisfy
  $
    [Q_k, P_k'] = i hbar delta_(k,k')
  $
  The Hamiltonian becomes
  $
    H = sum_k {1/(2 M) P_k P_(-k) + C Q_k Q_(-k) (1 - cos k a)}
  $
  defining
  $
    omega_k^2 equiv (2 C)/M (1 - cos k a)
  $
  we can write this as
  $
    H = sum_k {(P_k P_(-k))/(2 M) + 1/2 M omega_k^2 Q_k Q_(-k)}
  $
  The Heisenberg equations of motion give
  $
           i hbar dot(Q)_k & = [Q_k, H] = (i hbar P_(-k))/M \
    i hbar dot.double(Q)_k & = [dot(Q)_k, H] = i hbar omega_k^2 Q_k
  $
  We find $Q_k$ is just a harmonic oscillator
  $
    dot.double(Q)_k + omega_k^2 Q_k = 0
  $
  and we are done.

]

A phonon with $K$ interacts as if it had momentum $hbar K$ and energy $hbar omega_k$. We saw previously that scattering of a photon in a crystal is determined by
$
  bold(k)' = bold(k) + bold(G)
$
Given the scattering is inelastic it will create a phonon with $bold(K)$
$
  bold(k)' + bold(K) = bold(k) + bold(G)
$
or if a phonon is absorbed then
$
  bold(k)' = bold(k) + bold(K) + bold(G)
$

== Heat capacity
We define heat capacity as
$
  C_V equiv (pdv(U, T))_V
$
the phonon contribution is denoted by $C_"lat"$.

The total phonon energy at some $tau equiv k_B T$ can be written as
$
  U_"lat" = sum_(K,p) U_(K,p) = sum_(K,p) expval(n_(K,p)) hbar omega_(K,p)
$
with $expval(n_(K,p))$ given by the Planck distribution
$
  expval(n) = 1/(e^(hbar omega\/tau) -1)
$

#proof[

  By the Boltzmann distribution
  $
    N_(n+1)/N_n = e^(- hbar omega\/tau)
  $
  and
  $
    P_n = N_n/(sum_(s=0)^oo N_s) = e^(- n hbar omega\/tau)/(sum_(s=0)^oo e^(-s hbar omega\/tau))
  $
  so
  $
    expval(n) = sum_s s P_s = (sum_s s e^(- s hbar omega \/ tau))/(sum_s e^(-s hbar omega \/ tau)) = 1/(e^(hbar omega \/tau) -1)
  $
]

Then
$
  U_"lat" = sum_(K,p) (hbar omega_(K,p))/(e^(hbar omega_(K,p)\/tau)-1)
$
we replace the summation over $K$ by an integral
$
  U_"lat" = sum_p integral dd(omega) D_p (omega) (hbar omega)/(exp(hbar omega\/tau)-1)
$
where $D_p (omega)$ is the density of states. Then
$
  C_"lat" = k_B sum_p integral dd(omega) D_p (omega) (x^2 exp x)/(exp x - 1)^2
$
with $x = hbar omega tau^(-1)$.

=== Finding $D_p (omega)$.
==== One dimension
Consider a one-dimensional line of length $L$ carrying $N+1$ particles with separation $a$. We take the particles $s = 0$ and $s = N$ to be fixed. Each mode of polarization $p$ is a standing wave with
$
  u_s = u_0 e^(- i omega_(K,p) t) sin s K a
$
by the boundary condition
$
  K = pi/L, (2 pi)/L, dots, ((N-1)pi)/L
$
each of these is associated with a standing wave. There is one mode for each interval $Delta K = pi\/L$ so the number of modes is $N = K\/Delta K = (L\/pi) K$. We obtain
$
  D_p^"1D" (omega) dd(omega) = dv(N, omega) dd(omega) = L/pi dv(K, omega) dd(omega) = L/pi 1/v_g dd(omega)
$

==== Two dimensions
For a square lattice with periodic boundary conditions there is one allowed value of $K$ per area $(2 pi\/L)^2$ so within a circle of $pi K^2$ we have
$
  N(K) = L^2/(2pi)^2 pi K^2 = (L^2 K^2)/(4 pi) -> dv(N, K) = (L^2 K)/(2 pi)
$
so we obtain
$
  D_p^"2D" (omega) = L^2/(2 pi) 1/v_g K(omega)
$

==== Three dimensions
We apply periodic boundaries over $N^3$ primitive cells within a cube of side $L$. So
$
  e^(i(K_x x + K_y y + K_z z)) = e^(i (K_x (x+L) + K_y (y+L) + K_z (z+L)))
$
giving
$
  K_i = 0, plus.minus (2 pi)/L, dots, (N pi)/L
$
so we have one $bold(K)$ value per $(2 pi \/L)^3$. Then
$
  N = (4/3 pi K^3)/(Delta K)^3 = V/(6 pi^2) K^3
$
so we obtain
$
  D_p^"3D" (omega) = dv(N, omega) = (V)/(2 pi^2) 1/v_g K^2 (omega)
$

=== Debye model
We assume the Debye approximation is valid $omega = v K$. This is reasonable in the low temperature limit where long wavelengths dominate. Then
$
  D(omega) = (V omega^2)/(2 pi^2 v^3)
$
With $N$ primitive cells we have $N$ acoustic phonon modes, so we define the cutoff frequency
$
  N = V/(6 pi^2) K^3 = V/(6 pi^2 v^3) omega_D^3 => omega_D^3 = (6 pi^2 v^3 N)/V
$
this has a corresponding $K_D = omega_D v^(-1)$. The Debye model does not allow modes with $K > K_D$ since $K <= K_D$ exhausts the lattice.

Then
$
  U_"lat" &= integral dd(omega) D(omega) expval(n(omega)) hbar omega \
  &= integral_0^omega_D dd(omega) ((V omega^2)/(2 pi^2 v^3))((hbar omega)/(e^(hbar omega\/tau)-1))
$
for each $p$. Assuming the velocity does not care about $p$ we multiply by a factor $3$ to obtain
$
  U_"lat" = (3 V k_B^4 T^4)/(2 pi^2 v^3 hbar^3) integral_0^x_D dd(x) x^3/(exp x- 1)
$
with $x equiv hbar omega tau^(-1)$ and $x_D equiv hbar omega_D tau^(-1) equiv theta T^(-1)$. Where we introduce the Debye temperature
$
  theta = (hbar v)/k_B ((6 pi^2 N)/V)^(1\/3)
$
this separates the quantum and classical regime. Then the energy is
$
  U = 9 N k_B T (T/theta)^3 integral_0^x_D dd(x) x^3/(exp x - 1)
$
and the heat capacity can be found
$
  C_V = 9 N k_B (T/theta)^3 integral_0^x_D dd(x) (x^4 exp x)/(exp x -1)^2
$
At low temperatures we let $x_D -> oo$ giving
$
    U & = (3 pi^4 N k_B T^4)/(5 theta^3) \
  C_V & = (12 pi^4)/5 N k_B (T/theta)^3
$
this is the Debye heat capacity.

=== Einstein model
We consider $N$ oscillators with the same $omega_0$ in one dimension. The Einstein density of states is $D(omega) = N delta(omega-omega_0)$ giving
$
  U = N expval(n) hbar omega = (N hbar omega)/(e^(hbar omega\/tau)-1)
$
with $omega_0 -> omega$ for niceness. We find
$
  C_V = N k_B ((hbar omega)/tau)^2 e^(hbar omega\/tau)/(e^(hbar omega\/tau)-1)^2
$
for one dimension, for three dimensions we replace $N -> 3N$. The high temperature limit is $3 N k_B$.


== Thermal conductivity
The thermal conductivity $K$ is defined by
$
  j_U = - K dv(T, x)
$
where $j_U$ is the flux of thermal energy. From kinetic theory we know that
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

We now consider processes that limit $cal(l)$ for phonons. Consider two phonons scattering, by momentum conservation we have
$
  bold(K)_1 + bold(K)_2 + bold(G) = bold(K)
$
the scattering is normal if $bold(G) = 0$ and Umklapp if $bold(G) eq.not 0$. Consider the total phonon momentum
$
  bold(J) = sum_bold(K) n_bold(K) hbar bold(K)
$
in thermal equilibrium $expval(bold(J)) =^! 0$. By the condition above normal scattering events preserve $bold(J)$. But for Umklapp scattering we have $dd(bold(J), d: Delta) = hbar bold(G)$ and this manifests itself as thermal resistance.

For an Umklapp process to occur one of the involved phonons must have $bold(K) tilde bold(K)_D$ meaning it has energy $tilde bold(omega)_D$. The number of phonons with this energy is
$
  expval(n_(bold(K))) tilde e^(- theta/T)
$
so for $T gt.tilde theta$ most scattering will be Umklapp meaning to $cal(l) tilde T^(-1)$. While for $T lt.tilde theta$ we have exponential suppresion $cal(l) tilde tau tilde e^(- theta/T)$.

Then assuming the velocity is independent of temperature $K tilde C cal(l)$ we find for high temperatures $K tilde T^(-1)$ for intermediate temperatures $K tilde e^(-theta/T)$ and for low temperatures $K tilde T^3$.

Besides scattering $cal(l)$ can also be limited by the geometry of our sample. At low $T$ we have $cal(l) tilde$ size of sample. Then $cal(l)$ obviously becomes limited and the conductivity becomes $prop "size of sample"$. At the same time the Umklapp scattering becomes less effective at limiting conductivity. We obtain
$
  K tilde T^3 v D
$
for $cal(l) tilde D$ and $C tilde T^3$. For a real material with defects we instead have
$
  cal(l) tilde 1/root(3, N_D)
$
with $N_D$ being the density of defects.
