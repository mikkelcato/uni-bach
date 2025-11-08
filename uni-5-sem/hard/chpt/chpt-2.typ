//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#text(size: 35pt, strong("Elastic waves & Phonons"))
= Crystal Binding and Elastic Constants
What holds a crystal together?---in solids the electrostatic interaction between electrons and nuclei is entirely responsible for cohesion. There are four principal types of crystalline binding---van der Waals, ionic, metallic and covalent.

== Inert gases
These form the simplest crystals and the $e^-$ distribution is similar to the free atoms. In this case the atoms pack as tight as possible, so most are fcc. Tiny distortions of the $e^-$ distribution leads to the van der Waals interaction---small fluctuations lead to induced dipoles.

To model this consider two identical harmonic oscillators separated by $R$. Each has charges $plus.minus e$ with separations $x_1$ and $x_2$---these oscillate along $x$, let $p_1$ and $p_2$ be their momenta. Then the Hamiltonian
$
  cal(H)_0 = 1/(2 m) p_1^2 + 1/2 C x_1^2 + 1/(2 m) p_2^2 + 1/2 C x_2^2
$
we assume each has $omega_0$, then $C = m omega_0^2$. Let $cal(H)_1$ be the Coulomb interaction then
$
  cal(H)_1 = e^2/R + e^2/(R+x_1-x_2) - e^2/(R+x_1) - e^2/(R-x_2)
$
in the point-dipole limit $abs(x_1), abs(x_2) << R$ we find to lowest order
$
  cal(H)_1 tilde.equiv - (2 e^2 x_1 x_2)/R^3
$
writing
$
  x_s equiv 1/sqrt(2) (x_1 + x_2)",  " x_a equiv 1/sqrt(2) (x_1 - x_2)
$
we can obtain
$
  x_1 = 1/sqrt(2) (x_s + x_a)",  " x_2 = 1/sqrt(2) (x_s - x_a)
$
and equivalent expressions for $p_1$ and $p_2$---subbing these in the total Hamiltonian is
$
  cal(H) = [1/(2 m) p_s^2 + 1/2 (C-(2 e^2)/R_3) x_s^2] + [1/(2 m) p_a^2 + 1/2 (C+ (2 e^2)/R^3) x_a^2]
$
so have two coupled modes $s$ and $a$, by inspection
$
  omega = [(C plus.minus (2 e^2)\/R^3)/m]^(1\/2) = omega_0 [1 plus.minus 1/2 ((2 e^2)/(C R^3)) - 1/8 ((2 e^2)/(C R^3))^2 + dots]
$
the zero-point energy is $hbar\/2 (omega_s + omega_a)$ which is lower than the uncoupled value $hbar omega_0$ by
$
  dd(U, d: Delta) = 1/2 hbar(dd(omega_s, d: Delta)+dd(omega_a, d: Delta)) = - hbar omega_0 1/8 ((2 e^2)/(C R^3))^2 = - A/R^6
$
this is the van der Waals interaction or London interaction or induced dipole-dipole interaction$dots$ It is the principle attractive interaction in crystals of inert gases---and is solely a product of the dipole-dipole coupling---one can approximate $A tilde hbar omega_0 alpha^2$.

As two atoms get attracted they will begin to overlap, this lead to a repulsive force due to the Pauli exclusion principle---$e^-$ will be forced to higher states increasing the energy and leading to a repulsive interaction $tilde B\/R^12$. Combining these effects one can write the potential energy of two atoms as
$
  U(R) = 4 epsilon [(sigma/R)^12 - (sigma/R)^6]
$
this is the Lennard-Jones potential---other exponential forms $lambda exp(-R\/rho)$ with $rho$ being a measure of the interaction range can also be used and are usually easier to handle analytically.

We can approximate the total energy of a crystal with $N$ atoms by summing over the Lennard-Jones potential
$
  U_"tot" = 1/2 N(4 epsilon) [sum_j ' (sigma/(p_(i j)R))^12 - sum_j ' (sigma/(p_(i j)R))^6]
$
with $p_(i j) R$ being the distance between atom $i$ and any other atom $j$ in terms of the nearest neighbor distance $R$---the sums are known. For fcc we can find the equilibrium distance
$
  dv(U_"tot", R) = 0 => R_0/sigma = 1.09
$
using this one can then obtain
$
  U_"tot" (R_0) = -(2.15)(4 N epsilon)
$
which should hold for all inert gases---this is the cohesive energy when the atoms are at rest---of course this is also very naive and ignores quantum effects which become more evident for smaller atoms.

== Ionic crystals
Ionic crystals consist of positive and negative ions---the ionic bond comes from the electrostatic interaction of oppositely charged ions---as with inert gas atoms we expect that the ions have symmetric charge distributions, but with distortions where they touch.

The long range interaction between ions with $plus.minus q$ is the Coulomb interaction $plus.minus q^2\/r$---the ions will try to balance this with the repulsive interaction between ion cores, this is similar to the inert gas case. The van der Waals interaction is very weak for ionic bonds, instead the main contribution to the binding energy is the Madelung energy. We define
$
  U_i = sum_j ' U_(i j)
$
where we sum over all $j eq.not i$. We write $U_(i j)$ as a central field repulsive potential (instead of the Pauli potential) and a Coulomb interaction
$
  U_(i j) = lambda exp(- r_(i j)/rho) plus.minus q^2/r_(i j)
$
neglecting surface effects we write $U_"tot" = N U_i$ with $N$ being the amount of molecules---$2 N$ ions. We again use $r_(i j) equiv p_(i j) R$, only including the repulsive interaction for nearest neighbors we find
$
  U_(i j) = cases(
    lambda exp(-R\/rho)- q^2\/R #h(10pt) & "nearest",
    plus.minus 1\/p_(i j) q^2\/R & "else"
  )
$
so
$
  U_"tot" = N U_i = N (z lambda e^(-R\/rho) - (alpha q^2)/R)
$
where $z$ is the number of nearest neighbors of any ion and
$
  alpha equiv sum_j ' ((plus.minus))/p_(i j) equiv "Madelung constant"
$
the equilibrium separation can be written as
$
  N dv(U_i, R) = - (N z lambda)/rho e^(-R\/rho) + (N alpha q^2)/R^2 = 0 => R_0^2 e^(-R_0\/rho) = (rho alpha q^2)/(z lambda)
$
using this we obtain
$
  U_"tot" = - (N alpha q^2)/R_0 (1- rho/R_0)
$
with the Madelung energy being
$
  - (N alpha q^2)/R_0
$
for a crystal to be stable we require that $alpha$ is positive. If we take our reference ion to be negative then the plus sign applies to all positive ions---an equivalent definition of $alpha$ is
$
  alpha/R = sum_j ' ((plus.minus))/r_j
$
with $r_j$ being the distance from the $j^"th"$ ion to the reference ion. As an example consider an infinite line of alternating ions---let $R$ be the distance between adjacent ions, then
$
  alpha/R = 2 [1/R - 1/(2 R) + 1/(3 R) - 1/(4 R) + dots] => alpha = 2 (1-1/2+1/3-dots) = 2 ln 2
$
this is obviously way harder for three-dimensional structures.

== Other bonds
Aside from the mentioned bonds we quickly mention covalent bonds, metallic bonds and hydrogen bonds.

Covalent bonds are formed from two $e^-$---one from each atom. These tend to be partly localized between the two atoms and their spins when bonded are antiparallel. Notably this bond has strong directional properties and doesn't fill space as tightly as other bonds---it only allows four nearest neighbors. The binding energy depends on the spin orientation of the two $e^-$ due to the Pauli exclusion principle which will modify the charge distribution in accordance with the spin distribution. This spin-dependent Coulomb interaction is called the exchange interaction. The big point is that orbitals will hybridize, and sometimes this is favorably and a bond will form (bonding)---other times this is not at all favorably and a bond will not form (antibonding).

Metals have a large number of free $e^-$ zooming around---conduction $e^-$---the valence $e^-$ of the atoms become the conduction $e^-$ of the metal. In a metallic bond the energy of these valence $e^-$ is lower as compared with the free atom.

== Elasticity
We briefly discuss elasticity since we'll eventually discuss elastic waves in crystals.

We treat a crystal as a continuous homogeneous medium---this approximation is valid for elastic waves of $lambda$ longer that $tilde 10^(-8) "m"$. We'll only consider infinitesimal strains such that Hooke's law applies. Take three axis defined by $hat(bold(x))_i$ then after a small uniform deformation these get deformed
$
  bold(x)' &= (1 + epsilon.alt_(x x)) hat(bold(x)) + epsilon.alt_(x y) hat(bold(y)) + epsilon.alt_(x z) hat(bold(z)) \
  bold(y)' &= epsilon.alt_(y x) hat(bold(x)) + (1 + epsilon.alt_(y y)) hat(bold(y)) + epsilon.alt_(y z) hat(bold(z)) \
  bold(z)' &= epsilon.alt_(z x) hat(bold(x)) + epsilon.alt_(z y) hat(bold(y)) + (1 + epsilon.alt_(z z)) hat(bold(z))
$
the coefficients $epsilon.alt_(alpha beta)$ define the deformation---even though the old axes were of unit length this is no longer guaranteed. We define the displacement $bold(R)$ of the deformation as
$
  bold(R) equiv bold(r)'-bold(r) = x (bold(x)'-hat(bold(x))) + dots
$
or more generally
$
  bold(R) (bold(r)) = u(bold(r)) hat(bold(x)) + v(bold(r)) hat(bold(y)) + w(bold(r)) hat(bold(z))
$
with
$
  epsilon.alt_(x x) tilde.equiv pdv(u, x)",  " epsilon.alt_(y x) tilde.equiv pdv(u, y)",  etc."
$
we define the strain components by this
$
  e_(x x) equiv epsilon.alt_(x x) = pdv(u, x)
$
similarly for $e_(y y)$ and $e_(z z)$. For the other components we define
$
  e_(x y) equiv bold(x)' dot bold(y)' tilde.equiv epsilon.alt_(y x) + epsilon.alt_(x y) = pdv(u, y) + pdv(v, x)
$

We define the dilation by $V' = bold(x)' dot bold(y)' times bold(z)'$, and we obtain
$
  delta equiv (V'-V)/V tilde.equiv e_(x x)+e_(y y) + e_(z z)
$
we basically ignore terms of order $epsilon.alt^2$ in all of these.
=== Constants
We define the stress as the force acting on a unit area in the solid---we have nine components $X_x,X_y,X_z, dots$, with the capital letter denoting the direction, and the subscript denoting the normal to the plane to which the force is applied---we say that $Y_z = Z_y$ etc. giving us six stress components $X_x, Y_y, Z_z, Y_z, Z_x, X_y$. Assuming Hooke's law holds we can write
$
  e_(x x) = S_(1 1) X_x + S_(1 2) Y_y + S_(1 3) Z_z + S_(1 4) Y_z + S_(1 5) Z_x + S_(1 6) X_y
$
and so on for $e_(y y), e_(z z), e_(y z), e_(z x)$ and $e_(x y)$. Inversely we can write
$
  X_x = C_(11) e_(x x) + C_(12) e_(y y) + C_(13) e_(z z) + C_(14) e_(y z) + C_(15) e_(z x) + C_(16) e_(x y)
$
etc. etc. All the $S_(11) dots$ are called elastic compliance constants, and the $C_(11) dots$ are called elastic stiffness constants or moduli of elasticity.

We can write the energy density as (Hooke's law)
$
  U = 1/2 sum_(lambda=1)^6 sum_(mu=1)^6 tilde(C)_(lambda mu) e_lambda e_mu
$
with $1 -> 6 = {x x, y y, z z, y z, z x, x y}$, with
$
  X_x = pdv(U, e_(x x)) equiv pdv(U, e_1) = tilde(C)_(1 1) e_1 + 1/2 sum_(beta = 2)^6 (tilde(C)_(1 beta) + tilde(C)_(beta 1)) e_beta
$
now we see that
$
  C_(alpha beta) = 1/2 (tilde(C)_(alpha beta) + tilde(C)_(beta alpha)) = C_(beta alpha)
$
i.e. they are symmetric meaning the $36$ constants reduce to $21$.

If our crystal has some symmetry this number reduces further. We claim that for a cubic crystal
$
  U = 1/2 C_11 (e_(x x)^2 + e_(y y)^2 + e_(z z)^2) + 1/2 C_44 (e_(y z)^2 + e_(z x)^2 + e_(x y)^2) + C_12 (e_(y y)e_(z z) + e_(z z) e_(x x) + e_(x x) e_(y y))
$
for a cubic structure we require four three-fold rotation axes, this happens according to
$
  x -> y -> z -> x & "  " -x->z->y->-x \
    x-> z-> -y-> x & "  " -x->y->z->-x
$
all terms are invariant under these in the posed $U$---if a term were odd we could find a rotation which would change the sign since $e_(x y) = - e_(x(-y))$. So we have found that for a cubic crystal we have merely three elastic stiffness constants.

To conclude we can summarize all of this with Hooke's law
$
  bold(sigma) = bold(C) bold(epsilon)
$
with $bold(sigma)$ being all the $X_x dots$ the stress components, and $bold(epsilon)$ being the strain components $e_(x x) dots$. $bold(C)$ is then a tensor containing all the elastic constants.

=== Bulk modulus
Consider the case $e_(x x)=e_(y y)=e_(z z) = 1/3 delta$ then for a cubic crystal
$
  U = 1/6 (C_11 + 2 C_12) delta^2
$
we define the bulk modulus $B$ by
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

=== Elastic waves
For a volume element in a crystal we can obtain
$
  rho pdv(u, t, 2) = pdv(X_x, x)+pdv(X_y, y)+pdv(X_z, z)
$
or plugging stuff in
$
  rho pdv(u, t, 2) &= C_11 pdv(e_(x x), x) + C_12 (pdv(e_(y y), x)+pdv(e_(z z), x)) + C_44 (pdv(e_(x y), y)+pdv(e_(z x), z)) \
  &= C_11 pdv(u, x, 2) + C_44 (pdv(u, y, 2)+pdv(u, z, 2)) + (C_12+C_44) (pdv(v, x, y)+pdv(w, x, z))
$
by symmetry we have similar equations for $v$ and $w$. These equations are hell---we consider one simple solution
$
  u = u_0 exp[i(K x-omega t)]
$
with $u$ being the $x$-component of the particle displacement---with $K=2 pi\/lambda$ and $omega = 2pi nu$. Plugging this guy in gives
$
  omega^2 rho = C_11 K^2
$
so the velocity of a longitudinal wave in the $[1 0 0]$ direction is
$
  v_s = nu lambda = omega/K = (C_11/rho)^(1\/2)
$
take a transverse wave
$
  v = v_0 exp[i(K x- omega t)]
$
and we obtain the velocity of a transverse wave in the $[1 0 0]$ direction
$
  omega^2 rho = C_44 K^2 => v_s = (C_44/rho)^(1\/2)
$
similarly is obtained by $w$. So for $bold(K)$ parallel to $[1 0 0]$ the two shear waves have equal velocities---this is not generally true.

To find the velocities one can also solve a system of equations see the problems---further note that for longitudinal wave the motion is parallel to $bold(K)$ and for transverse motion it is perpendicular to $bold(K)$---in general there are three modes, see problems.

#pagebreak()
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


