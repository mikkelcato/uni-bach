//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *
#import "@preview/mannot:0.3.0": *
#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
= Free electron models
As the name suggests these models assume electrons move freely. We take all scattering to occur with phonons or defects.

== Drude model
We take electrons to be accelerated by a field
$
  m evaluated(dv(bold(v), t))_"field" = - e bold(E)
$
As they move they scatter with other electrons. The amount of scattering is determined by $tau$ and
$
  evaluated(dv(bold(v), t))_"collision" = - bold(v)/tau
$
in steady state
$
  - e bold(E) - (m expval(bold(v)))/tau = 0
$
giving
$
  expval(bold(v)) = - (e bold(E) tau)/m
$
we define the current density by
$
  bold(j) = -n e expval(bold(v)) = sigma bold(E)
$
where we define the Drude conductivity
$
  sigma equiv (n e^2 tau)/m = rho^(-1)
$
this is Ohm's law!

We have $rho = rho_"ph" + rho_i$ since this guy depends on thermal phonons $rho_"ph" (T)$ and impurities $rho_i$. This is called Matthiessen's rule.

Typically $rho$ is dominated by collisions of electrons with lattice phonons giving $rho prop T$.



== Sommerfeld theory
The free electron Hamiltonian is given by
$
  H psi_k = - hbar^2/(2 m) nabla^2 psi_k = epsilon.alt_k psi_k
$
with solutions
$
  psi_bold(k) (bold(r)) = e^(i bold(k) dot bold(r))
$
given that all components of $bold(k)$ satisfy
$
  k_i = 0, plus.minus (2 pi)/L, dots, (2 n pi)/L, dots
$
these ensure periodicity. Then
$
  epsilon.alt_bold(k) = hbar^2/(2 m) k^2 = hbar^2/(2m) (k_x^2 +k_y^2+ k_z^2 )
$
and
$
  bold(p) psi_bold(k) (bold(r)) = - i hbar nabla psi_bold(k) (bold(r)) = hbar bold(k) psi_bold(k) (bold(r))
$
so $psi_bold(k)$ is an eigenfunction of $bold(p)$ with eigenvalue $hbar bold(k)$.

Consider filling all orbitals of a system. This defines some wavevector $k_F$ and sphere in $bold(k)$-space the Fermi surface. The wavevector $k_F$ is defined by
$
  N = overbracket(2, "spin") (4 pi k_F^3\/3)/(2pi \/L)^3 = V/(3 pi^2) k_F^3
$
so all states are filled. The energy at the surface of this sphere is the Fermi energy $epsilon.alt_F$
$
  epsilon.alt_F & = hbar^2/(2 m) k_F^2 \
  epsilon.alt_F & = hbar^2/(2 m) ((3 pi^2 N)/V)^(2\/3)
$
we can further define the electron velocity on the Fermi surface $v_F$
$
  v_F = (hbar k_F)/m = hbar/m ((3 pi^2 N)/V)^(1\/3)
$
and the Fermi temperature $T_F = epsilon.alt_F\/k_B$.

We can also find the density of states
$
  D(epsilon.alt) equiv dv(N, epsilon.alt) = V/(2 pi^2) ((2 m)/hbar^2)^(3\/2) epsilon.alt^(1\/2) = (3 N)/(2 epsilon.alt)
$
where we use
$
  N = V/(3 pi^2) k^3";  " k = sqrt(2 m epsilon.alt)/hbar
$

=== Heat capacity
The Fermi-Dirac distribution tells us that when we heat something up not every electron gains energy $tilde k_B T$. Only those within $tilde k_B T$ of $mu$ can be excited thermally.

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
This gives
$
  dd(U, d: Delta) = integral_(epsilon.alt_F)^oo dd(epsilon.alt) (epsilon.alt-epsilon.alt_F) f(epsilon.alt) D(epsilon.alt) + integral_0^epsilon.alt_F dd(epsilon.alt) (epsilon.alt_F-epsilon.alt) (1-f(epsilon.alt)) D(epsilon.alt)
$
which simplifies to
$
  C_"el" = dv(U, T) = integral_0^oo dd(epsilon.alt) (epsilon.alt-epsilon.alt_F) dv(f, T) D(epsilon.alt)
$
We approximate $D(epsilon.alt) tilde.eq D(epsilon.alt_F)$
$
  C_"el" tilde.eq D(epsilon.alt_F) integral_0^oo dd(epsilon.alt) (epsilon.alt - epsilon.alt_F) dv(f, T)
$
when $k_B T << epsilon.alt_F$ we can ignore the $T$ dependence of $mu$ so we replace $mu -> epsilon.alt_F$. With $x equiv (epsilon.alt-epsilon.alt_F)\/tau$ we write
$
  C_"el" tilde.eq k_B^2 T D(epsilon.alt_F) integral_(-epsilon.alt_F\/tau)^oo dd(x) x^2 e^x/(e^x+1)^2
$
extending the lower limit to $- oo$ we find
$
  C_"el" tilde.eq 1/3 pi^2 D(epsilon.alt_F) k_B^2 T =^"f.e.g" (pi^2 N k_B)/(2) ( T/T_F )
$

For temperatures below both $theta$ and $T_F$ we can then write the heat capacity as
$
  C = C_"lat" + C_"el" tilde.equiv gamma T + A T^3
$
or equivalently
$
  C/T tilde.equiv gamma + A T^2
$

=== Thermal conductivity
We previously showed $K = 1/3 C v cal(l)$. For an electron gas we use $epsilon.alt_F = 1/2 m v_F^2$ and $cal(l) = v_F tau$ giving
$
  K_"el" = pi^2/3 (n k_B^2 T)/(m v_F^2) v_F cal(l) = (pi^2 n k_B^2 T tau)/(3 m)
$

We can compute
$
  K/sigma = ((pi^2 n k_B^2 T tau)/(3 m))/((n e^2 tau)/m) = pi^2/3 (k_B/e)^2 T
$
where we use the $sigma$ and $K$ previously found. We define the Lorenz number
$
  L = K/(sigma T)
$
notably this does not depend on $n$ or $m$. So the ratio of thermal conductivity to electrical conductivity is $prop T$.

= Nearly free electron model
The free electron model is very naive and we miss a lot of detail. We easily get a more sophisticated model by including a periodic crystal potential.

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
  psi_bold(k) (bold(r)) = e^(i bold(k) dot bold(r))
$
these represent free waves with $bold(p) = hbar bold(k)$.

== Bloch's theorem and the central equation
Bloch proved that solutions to the Schrödinger equation for a periodic potential takes the form
$
  psi_bold(k) (bold(r)) = u_bold(k) (bold(r)) e^(i bold(k) dot bold(r))
$
with $u_bold(k) (bold(r)) = u_bold(k) (bold(r)+bold(T))$. So any solution $psi_bold(k)$ takes the form of a free wave with some periodic modulation $u_bold(k)$. These are called Bloch functions. A proof of Bloch's theorem is given below.

We now treat the wave equation for a general potential at general $bold(k)$. Let $U(bold(r))$ be the potential energy of an electron in the lattice, then we know $U(bold(r)+bold(T))=U(bold(r))$. So we can write the potential as a Fourier series
$
  U(bold(r)) = sum_bold(G) U_bold(G) e^(i bold(G) dot bold(r))
$
where the $U_bold(G)$ tend to decrease rapidly.

We now assume $U(bold(r))$ is symmetric about $bold(r)=bold(0)$ and $U_bold(0) = bold(0)$. The Schrödinger equation then becomes
$
  (bold(p)^2/(2 m) + sum_bold(G) U_bold(G) e^(i bold(G)dot bold(r))) psi(bold(r)) = epsilon.alt psi(bold(r))
$
This describes an electron in the potential of ion cores, and in the average potential of other electrons. Note that $U(bold(r)) in RR$ since $U_(-bold(G))=U_bold(G)^*$ and any Bravais lattice has inversion symmetry $U_(-bold(G))=U_bold(G)$.

We can also expand $psi(bold(r))$ as
$
  psi(bold(r)) = sum_bold(k) C(bold(k)) e^(i bold(k)dot bold(r))
$
where $bold(k)$ ensures this satisfies boundary conditions---$k_i = 2 pi n_i\/L$. Not all $bold(k)$ of the given form enter the expansion, since if any $bold(k)$ is in the expansion, then all $bold(k) + bold(G)$ enter aswell. So we label the Bloch functions by the $bold(k)$ in the first Brillouin zone. Note that $hbar bold(k)$ is called the crystal momentum, since $bold(k)$ represents a phase factor picked up by $psi(bold(r))$ after a translation $bold(r)->bold(r)+bold(T)$, and since this is what enters in conservation laws.

Substituting both expansions into the Schrödinger equation gives
$
  sum_bold(k) (hbar^2 k^2)/(2 m) C(bold(k)) e^(i bold(k) dot bold(r)) + sum_(bold(G),bold(k)) U_bold(G) C (bold(k)) e^(i(bold(k)+bold(G))dot bold(r)) = epsilon.alt sum_bold(k) C(bold(k)) e^(i bold(k)dot bold(r) )
$
by matching exponential factors one can see this implies
$
  (epsilon.alt_k^0 - epsilon.alt) C(bold(k)) + sum_bold(G) U_bold(G) C (bold(k)-bold(G)) = 0"  with  " epsilon.alt_k^0 = (hbar^2 k^2)/(2 m)
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
We consider a potential of the form
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
Define
$
  f(k) equiv sum_n C(k - (2 pi n)/a)
$
then
$
  C_k = - ((2 m U_0\/hbar^2) f(k))/(k^2 - (2 m epsilon.alt\/hbar^2))
$
Using
$
  f(k) = f(k - (2 pi n)/a)
$
we can write
$
  C(k - (2 pi n)/a) = - ((2 m U_0)/hbar^2) f(k) [(k-(2 pi n)/a)^2 - (2 m epsilon.alt)/hbar^2]^(-1)
$
and summing over all $n$ we obtain
$
  hbar^2/(2 m U_0) = - sum_n [(k-(2pi n)/a)^2 - ((2 m epsilon.alt)/hbar^2)]^(-1)
$
Using
$
  "ctn" x = sum_n 1/(n pi + x)
$
and defining $(k_epsilon.alt^0)^2 equiv 2 m epsilon.alt\/hbar^2$ we can write the result as
$
  -underbracket((m U_0 a^2)/(2 hbar^2), equiv P) 1/(k_epsilon.alt^0 a) sin k_epsilon.alt^0 a + cos k_epsilon.alt^0 a = cos k a
$
this shows band gaps since $abs(cos k a)<= 1$ so only certain $K tilde epsilon.alt$ are allowed.

== For nearly free electrons
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
  epsilon.alt_plus.minus = epsilon.alt_bold(k) plus.minus abs(U_bold(G)) => E_g = 2U_bold(G)
$
Consider now a $bold(k)$ close to the edge of the first Brillouin zone, with $bold(K) = bold(k) - 1/2 bold(G) tilde.eq 0$. This eventually gives
$
  epsilon.alt_plus.minus &tilde.eq^("Taylor") epsilon.alt_plus.minus^0 + (hbar^2 bold(K)^2)/(2 m) (1 plus.minus (2 epsilon.alt_(1/2 bold(G))^0)/U_bold(G)) " with " epsilon.alt_plus.minus^0 = underbrace(epsilon.alt_(1/2 bold(G))^0, "free energy" #linebreak() "at edge") plus.minus U_bold(G)
$

= Tight-binding model
The tight-binding model is the last electron model we treat in this course. This is in many ways the simplest sophisticated model since we treat the crystal potential as some unknown pertubation and apply pertubation theory.

We start with two separated neutral atoms and consider what happens when they are brought together to form a crystal. Say the atoms have wavefunctions $psi_A$ and $psi_B$. As they are brought together we consider the linear combinations $psi_A plus.minus psi_B$. An electron in the $+$ state will have lower energy than an electron in the $-$ state. This approximation based on the free atom wavefunctions is called linear combination of atomic orbitals or just tight-binding.

We assume the ground state of an electron moving in $U(bold(r))$ of a free atom is $phi(bold(r))$. Taking the influence of one atom on another to be small we approximate the wavefunction of the crystal to be a sum of these
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

To find the first order energy we use the variational principle by computing $braket(bold(k), H, bold(k))$ (for $s$-orbitals)
$
  braket(bold(k), H, bold(k)) &= 1/N sum_(j,m) e^(i bold(k) dot overbrace((r_j-r_m), -bold(rho)_m)) braket(phi_m, H, phi_j) \
  &=^(sum_j = N) sum_m e^(- i bold(k) dot bold(rho)_m) braket(phi(bold(r)-bold(rho)_m-r_j), H, phi(bold(r)-bold(r)_j)) \
  &=^"translational invariance" sum_m e^(- i bold(k) dot bold(rho)_m) braket(phi(bold(r)-bold(rho)_m), H, phi(bold(r)))
$
this is quite nice. To evaluate the matrix elements we assume only nearest neighbors matter
$
  braket(phi(bold(r)), H, phi(bold(r))) = -alpha";  " braket(phi(bold(r)-bold(rho)), H, phi(bold(r))) = - gamma
$
with $bold(rho)$ being nearest neighbor distance. Then
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
and the energies are confined to a band of width $12 gamma$.

#alt[

  Above we used the variational principle
  $
    braket(bold(k), H, bold(k)) >= E_"gs"
  $
  and taken the interaction with nearest neighbors as a pertubative term
  $
    H = underbracket(H_0, "free-electron sum") + V
  $
  We could have written
  $
    braket(bold(k), H, bold(k)) = underbracket(braket(phi, H_0, phi), epsilon.alt_0)+underbracket(braket(phi, V, phi), -beta) + sum_"nearest" e^(-i bold(k) dot bold(rho)_m) [underbracket(braket(phi_bold(rho), H_0, phi), epsilon.alt_0 braket(phi_bold(rho), phi))+underbracket(braket(phi_bold(rho), V, phi), -t)]
  $
  since we take the wavefunctions to be highly localized $braket(phi_bold(rho), phi) tilde 0$ we obtain
  $
    braket(bold(k), H, bold(k)) = underbracket(epsilon.alt_0 - beta, -alpha) - t sum_"nearest" e^(-i bold(k) dot bold(rho)_m) = epsilon.alt_bold(k)
  $
  which is equivalent to the above. This more clearly shows that in the limit $V -> 0$ we just get back
  $epsilon.alt_bold(k) = epsilon.alt_0$.
]
