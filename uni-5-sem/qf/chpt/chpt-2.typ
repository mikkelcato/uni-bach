//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Quantum Dynamics
== Time evolution operator
We want to describe how a system starting in $ket(alpha)$ at $t_0$ evolves in time. We define the time-evolution operator $U (t,t_0)$ by
$
  ket(alpha\,t_0\;t) = U (t,t_0) underbracket(ket(alpha\,t_0), equiv ket(alpha))
$
At $t_0$ we can write
$
  ket(alpha) = sum_a' c_(a') (t_0) ket(a')
$
and similarly for $t$
$
  ket(alpha\,t_0\;t) = sum_a' c_(a') (t) ket(a')
$
We must have
$
  sum_a' abs(c_(a') (t_0))^2 = sum_a' abs(c_(a') (t))^2
$
meaning a state remains normalized in time. This is true when $U^dagger U = bb(1)$. We also want
$
  U (t_2, t_0) &= U (t_2, t_1) U (t_1, t_0) #h(2em) (t_2 > t_1 > t_0) \
  lim_(dd(t)->0) U(t_0+dd(t),t_0) &= bb(1)
$
We claim
$
  U (t_0 + dd(t),t_0) = 1 - i Omega dd(t)
$
with $Omega$ being Hermitian satisfies the above. This is simple to check.

We write
$
  U (t_0 + dd(t),t_0) = 1 - (i H dd(t))/hbar
$
with $H$ being the Hamiltonian.

== The Schrödinger equation
Consider
$
  U (t + dd(t),t_0) & = U (t + dd(t),t) U (t,t_0) \
                    & = (1 - (i H dd(t))/hbar) U (t,t_0)
$
which can be written as
$
  i hbar pdv(, t) U (t,t_0) = H U(t,t_0)
$
this is the Schrödinger equation for $U (t,t_0)$! We can immediately get the equation for $ket(alpha\, t_0\;t)$ by
$
  i hbar pdv(, t) U (t,t_0) ket(alpha) & = H U (t,t_0) ket(alpha) \
    i hbar pdv(, t) ket(alpha\,t_0\;t) & = H ket(alpha\,t_0\;t)
$
if we know $U(t,t_0)$ this is overkill. Assuming $H$ is time-independent we find
$
  U (t,t_0) = exp[(-i H (t-t_0))/hbar]
$
which can be proved in many ways.

We want the action on $ket(alpha)$. Assuming $[A,H] = 0$ then
$
  H ket(a') = E_a' ket(a')
$
where we would call $ket(a')$ energy eigenkets. Letting $t_0 = 0$ we then have
$
  exp((-i H t)/hbar) &= sum_(a') sum_(a'') ketbra(a'') exp((-i H t)/hbar) ketbra(a') \
  &= sum_a' ket(a') exp((-i E_a' t)/hbar) bra(a')
$
Assuming we know the expansion
$
  ket(alpha) = sum_a' ket(a') braket(a', alpha)
$
then
$
  ket(alpha\,t_0 = 0\;t) = exp((-i H t)/hbar) ket(alpha) = sum_a' ket(a') underbracket(braket(a', alpha), c_a' (0)) exp((-i E_a' t)/hbar)
$
meaning we have
$
  c_a' (t) = c_a' (0) exp((-i E_a' t)/hbar)
$
Let $ket(alpha) = ket(a')$ then
$
  ket(a'\,t_0 = 0\;t) = ket(a') exp((-i E_a' t)/hbar)
$
So if the system is in a simultaneous eigenstate of $A$ and $H$ then it will stay in it. We call observables $A$ with $[A, H] = 0$ _constants of motion_.



Assume $ket(alpha) = ket(a')$ and $[A,H] = 0$. Then
$
  expval(B) & = bra(a') U^dagger (t,0) B U (t,0) ket(a') \
            & = braket(a', B, a')
$
which is independent of time! So we call energy eigenstates _stationary_.

As an example consider
$
  H = - e/(m_e c) bold(S) dot bold(B)
$
We assume $bold(B) = B hat(bold(z))$
$
  H = - (e B)/(m_e c) S_z equiv omega S_z
$
meaning $[H,S_z] = 0$ so the $S_z$ eigenstates are energy eigenstates. Then
$
  H ket(plus.minus) = E_plus.minus ket(plus.minus) = plus.minus (hbar omega)/2 ket(plus.minus)
$
The time-evolution operator is
$
  U (t,0) = exp((- i omega S_z t)/hbar)
$
take an arbitrary initial state
$
  ket(alpha) = c_+ ket(+) + c_- ket(-)
$
then
$
  ket(alpha\,t_0=0\;t) = c_+ exp((-i omega t)/2) ket(+) + c_- exp((i omega t)/2) ket(-)
$

== Heisenberg picture
The above takes the state as dynamic. This is the Schrödinger picture in the Heisenberg picture the observables are dynamic.

Consider
$
  braket(beta, X, alpha) -> braket(beta, U^dagger X U, alpha)
$
with $U$ unitary. Two approaches follow
$
  ket(alpha) & -> U ket(alpha) --> "Schrödinger" \
           X & -> U^dagger X U --> "Heisenberg"
$
We fix all states $ket(alpha)$ at $t_0 = 0$. Then
$
  U (t, 0) eq U (t) = exp((-i H t)/hbar)
$
and we define
$
  A^((H)) (t) equiv U^dagger (t) A^((S)) U(t)
$
at $t=0$ the pictures coincide $A^((H)) (0) = A^((S))$. We treat the state as _frozen_
$
  ket(alpha\,t_0=0\;t)_H = ket(alpha)
$
unlike before where
$
  ket(alpha\,t_0=0\;t)_S = U (t) ket(alpha)
$
Consider the expectation value
$
  expval(A^((S)), alpha\,t_0=0\;t)_S & = expval(U^dagger A^((S)) U, alpha\,t_0=0) \
                                     & = expval(A^((H)) (t), alpha\,t_0 = 0)_H
$
it does not care! Which is important since this is what we measure.

Assume $A^((S))$ does not depend on time. Then
$
  dv(A^((H)), t) & = pdv(U^dagger, t) A^((S)) U + U^dagger A^((S)) pdv(U, t) \
                 & = 1/(i hbar) [A^((H)), U^dagger H U]
$
with $[H,U] = 0$ we find
$
  dv(A^((H)), t) = 1/(i hbar) [A^((H)), H]
$
this is the Heisenberg equation of motion!

We have
$
  A^((H)) (t) = U^dagger A(0) U
$
so at $t = 0$
$
  U^dagger A(0) U U^dagger ket(a') & = a' U^dagger ket(a') \
        A^((H)) (U^dagger ket(a')) & = a' (U^dagger ket(a'))
$
meaning ${cal(U)^dagger ket(a')}$ are the Heisenberg basekets. They evolve as
$
  ket(a'\,t)_H = U^dagger ket(a')
$
and satisfy
$
  i hbar pdv(, t) ket(a'\,t)_H = - H ket(a'\,t)_H
$
so basekets evolve oppositely!

== Ehrenfest's theorem
We assume the Hamiltonian takes the same form as in classical mechanics with operators replacing $x_i$ and $p_i$. We typically need
$
  [x_i, F(bold(p))] = i hbar pdv(F, p_i)";  " [p_i, G(bold(x))] = - i hbar pdv(G, x_i)
$
where $F$ and $G$ are power series in $p_j$ and $x_j$ respectively.

Consider a free particle. The Hamiltonian is
$
  H = (p_x^2 + p_y^2 + p_z^2)/(2 m)
$
We have
$
  dv(p_i, t) = 1/(i hbar) [p_i, H] = 0
$
so $p_i$ is a constant of motion. Also
$
  dv(x_i, t) = 1/(i hbar) [x_i, H] =^"by identity" 1/(2 m) pdv(, p_i) sum_j^3 p_j^2 = p_i/m =^"constant of motion" (p_i (0))/m
$
solving for $x_i (t)$ we find
$
  x_i (t) = x_i (0) + (p_i (0))/m t
$

$
  [x_i (t), x_i (0)] = [(p_i (0) t)/m, x_i (0)] = - (i hbar t)/m
$
We add a potential $V(bold(x))$
$
  H = bold(p)^2/(2m) + V(bold(x))
$
then
$
  dv(p_i, t)= 1/(i hbar) [p_i, V(bold(x))] =^"identity" - pdv(, x_i) V(bold(x))
$
and
$
  dv(x_i, t, 2) = 1/(i hbar)[dv(x_i, t),H] = 1/(i hbar) [p_i/m,H] = 1/m dv(p_i, t)
$
Combing these we obtain
$
  m dv(bold(x), t, 2) = - grad V(bold(x))
$
which is the quantum mechanical Newton's second law! Taking the expectation value with resepct to a Heisenberg state
$
  m dv(, t, 2) expval(bold(x)) = dv(expval(bold(p)), t) = - expval(grad V(bold(x)))
$
which is Ehrenfest's theorem! This holds in both pictures.

== The simple harmonic oscillator
We summarize the main results.

The Hamiltonian is
$
  H = p^2/(2m) + (m omega^2 x^2)/2
$
We define the ladder operators
$
         a & = sqrt((m omega)/(2 hbar)) (x + (i p)/(m omega)) \
  a^dagger & = sqrt((m omega)/(2 hbar)) (x - (i p)/(m omega))
$
they satisfy
$
  [a,a^dagger] = 1
$
We define the number operator $N = a^dagger a$ giving
$
  H = hbar omega ( N + 1/2)
$
so $[H,N]=0$. We denote the eigenket of $N$ by $ket(n)$ with $N ket(n) = n ket(n)$. Then
$
  H ket(n) = hbar omega (n + 1/2) ket(n) equiv E_n ket(n)
$
We compute
$
         [N,a] & = - a \
  [N,a^dagger] & = a^dagger
$
from which
$
         a ket(n) & = sqrt(n) ket(n-1) \
  a^dagger ket(n) & = sqrt(n+1) ket(n+1)
$
follow.

We define the ground state by $a ket(0) = 0$ and all higher states can be derived by applying $a^dagger$
$
  ket(n) & = ((a^dagger)^n/sqrt(n!)) ket(0)
$
see Sakurai for more.

== The wave equation
We want the time-evolution of the wavefunction $psi(bold(x)', t)$.

Let the Hamiltonian be
$
  H = bold(p)^2/(2 m) + V(bold(x))
$
and take $V(bold(x)')$ as local $braket(bold(x)'', V(bold(x)), bold(x)') = V(bold(x)') delta^3 (bold(x)'-bold(x)'')$.

Then acting with $bra(bold(x)') dot (dots)$ on the Schrödinger equation gives
$
  i hbar pdv(, t) braket(bold(x)', alpha\,t_0\;t) = braket(bold(x)', H, alpha\,t_0\;t)
$
which can be written as
$
  i hbar pdv(, t) underbracket(braket(bold(x)', alpha\,t_0\;t), psi(bold(x)', t)) = - hbar^2/(2m) nabla'^2 braket(bold(x)', alpha\,t_0\;t) + V(bold(x)') braket(bold(x)', alpha\,t_0\;t)
$
this should be very familiar!

Take $ket(alpha)$ to be an energy eigenstate. Then
$
  braket(bold(x)', a'\,t_0\;t) = braket(bold(x)', a') exp((- i E_a' t)/hbar)
$
with $ket(a')$ an eigenstate of $A$ with $[A,H]=0$. We obtain
$
  - hbar^2/(2m) nabla'^2 underbracket(braket(bold(x)', a'), "energy" #linebreak() "eigenfunction") + V(bold(x)') braket(bold(x)', a') = E_a' braket(bold(x)', a')
$
We let $A$ be the function of $bold(x)$ and $bold(p)$ that coincides with $H$. Then we can write
$
  - hbar^2/(2m) nabla'^2 u_E (bold(x)') + V(bold(x)') u_E (bold(x)') = E u_E (bold(x)')
$
which is the time-independent Schrödinger wave equation!

To solve it we take $ E < lim_(abs(bold(x)') arrow oo) V(bold(x)') $
meaning
$
  u_E (bold(x)') arrow 0 "as" abs(bold(x)') arrow oo
$
so the particle must be bound. This leads to a discrete spectrum of energies. Likewise if the condition is not satisfied we have scattering states with a continuum of energies.

== Free particle
Consider $V(bold(x)) = 0$.

The time-independent Schrödinger wave equation becomes
$
  nabla^2 u_E (bold(x)) = - (2 m E)/hbar^2 u_E (bold(x))
$
We define
$
  bold(k)^2 = k_x^2 +k_y^2 + k_z^2 equiv (2 m E)/hbar^2 = bold(p)^2/hbar^2
$
and make the ansatz $u_E (bold(x)) = u_x (x) u_y (y) u_z (z)$. We obtain
$
  (1/u_x dv(u_x, x, 2) + k_x^2) + (1/u_y dv(u_y, y, 2) + k_y^2) + (1/u_z dv(u_z, z, 2) + k_z^2) = 0
$
which has solutions
$
  u_w (w) = c_w e^(i k_w w)" for "w={x,y,z}
$
meaning
$
  u_E (bold(x)) = c_x c_y c_z e^(i k_x x + i k_y z + i k_z z) = C e^(i bold(k) dot bold(x))
$
This is not normalizable usually. We use big-box normalization and say all space is within a cube of side $L$ and assume periodic boundaries. This gives
$
  u_x (x + L) = u_x (x) => k_x L = 2 pi n_x => k_x = (2 pi)/L n_x
$
and similarly for $y$ and $z$. We normalize to find
$
  1 = integral_0^L dd(x) integral_0^L dd(y) integral_0^L dd(z) u_E^* (bold(x)) u_E (bold(x)) = L^3 abs(C)^2
$
meaning $C = L^(-3\/2)$ so
$
  u_E (bold(x)) = 1/(L^(3\/2)) e^(i bold(k) dot bold(x))
$
with energies
$
  E = bold(p)^2/(2m) = hbar^2/(2m) ((2pi)/L)^2 (n_x^2 + n_y^2 + n_z^2)
$

== Linear potential
Consider $V(x) = k abs(x)$. This potential has a classical turning point at some $x = a$ with $E = k a$.

The time-independent Schrödinger wave equation becomes
$
  - hbar^2/(2m) dv(u_E, x, 2) + k abs(x) u_E (x) = E u_E (x)
$
due to $V(-x) = V(x)$ we care about $x >= 0$.

We find two types of solutions $u_E (-x) = plus.minus u_E (x)$. Assuming $u_E (-x) = - u_E (x)$ (odd) we require $u_E (0) = 0$. While for $u_E (-x) = u_E (x)$ (even) we require $u'_E (0) = 0$ since $u_E (epsilon) - u_E (-epsilon)$ vanishes for $epsilon arrow 0$.

We define
$
  x_0 & = ((hbar^2)/(m k))^(1\/3) \
  E_0 & = k x_0 = ((hbar^2 k^2)/m)^(1\/3)
$
giving dimensionless $y equiv x\/x_0$ and $epsilon equiv E\/E_0$. We obtain
$
  dv(u_E, y, 2) - 2 (y-epsilon) u_E (y) = 0 "for" y>=0
$
We introduce $z equiv 2^(1\/3) (y-epsilon)$ giving
$
  dv(u_E, z, 2) - z u_E (z) = 0
$
which is the Airy equation. The solution is the Airy function $"Ai"(z)$. The boundary conditions become zeroes for $"Ai"'(z)$ and $"Ai" (z)$ with $z = - 2^(1\/3) epsilon$. These determine the quantized energies.

Though the linear potential seems non-physical it appears in many places and is always present due to the gravitational force.

== The WKB approximation (briefly)
The WKB approximation uses the linear potential to approximate any potential.

We can write the time-independent Shcrödinger wave equation as
$
  dv(u_E, x, 2) + (2 m)/hbar^2 (E - V(x)) u_E (x) = 0
$
We define
$
                     k(x) & equiv [(2 m)/hbar^2 (E-V(x))]^(1\/2) "for" E > V(x) \
  k(x) equiv - i kappa(x) & eq - i [(2 m)/hbar^2 (V(x)-E)]^(1\/2) "for" E < V(x)
$
then
$
  dv(u_E, x, 2) + k^2(x) u_E (x) = 0
$
Assuming $V(x)$ varies slowly we make the ansatz
$
  u_E (x) equiv e^(i W(x) hbar^(-1))
$
giving an equation for $W$
$
  i hbar dv(W, x, 2) - (dv(W, x))^2 + hbar^2 k^2(x) = 0
$
_Varying slowly_ means
$
  hbar abs(dv(W, x, 2)) << abs(dv(W, x))^2
$
so the first term is small giving a zeroth-order approximation for $W(x)$
$
  W'_0 (x) = plus.minus hbar k(x)
$
A first-order approximation is then obtained
$
  (dv(W_1, x))^2 = hbar^2 k(x)^2 plus.minus i hbar^2 k' (x)
$
Taking $W tilde W_1$ we find
$
  W(x) approx W_1 (x) &= plus.minus hbar integral^x dd(x)' [k^2 (x') plus.minus i k' (x')]^(1\/2) \
  &approx^"binom" plus.minus hbar integral^x dd(x)' k(x') [1 plus.minus i/2 (k'(x'))/(k^2 (x'))] \
  &= plus.minus hbar integral^x dd(x') k(x') + i/2 hbar ln k(x)
$
Then
$
  u_E (x) approx e^(i W(x) hbar^(-1)) = 1/sqrt(k(x)) exp(plus.minus i integral^x dd(x)' k(x'))
$
this gives solutions for $E > V$ and $E < V$. We will not cover the joining procedure.

Consider a potential well with turning points $x_1$ and $x_2$ creating three regions. In the middle region the wave function behaves like our approximation with the first $k(x)$ and in the two outer regions with the second $k(x)$. Around the turning points the solutions are given by Airy functions since we assume a linear approximation in those regions. This leads to a consistency check
$
  integral_(x_1)^(x_2) dd(x) sqrt(2 m (E-V(x))) = (n+1/2) pi hbar
$
giving approximate expressions for the energies.

Consider a bouncing ball with
$
  V = cases(m g x "  for" x > 0, oo "     for" x < 0)
$
where $x$ is the height from the surface. We could use $x_1 = 0$ and $x_2 = E\/m g$. These are the classical turning points. But our wave function leaks into the $x < x_1$ region even though we require it vanishes. For this reason we use the odd solution which vanishes at $x = 0$. So we have
$
  V(x) = m g abs(x)
$
with turning points $x_1 = - E\/m g$ and $x_2 = E\/m g$. Then
$
  integral_(- E\/m g)^(E\/m g) dd(x) sqrt(2m (E-m g abs(x))) &= (n_"odd" + 1/2) pi hbar \
  integral_0^(E\/m g) dd(x) sqrt(2 m (E- m g x)) &= (n - 1/4) pi hbar
$
giving
$
  E_n = {[3(n - 1\/4) pi]^(2\/3)/2} (m g^2 hbar^2)^(1\/3)
$

The WKB limit is equivalent to
$
  lambda = hbar/sqrt(2 m [E-V(x)]) << (2 (E-V(x)))/abs(dd(V)\/dd(x))
$

== The propagator
Consider
$
  ket(alpha\,t_0\;t) = sum_a' ket(a') braket(a', alpha) exp[(-i E_a' (t-t_0))/hbar]
$
or in terms of the wavefunction
$
  psi(bold(x)'', t) &= sum_a' braket(a', alpha) braket(bold(x)'', a') exp[(-i E_a' (t-t_0))/hbar] \
  &= integral dd(x', 3) K (bold(x)'', t; bold(x)',t_0) underbracket(psi(bold(x)', t_0), braket(bold(x)', alpha))
$
where we define the propagator $K$ by
$
  K(bold(x)'',t;bold(x)',t_0) equiv sum_a' braket(bold(x)'', a') braket(a', bold(x)') exp[(-i E_a' (t-t_0))/hbar]
$
This guy has many properties including
$
  lim_(t arrow t_0) K(bold(x)'',t;bold(x)',t_0) = delta^3 (bold(x)''-bold(x)')
$
and it is a Green's function for the Schrödinger wave equation
$
  [- hbar^2/(2m) nabla''^2 + V(bold(x)'') - i hbar pdv(, t)] K(bold(x)'',t;bold(x)',t_0) = - i hbar delta^3 (bold(x)''-bold(x)') delta(t-t_0)
$
given $K = 0$ for $t < t_0$.

Let $t_0 = 0$ and $bold(x)''=bold(x)'$. Then
$
  G(t) &equiv integral dd(x', 3) K(bold(x)',t;bold(x)',0) \
  &= integral dd(x', 3) sum_a' abs(braket(bold(x)', a'))^2 exp((-i E_a' t)/hbar) \
  &= sum_a' exp((-i E_a' t)/hbar)
$
this looks like a partition function! And by a Wick rotation we find
$
  Z = sum_a' e^(- beta E_a') " with " beta=(i t)/hbar
$
we can also find
$
  tilde(G) (E) = sum_a' 1/(E-E_a')
$
where $tilde(G)$ is the Laplace-Fourier transform of $G$.

We can write the propagator as
$
  K(bold(x)'',t;bold(x)',t_0) &= sum_a' braket(bold(x)'', exp((-i H t)/hbar), a') braket(a', exp((i H t_0)/hbar), bold(x)') \
  &= braket(bold(x)'', U (t,t_0), bold(x)') \
  &= braket(bold(x)''\,t, bold(x)'\,t_0)
$
with $ket(bold(x)'\,t_0)$ and $bra(bold(x)''\,t)$ in the Heisenberg picture. So the propagator is the transition amplitude for a particle to go from $(bold(x)',t_0)$ to $(bold(x)'',t)$.

We can also split some interval $(t',t''')$ into $(t',t'')$ and $(t'',t''')$. Then
$
  braket(bold(x)'''\,t''', bold(x)'\,t') = integral dd(x'', 3) braket(bold(x)'''\,t''', bold(x)''\,t'') braket(bold(x)''\,t'', bold(x)'\,t')
$
since
$
  integral dd(x'', 3) braket(bold(x)''\,t'', bold(x)''\,t'') = 1
$
this can of course be done for more subdivisions.

== The path integral
Consider the transition amplitude for $(x_1,t_1) arrow (x_N, t_N)$ with $N-1$ equal intervals of length
$
  Delta t = (t_N-t_1)/(N-1)
$
By the composition property
$
  braket(x_N\,t_N, x_1\,t_1) &= integral dd(x_(N-1)) dots integral dd(x_2) braket(x_N\,t_N, x_"N-1"\,t_"N_1") dots braket(x_2\,t_2, x_1\,t_1)
$
this means that for each $t_(n-1) arrow t_n$ we integrate over all paths $x_2,x_3,dots, x_(N-1)$.

We define
$
  S(n,n-1) equiv integral_(t_(n-1))^t_n dd(t) L_"classical" (x,dot(x))
$
Consider some small segment $(x_(n-1),t_(n-1)) -> (x_n,t_n)$ and associate $exp^(i S(n,n-1) hbar^(-1))$ to it. Along any path we then get the contribution
$
  product_(n=2)^N exp[(i S(n,n-1))/hbar] = exp[i/hbar sum_(n=2)^N S(n,n-1)] = exp[(i S(N,1))/hbar]
$
to $braket(x_N\,t_N, x_1\,t_1)$. We obtain
$
  braket(x_N\,t_N, x_1\,t_1) tilde sum_"all paths" exp[(i S(N,1))/hbar]
$
in the limit $hbar arrow 0$ the exponential will oscillate a bunch leading to many paths cancelling. However given a path satisfies $delta S(N,1) = 0$ i.e. it would be the classical path, then slight deformations do not change $exp(i S\/hbar)$. So in the $hbar -> 0$ limit this leads to constructive interference with the classical path being singled out and all others vanishing.

To make this precise we write
$
  braket(x_n\,t_n, x_(n-1)\,t_(n-1)) = 1/w(Delta t) exp[(i S(n,n-1))/hbar]
$
where we assume $t_n-t_(n-1)$ is an infinitesimal. We now determine the weight $w(Delta t)$. Due to $Delta t -> 0$ we treat the path between $t_(n-1) -> t_n$ as a straight line
$
  S(n,n-1) & = integral_(t_(n-1))^(t_n) dd(t) [(m dot(x)^2)/2 - V(x)] \
           & = Delta t { m/2 [(x_n-x_(n-1))/(Delta t)]^2 - V((x_n +x_(n-1))/2)}
$
$w(Delta t)$ does not depend on $V$ so we let $V = 0$ giving
$
  braket(x_n\,t_n, x_(n-1)\,t_(n-1)) = 1/w(Delta t) exp[(i m(x_n-x_(n-1))^2)/(2 hbar Delta t)]
$
in the limit $t_n = t_(n-1)$ this should be $delta(x_n - x_(n-1))$ giving
$
  1/w(Delta t) = sqrt(m/(2 pi i hbar Delta t))
$
So we find
$
  braket(x_n\,t_n, x_(n-1)\,t_(n-1)) = sqrt(m/(2 pi i hbar Delta t)) exp[(i S(n,n-1))/hbar]
$
meaning
$
  braket(x_N\,t_N, x_(N-1)\,t_(N-1)) &= lim_(N arrow oo) (m/(2 pi i hbar Delta t))^((N-1)\/2) integral dd(x_(N-1)) dots integral dd(x_2) product_(n=2)^N exp[(i S(n,n-1))/hbar]
$
We define the operator
$
  integral_(x_1)^(x_N) D [x(t)] equiv lim_(N -> oo) (m/(2 pi i hbar Delta t))^((N-1)\/2) integral dd(x_(N-1)) dots integral dd(x_2)
$
giving
$
  braket(x_N\,t_N, x_(N-1)","t_(N-1)) &= integral_(x_1)^(x_N) D [x(t)] exp[i integral_(t_1)^(t_N) dd(t) (L_"classical" (x,dot(x)))/hbar] \
  &= integral_(x_1)^(x_N) D [x(t)] exp[(i S(N,1))/hbar]
$
this is the Feynman path integral.

/*
== Alternate derivation of the path integral
We define
$
  Delta t = (t_f - t_i)/N
$
where $t_i = t_0$ and $t_f = t$, similarly $x' = x_i$ and $x''=x_f$. We can write
$
  braket(x_f","t_f, x_i","t_i) &= braket(x_f, exp[(-i H (t_f-t_i))/hbar], x_i) \
  &= braket(x_f, exp[(-i H(t_f-t_i))/(N hbar)]^N, x_i) \
  &= braket(x_f, exp[(-i H Delta t)/hbar]^N, x_i) \
  &=^"insert identity" braket(x_N, exp[(-i H Delta t)/hbar] integral dd(x_(N-1)), x_(N-1))bra(x_(N-1)) dots \ & dots integral dd(x_1) ket(x_1) braket(x_1, exp[(-i H Delta t)/hbar], x_0) \
  &= integral product_(i=1)^(N-1) dd(x_i) product_(j=1)^N braket(x_j, exp[(- i H Delta t)/hbar], x_(j-1)) \
  &=^"insert identity" integral product_(i=1)^(N-1) dd(x_i) product_(j=1)^N dd(p_i) product_(k=1)^N braket(x_k, p_k) braket(p_k, exp[(-i H Delta t)/hbar], x_(k-1))
$
this is exact. We assume
$
  H = p^2/(2m) + V(x)
$
we rewrite the last term
$
  braket(p_j, exp[(-i H Delta t)/hbar], x_(j-1)) &= braket(p_j, 1-i/hbar H Delta t + dots, x_(j-1)) \
  &= braket(p_j, 1-i/hbar [p^2/(2 m) + V(x)] Delta t, x_(j-1)) \
  &= [1-i/hbar (p_j^2/(2 m) + V(x_(j-1))) Delta t ] 1/(sqrt(2 pi hbar)) exp[(- i p_j x_(j-1))/hbar] \
  &= 1/sqrt(2 pi hbar) exp[- i/hbar (p_j x_(j-1) + {p_j^2/(2m) + V(x_(j-1))} Delta t)] \
  &= 1/sqrt(2 pi hbar) exp[- i/hbar (p_j x_(j-1) + H(x_(j-1),p_j) Delta t)]
$
so
$
  braket(x_j, p_j) braket(p_j, exp[(-i H Delta t)/hbar], x_(j-1)) &= 1/(2 pi hbar) exp[i/hbar p_j (x_j-x_(j-1)) - i/hbar H(x_(j-1),p_j) Delta t]
$
and
$
  product_(j=1)^N dots = 1/(2 pi hbar)^N exp[i/hbar Delta t sum_(j=1)^N {p_j (x_j - x_(j-1))/(Delta t) - H(x_(j-1),p_j)}]
$
in the limit $N arrow oo$ we find
$
  braket(x_f","t_f, x_i","t_i) &= lim_(N -> oo) integral product_(i=1)^(N-1) dd(x_i) product_(j=1)^N dd(p_j)/(2pi hbar) exp[i/hbar Delta t sum_(k=1)^N {p_k (x_k-x_(k-1))/Delta t - H(x_(k-1),p_k)}] \
  &equiv integral dd(x, p, d: D) exp[i/hbar integral_(t_i)^(t_f) dd(t) {p dot(x) - H(x,p)}] \
  &= integral dd(x, p, d: D) exp[i/hbar S_H (x,p)]
$
with the Hamiltonian action being
$
  S_H (x,p) = p dot(x) - H(x,p)
$
and
$
  integral dd(x, d: D) = lim_(N -> oo) product_(i=1)^(N-1) dd(x_i)",  " integral dd(p, d: D) = lim_(N -> oo) product_(j=1)^N dd(p_j)/(2 pi hbar)
$
this is the path integral in phase space---we perform the integration over $p$ to get the path integral in configuration space.
$
  product_(j=1)^N & integral dd(p_j)/(2pi hbar) exp[i/hbar Delta t sum_(k=1)^N {p_k (x_k-x_(k-1))/Delta t - H(x_(k-1),p_k)}] \
  &= product_(j=1)^N integral dd(p_j)/(2 pi hbar) exp[i/hbar Delta t sum_(k=1)^N {-1/(2 m) (p_k - m (x_k - x_(k-1))/(Delta t))^2 + m/2 ((x_k-x_(k-1))/(Delta t))^2 - V(x_(k-1))}] \
  &= exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k-x_(k-1))/(Delta t))^2 - V(x_(k-1))}] product_(j=1)^N integral dd(p_j)/(2 pi hbar) exp[i/hbar Delta t sum_(k=1)^N {-1/(2m) (p_k - m (x_k-x_(k-1))/(Delta t))^2}]
$
letting
$
  p'_k = p_k - m (x_k-x_(k-1))/(Delta t)
$
we find
$
  product_(j=1)^N & integral dd(p_j)/(2pi hbar) exp[i/hbar Delta t sum_(k=1)^N {p_k (x_k-x_(k-1))/Delta t - H(x_(k-1),p_k)}] \
  &= exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k-x_(k-1))/(Delta t))^2 - V(x_(k-1))}] product_(j=1)^N (integral dd(p'_j)/(2 pi hbar)) exp[i/hbar Delta t sum_(k=1)^N {-1/(2m)p'_j^2}] \
  &= exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k - x_(k-1))/(Delta t))^2 - V(x_(k-1))}] (integral dd(p')/(2 pi hbar) exp[- 1/2 {(i Delta t)/(m hbar)} p'^2])^N \
  &= exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k - x_(k-1))/(Delta t))^2 - V(x_(k-1))}] (1/(2 pi hbar) sqrt((2 pi m hbar)/(i Delta t)))^N \
  &= (m/(2 pi i hbar Delta t))^(N\/2) exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k - x_(k-1))/(Delta t))^2 - V(x_(k-1))}]
$
so we find
$
  braket(x_f","t_f, x_i","t_i) &= lim_(N -> oo) integral (product_(i=1)^(N-1) dd(x_i)) (m/(2 pi i hbar Delta t))^(N\/2) exp[i/hbar Delta t sum_(k=1)^N {m/2 ((x_k - x_(k-1))/(Delta t))^2 - V(x_(k-1))}] \
  & equiv integral dd(x, d: D) exp[i/hbar integral_(t_i)^(t_f) dd(t) {1/2 m dot(x)^2 - V(x)}] \
  &= integral dd(x, d: D) exp[i/hbar S(x)]
$
with the Lagrangian action being
$
  S(x) = integral_(t_i)^(t_f) L(x,dot(x)) dd(t) = integral_(t_i)^(t_f) dd(t) (1/2 m dot(x)^2 - V(x))
$
and
$
  integral dd(x, d: D) = lim_(N -> oo) c(N) integral product_(i=1)^(N-1) dd(x_i)",  " c(N) = (m/(2 pi i hbar Delta t))^(N\/2)
$
this is the expression the book found.
*/

== Gauge transformations
Consider adding a constant $V_0$ to some $V(bold(x))$.

We want to relate $overline(ket(alpha\,t_0\;t))$ for $overline(V)(bold(x)) = V(bold(x))+V_0$ to $ket(alpha\,t_0\;t)$. Assuming they coincide at $t = t_0$ we find
$
  overline(ket(alpha\,t_0\;t)) &= exp[-i (bold(p)^2/(2m) + V(x) + V_0) (t-t_0)/hbar] ket(alpha) \
  &= exp[(-i V_0 (t-t_0))/hbar] ket(alpha\,t_0\;t)
$
so $V_0$ adds a phase! This is an example of a gauge transformation.

For $V(bold(x))+V_0 (t)$ we find
$
  ket(alpha\,t_0\;t) -> exp[- i integral_(t_0)^(t) dd(t') (V_0 (t'))/hbar] ket(alpha\,t_0\;t)
$

Recall from electromagnetism
$
  bold(E) = - grad phi",  " bold(B) = curl bold(A)
$
and
$
  H = 1/(2 m) (bold(p) - (e bold(A))/c)^2 + e phi
$
We say $phi$ and $bold(A)$ are functions of the operator $bold(x)$ in quantum mechanics. Meaning $bold(p)$ and $bold(A)$ do not commute and we write
$
  (bold(p) - (e bold(A))/c)^2 -> p^2 - (e/c) (bold(p) dot bold(A) + bold(A) dot bold(p)) + (e/c)^2 bold(A)^2
$
in the Heisenberg picture
$
  dv(x_i, t) = [x_i,H]/(i hbar) = (p_i - e A_i\/c)/m
$
We call $bold(p)$ the canonical momentum and the new quantity the kinematical momentum
$
  bold(Pi) equiv m dv(bold(x), t) = bold(p) - (e bold(A))/c
$
which satisfies
$
  [Pi_i,Pi_j] = (i hbar e)/c epsilon_(i j k) B_k
$
We can also obtain the quantum mechanical Lorentz law
$
  m dv(bold(x), t, 2) = dv(bold(Pi), t) = e [bold(E) + 1/(2 c) (dv(bold(x), t) times bold(B)-bold(B) times dv(bold(x), t))]
$
Sandwiching $H$ between $bra(bold(x)')$ and $ket(alpha\,t_0\;t)$ gives
$
  1/(2m) [- i hbar nabla' - (e bold(A) (bold(x)'))/c ] [- i hbar nabla' - (e bold(A) (bold(x)'))/c] braket(bold(x)', alpha","t_0";"t) + e phi (bold(x)') braket(bold(x)', alpha","t_0";"t) = i hbar pdv(, t) braket(bold(x)', alpha","t_0";"t)
$
which we call minimal coupling and amounts to replacing $bold(p) -> bold(Pi)$ in the Schrödinger equation!

The usual gauge transformation is given by
$
  phi -> phi - 1/c pdv(Lambda, t)",  " bold(A) -> bold(A) + grad Lambda
$
which leaves $bold(E)$ and $bold(B)$ invariant given
$
  bold(E) = - grad phi - 1/c pdv(bold(A), t)",  " bold(B) = curl bold(A)
$
we only worry about time-independent fields. Classically $bold(p)$ is not gauge-invariant while $bold(Pi)$ is. We demand that expectation values behave similarly. So $expval(bold(x))$ and $expval(bold(Pi))$ should be gauge-invariant while $expval(bold(p))$ is expected to change.

Let $tilde(bold(A)) = bold(A) + grad Lambda$ then
$
  braket(alpha, bold(x), alpha) = braket(tilde(alpha), bold(x), tilde(alpha))
$
and
$
  expval((bold(p)-(e bold(A))/c), alpha) = expval((bold(p)-(e tilde(bold(A)))/c), tilde(alpha))
$
the norm should also be preserved $braket(alpha) = braket(tilde(alpha))$. We define the gauge operator $cal(G)$ by
$
  ket(tilde(alpha)) = cal(G) ket(alpha)
$
meaning it must satisfy
$
  cal(G)^dagger bold(x) cal(G) &= bold(x) \
  cal(G)^dagger (bold(p) - (e bold(A))/c - (e grad Lambda)/c) cal(G) &= bold(p) - (e bold(A))/c
$
We claim
$
  cal(G) = exp[(i e Lambda(bold(x)))/(hbar c)]
$
to see this satisfies the above consider
$
  exp[(-i e Lambda)/(hbar c)] bold(p) exp[(i e Lambda)/(hbar c)] &= exp[(-i e Lambda)/(hbar c)] (bold(p), exp[(i e Lambda)/(hbar c)]) + bold(p) \
  &= - exp[(-i e Lambda)/(hbar c)] i hbar grad (exp[(i e Lambda)/(hbar c)]) + bold(p) \
  &= bold(p) + (e grad Lambda)/c
$
So any gauge transformation corresponds to adding a phase $cal(G) = exp[i e Lambda (bold(x))\/hbar c]$! We find
$
  tilde(psi) (bold(x)',t) = exp[(i e Lambda(bold(x)'))/(hbar c)] psi(bold(x)', t)
$
as we wanted.
