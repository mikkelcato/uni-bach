//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Quantum Dynamics
== Time evolution
We want to describe how some system starting in $ket(alpha)$ at $t_0$ evolves over time. We want
$
  ket(alpha","t_0";"t)"  "(t > t_0)
$
since time is continuous we expect $lim_(t arrow t_0) ket(alpha","t_0";"t) = ket(alpha)$, and we denote $ket(alpha","t_0";"t_0) = ket(alpha","t_0) = ket(alpha)$. We are interested in how $ket(alpha)$ changes under $t_0 arrow t$.

=== The operator
We define a time-evolution operator $cal(U) (t,t_0)$,
$
  ket(alpha","t_0";"t) = cal(U) (t,t_0) ket(alpha","t_0)
$
at $t_0$ we can write
$
  ket(alpha","t_0) = sum_a' c_(a') (t_0) ket(a')
$
similarly for $t$
$
  ket(alpha","t_0";"t) = sum_a' c_(a') (t) ket(a')
$
we don't expect in general that $abs(c_(a') (t)) eq.not abs(c_(a') (t_0))$, but we must have
$
  sum_a' abs(c_(a') (t_0))^2 = sum_a' abs(c_(a') (t))^2
$
so if a state is normalized then it remains normalized,
$
  braket(alpha","t_0) = 1 => braket(alpha","t_0";"t) = 1
$
this is quaranteed if $cal(U) (t,t_0)$ is unitary---so we take this as a fundamental property. We also require that $cal(U) (t,t_0)$ has the composition property,
$
  cal(U) (t_2, t_0) = cal(U) (t_2, t_1) cal(U) (t_1, t_0)"  "(t_2 > t_1 > t_0)
$
We now consider an infinitesimal time-evolution,
$
  ket(alpha","t_0";"t_0 + dd(t)) = cal(U) (t_0 + dd(t),t_0) ket(alpha","t_0)
$
where $lim_(dd(t) arrow 0) cal(U) (t_0 + dd(t),t_0) = 1$. Everything is satisfied by
$
  cal(U) (t_0 + dd(t),t_0) = 1 - i Omega dd(t)
$
with $Omega$ being Hermitian---composition and unitarity can easily be checked.

$Omega$ has units of $s^(-1)$, recalling $E = hbar omega$ and borrowing from classical mechanics that the Hamiltonian is the generator of time evolution we let $Omega = H\/hbar$. Giving,
$
  cal(U) (t_0 + dd(t),t_0) = 1 - (i H dd(t))/hbar
$
where we assume the Hamiltonian $H$ is Hermitian.

=== Schrödinger equation
By composition we can write
$
  cal(U) (t + dd(t),t_0) = cal(U) (t + dd(t),t) cal(U) (t,t_0) = (1 - (i H dd(t))/hbar) cal(U) (t,t_0)
$
or
$
  cal(U) (t + dd(t),t_0) - cal(U) (t,t_0) = -i H/hbar dd(t) cal(U) (t,t_0) => i hbar pdv(, t) cal(U) (t,t_0) = H cal(U) (t,t_0)
$
which is the Schrödinger equation for $cal(U) (t,t_0)$. We can immediately get the equation for a state ket by
$
  i hbar pdv(, t) cal(U)(t,t_0) ket(alpha","t_0) &= H cal(U) (t,t_0) ket(alpha","t_0) \
  i hbar pdv(, t) ket(alpha","t_0";"t) &= H ket(alpha","t_0";"t)
$
of course if we know $cal(U) (t,t_0)$ we don't need the Schrödinger equation, since we can just apply it directly---so we are interested in formal solutions to the Schrödinger equation. We'll just consider the simplest case; a time-independent $H$:
$
  cal(U) (t,t_0) = exp((-i H (t-t_0))/hbar)
$
this can be proven by expanding it, or by considering infinitesimal timesteps $dd(t)$:
$
  lim_(N arrow oo) (1 - ((i H \/hbar)(t-t_0))/N)^N = exp((- i H(t-t_0))/hbar)
$

=== Eigenkets of $H$
We'd like to know how the time-evolution operator acts on $ket(alpha)$, to do this we must know how it acts on the base kets used to expand $ket(alpha)$. This is simple if the base kets are eigenkets of $A$ with $[A,H] = 0$. Then they are also eigenkets of $H$---energy eigenkets:
$
  H ket(a') = E_a' ket(a')
$
taking $t_0 = 0$ we can write
$
  exp((-i H t)/hbar) &= sum_(a') sum_(a'') ketbra(a'') exp((-i H t)/hbar) ketbra(a') \
  &= sum_a' ket(a') exp((-i E_a' t)/hbar) bra(a')
$
now suppose we know $ ket(alpha","t_0=0) = sum_a' ket(a') braket(a', alpha) = sum_a' c_a' ket(a') $ then
$
  ket(alpha","t_0 = 0";"t) = exp((-i H t)/hbar) ket(alpha","t_0=0) = sum_a' ket(a') braket(a', alpha) exp((-i E_a' t)/hbar)
$
so the expansion coefficient evolves as
$
  c_a' (t=0) -> c_a' (t) = c_a' (t=0) exp((-i E_a' t)/hbar)
$
in the case that $ket(alpha","t_0 = 0) = ket(a')$ we get
$
  ket(a","t_0 = 0";"t) = ket(a') exp((-i E_a' t)/hbar)
$
so if the system is initially a simultaneous eigenstate of $A$ and $H$, then it stays like that. So observables compatible with $H$ are _constants of motion_.

This generalizes in the case of multiple mutually compatible observables, all of which also commute with $H$:
$
  exp((-i H t)/hbar) = sum_K' ket(K') exp((-i E_K' t)/hbar) bra(K')
$

=== Expectation values
Assume the initial state is an eigenstate of $A$ which commutes with $H$, we want to find $expval(B)$. We have $ket(a'","t_0 = 0";"t) = cal(U)(t,0) ket(a')$ so
$
  expval(B) & = (bra(a')cal(U)^dagger (t,0)) dot B dot (cal(U) (t,0) ket(a')) \
  & = braket(a', exp((i E_a' t)/hbar) B exp((-i E_a' t)/hbar), a') = braket(a', B, a')
$
so taking the expectation value with respect to an energy eigenstate is independent of time---stationary state.

If instead $ket(alpha";"t_0 = 0) = sum_a' c_a' ket(a')$ we get
$
  expval(B) = sum_a' sum_a'' c_a'^* c_a'' braket(a', B, a'') exp((-i (E_a''-E_a')t)/hbar)
$

=== Spin example
We consider
$
  H = - (e/(m_e c)) bold(S) dot bold(B)
$
assuming a static $bold(B)$ along $hat(bold(z))$ we get
$
  H = - ((e B)/(m_e c)) S_z
$
so $H$ and $S_z$ obviously commute, so the $S_z$ eigenstates are also energy eigenstates, with eigenvalues
$
  E_plus.minus = minus.plus (e hbar B)/(2 m_e c)"  for " S_z plus.minus
$
defining $omega equiv abs(e) B\/m_e c$ we can write $H = omega S_z$. All time evolution is contained in
$
  cal(U) (t,0) = exp((- i omega S_z t)/hbar)
$
we want to apply this to the initial state, which we can write as (for $t=0$),
$
  ket(alpha) = c_+ ket(+) + c_- ket(-)
$
giving
$
  ket(alpha","t_0=0";"t) = c_+ exp((-i omega t)/2) ket(+) + c_- exp((i omega t)/2) ket(-)
$
where we use
$
  H ket(plus.minus) = (plus.minus hbar omega)/2 ket(plus.minus)
$
#pagebreak()
== Schrödinger $arrow.l.r.long.double$ Heisenberg
The quantum dynamics just described is called the Schrödinger picture. Another approach is the Heisenberg picture, where our observables evolve in time.

We have introduced two unitary operators $cal(J) (dd(bold(x)'))$ and $cal(U) (t,t_0)$ which translate our state $ket(alpha) -> U ket(alpha)$. As with all unitary transformations the inner product is not affected,
$
  braket(beta, alpha) -> braket(beta, alpha)
$
we can also consider
$
  braket(beta, X, alpha) -> braket(beta, U^dagger X U, alpha)
$
from the associative axiom it follows that we have two approaches,
$
  ket(alpha) & -> U ket(alpha)", with operators unchanged" \
           X & -> U^dagger X U", with state kets unchanged"
$
so far with the Schrödinger picture we've followed the first approach. In the Heisenberg picture we treat state kets as fixed like they were at some $t_0$---for convenience $t_0 = 0$. We define
$
  cal(U) (t, t_0 = 0) equiv cal(U) (t) = exp((-i H t)/hbar)
$
we then define the Heisenberg picture observable by
$
  A^((H)) (t) equiv cal(U)^dagger (t) A^((S)) cal(U)(t)
$
at $t=0$ the pictures are equivalent
$
  A^((H)) (0) = A^((S))
$
in the Heisenberg picture the state ket is frozen like it were at $t=0$
$
  ket(alpha","t_0=0";"t)_H = ket(alpha","t_0 = 0)
$
as opposed to the Schrödinger picture where
$
  ket(alpha","t_0=0";"t)_S = cal(U)(t) ket(alpha","t_0=0)
$
importantly the expectation value doesn't care about which picture we use
$
  expval(A^((S)), alpha","t_0=0";"t)_S &= expval(cal(U)^dagger A^((S)) cal(U), alpha","t_0=0) \
  &= expval(A^((H)) (t), alpha","t_0 = 0)_H
$

=== Equation of motion
We assume $A^((S))$ does not explicitely depend on time. Then
$
  dv(A^((H)), t) &= pdv(cal(U)^dagger, t) A^((S)) cal(U) + cal(U)^dagger A^((S)) pdv(cal(U), t) \
  &= - 1/(i hbar) cal(U)^dagger H cal(U) cal(U)^dagger A^((S)) cal(U) + 1/(i hbar) cal(U)^dagger A^((S)) cal(U) cal(U)^dagger H cal(U) \
  &= 1/(i hbar) [A^((H)), cal(U)^dagger H cal(U)]
$
in our case $cal(U)$ and $H$ commute so
$
  H = cal(U)^dagger H cal(U)
$
though we might have been tempted to write $H^((H))$, but this is not necessary in what we're doing. So we find
$
  dv(A^((H)), t) = 1/(i hbar) [A^((H)), H]
$
this is the Heisenberg equation of motion, though it was first written by Dirac who used
$
  [,]/(i hbar) -> [,]_"class"
$
however this is of course limited since stuff only in quantum mechanics also satisfies the equation of motion, e.g. spin. In some sense we can _derive_ classical mechanics from the Heisenberg picture---which we derived from the Schrödinger picture, i.e. we derived it using the properties of $cal(U)(t)$ and the defining equation for $A^((H))$.

=== Ehrenfest's theorem
To be able to use any equation of motion we need to know how to build our Hamiltonian. We assume it takes the same form as in classical physics just with operators replacing $x_i$ and $p_i$---this is not always robust, but it works in many cases. We'll typically need (see exercises)
$
  [x_i, F(bold(p))] = i hbar pdv(F, p_i)",  " [p_i, G(bold(x))] = - i hbar pdv(G, x_i)
$
where $F$ and $G$ are expanded in powers of $p_j$ and $x_j$ respectively.

Now we apply the Heisenberg equation of motion to a free particle of mass $m$. The Hamiltonian is taken to be
$
  H = (p_x^2 + p_y^2 + p_z^2)/(2 m)
$
where the operators are assumed to be in the Heisenberg picture. Since $p_i$ commutes with any function of $p_j$ we have
$
  dv(p_i, t) = 1/(i hbar) [p_i, H] = 0
$
so $p_i$ is a constant of motion---more generally whenever $A^((H))$ commutes with $H$ then it is a constant of motion. We also have
$
  dv(x_i, t) = 1/(i hbar) [x_i, H] = 1/(2 m) pdv(, p_i) sum_j^3 p_j^2 = p_i/m = (p_i (0))/m
$
so
$
  x_i (t) = x_i (0) + (p_i (0))/m t
$
note that
$
  [x_i (t), x_i (0)] = [(p_i (0) t)/m, x_i (0)] = - (i hbar t)/m
$
the uncertainty principle gives
$
  expval((Delta x_i)^2)_t expval((Delta x_i)^2)_(t=0) >= (hbar^2 t^2)/(4 m^2)
$
so the position of some particle becomes more and more uncertain with time.

Now we add some potential $V(bold(x))$,
$
  H = bold(p)^2/(2m) + V(bold(x))
$
we obtain
$
  dv(p_i, t)= 1/(i hbar) [p_i, V(bold(x))] = - pdv(, x_i) V(bold(x))
$
and
$
  dv(x_i, t) = p_i/m => dv(x_i, t, 2) = 1/(i hbar)[dv(x_i, t),H] = 1/(i hbar) [p_i/m,H] = 1/m dv(p_i, t)
$
so
$
  m dv(bold(x), t, 2) = - grad V(bold(x))
$
this the quantum version of Newton's second law! Taking the expectatition value with resepct to a Heisenberg state ket---which doesn't move in time---we get
$
  m dv(, t, 2) expval(bold(x)) = dv(expval(bold(p)), t) = - expval(grad V(bold(x)))
$
which is Ehrenfest's theorem---we take the expectation value since now it holds in both pictures.

=== Base kets
In the Schrödinger picture we had the defining equation
$
  A ket(a') = a' ket(a')
$
and notably $A$ doesn't change so neither does the base kets---even though state kets do change. This changes is the Heisenberg picture where
$
  A^((H)) (t) = cal(U)^dagger A(0) cal(U)
$
at $t=0$ the pictures coincide so
$
  cal(U)^dagger A(0) cal(U) cal(U)^dagger ket(a') = a' cal(U)^dagger ket(a') => A^((H)) (cal(U)^dagger ket(a')) = a' (cal(U)^dagger ket(a'))
$
so ${cal(U)^dagger ket(a')}$ form the base kets in the Heisenberg picture, so the base kets change like
$
  ket(a'","t)_H = cal(U)^dagger ket(a')
$
so they satisfy
$
  i hbar pdv(, t) ket(a'","t)_H = - H ket(a'","t)_H
$
but the eigenvalues themselves don't change. As a sanity check consider
$
  A^((H)) (t) & = sum_a' ketbra(a'","t)_H A^((H)) \
              & = sum_a' ket(a'","t)_H a' bra(a'","t)_H \
              & = sum_a' cal(U)^dagger ket(a') a' bra(a') cal(U) \
              & = cal(U)^dagger A^((S)) cal(U)
$
so base kets transforming like this is consistent.

Expansion coefficients are also the same
$
  c_a' (t) & = bra(a') dot (cal(U) ket(alpha","t_0=0)) " Schrödinger picture" \
  c_a' (t) & = (bra(a') cal(U)) dot ket(alpha","t_0=0) " Heisenberg picture"
$
in the Schrödinger picture we take the inner product of a stationary eigenbra with a moving state ket, and in the Heisenberg picture we take the inner product of a moving eigenbra with a stationary state ket.

Similarly we can find the transition amplitude, so the probability for a system in an eigenstate of some $A$ with eigenvalue $a'$, to be in an eigenstate of $B$ with eigenvalue $b'$. In the Schrödinger picture the state ket at $t$ is $cal(U) ket(a')$, while the base kets are constant, so
$
  bra(b') dot (cal(U) ket(a'))
$
in the Heisenberg picture the state ket is constant, but the base kets move oppositely so
$
  (bra(b')cal(U)) dot ket(a')
$
these are obviously the same $braket(b', cal(U) (t,0), a')$.

The basic difference between the two pictures is that in the Schrödinger picture only the state kets evolve, observables and base kets are stationary---while in the Heisenberg picture the state kets are stationary and both observables and base kets evolve, and here base kets evolve "backwards".

#pagebreak()
== Simple Harmonic Oscillator
We want to solve
$
  H = p^2/(2m) + (m omega^2 x^2)/2
$
it turns out to be convenient to use
$
  a = sqrt((m omega)/(2 hbar)) (x + (i p)/(m omega))",  " a^dagger = sqrt((m omega)/(2 hbar)) (x - (i p)/(m omega))
$
these are non-Hermitian, but are each others Hermitian adjoint. The first is the annihilation operator, and the second is the creation operator. It's easy to find
$
  [a,a^dagger] = 1
$
we also define the number operator $N = a^dagger a$, which an be calculated explictely to give
$
  H = hbar omega ( N + 1/2)
$
so $[H,N]=0$ meaning we have simultaneous eigenkets. Denote the energy eigenket of $N$ by $ket(n)$ and
$
  N ket(n) = n ket(n)
$
so
$
  H ket(n) = (n + 1/2) hbar omega ket(n)
  =>
  E_n = (n+1/2) hbar omega
$
note
$
  [N,a] = [a^dagger a, a] = a^dagger [a,a] + [a^dagger, a]a = - a
$
likewise $[N,a^dagger]=a^dagger$. It follows that
$
  N a^dagger ket(n) = ([N,a^dagger] + a^dagger N) ket(n) = (n+1) a^dagger ket(n)
$
and
$
  N a ket(n) = ([N,a]+a N) ket(n) = (n-1) a ket(n)
$
so $a^dagger ket(n)$ is an eigenket of $N$ but with eigenvalue $n+1$, and $a ket(n)$ is an eigenket of $N$ but with eigenvalue $n-1$. Increasing or decreasing $n$ by one amounts to creating or annihilating one unit of energy $hbar omega$---why they are given their respective names.

The previous implies that
$
  a ket(n) = c ket(n-1)
$
note $braket(n, a^dagger a, n) = abs(c)^2$, but $N = a^dagger a$, so
$
  n = abs(c)^2
$
and we obtain
$
  a ket(n) = sqrt(n) ket(n-1)
$
similarly
$
  a^dagger ket(n) = d ket(n+1)
$
note $braket(n, a a^dagger, n) = abs(d)^2$, but $a a^dagger = N+1$, so
$
  n+1 = abs(d)^2
$
and we obtain
$
  a^dagger ket(n) = sqrt(n+1) ket(n+1)
$
now suppose we kept applying the annihilation operator
$
  a^m ket(n) = sqrt(n(n-1)(n-2)dots (n-m)) ket(n - m)
$
We get smaller and smaller eigenkets of $N$ and if $n$ is a positive integer then it terminates. If we start with a non-integer $n$, then the sequence won't terminate and we get negative values of $n$ at some point. This can't happen since
$
  n = braket(n, N, n) = (bra(n)a^dagger) dot (a,n) >= 0
$
so $n$ can't be negative, and therefore it must be a positive integer, and the sequence terminates with $n=0$. Therefore the ground state has
$
  E_0 = 1/2 hbar omega
$
and we can get all states
$
  ket(1) & = a^dagger ket(0) \
  ket(2) & = (a^dagger)/sqrt(2) ket(1) = ((a^dagger)^2/sqrt(2)) ket(0) \
  ket(n) & = ((a^dagger)^n/sqrt(n!)) ket(0)
$
and all energies
$
  E_n = (n +1/2) hbar omega
$
requiring orthonomality for ${ket(n)}$ we also get
$
  braket(n', a, n) = sqrt(n) delta_(n',n-1)",  " braket(n', a^dagger, n) = sqrt(n+1) delta_(n',n+1)
$
with
$
  x = sqrt(hbar/(2 m omega)) (a + a^dagger)",  " p = i sqrt((m hbar omega)/2) (-a + a^dagger)
$
we get
$
  braket(n', x, n) &= sqrt(hbar/(2 m omega)) (sqrt(n) delta_(n',n-1) + sqrt(n+1) delta_(n',n+1)) \
  braket(n', p, n) &= i sqrt((m hbar omega)/2) (- sqrt(n)delta_(n',n-1) + sqrt(n+1) delta_(n',n+1))
$

The ground state is defined by $a ket(0) = 0$, in the $x$-representation
$
  braket(x', a, 0) = sqrt((m omega)/(2 hbar)) bra(x') (x + (i p)/(m omega)) ket(0) = 0
$
or
$
  0 & = braket(x', x, 0) + i/(m omega) braket(x', p, 0) \
    & = x' braket(x', 0) + i /(m omega) (- i hbar) dv(, x') braket(x', 0) \
    & = (x' + x_0^2 dv(, x')) braket(x', 0) = 0
$
with
$
  x_0 equiv sqrt(hbar/(m omega))
$
the solution is
$
  braket(x', 0) = (1/(pi^(1\/4) sqrt(x_0))) exp(-1/2 (x'/x_0)^2)
$
the excited states are then
$
  braket(x', 1) &= braket(x', a^dagger, 0) = 1/(sqrt(2)x_0) (x' - x_0^2 dv(, x)) braket(x', 0) \
  braket(x', 2) &= 1/sqrt(2) braket(x', (a^dagger)^2, 0) = (1/sqrt(2!)) (1/(sqrt(2)x_0))^2 (x'-x_0^2 dv(, x'))^2 braket(x', 0) \
  braket(x', n)&= (1/(pi^(1\/4) sqrt(2^n n!))) 1/(x^(n+1/2)) (x' - x_0^2 dv(, x'))^n exp(-1/2 (x'/x_0)^2)
$

we can find
$
  expval(x^2) = hbar/(2m omega) = x_0^2/2",  " expval(p^2) = (hbar m omega)/2
$
for $n=0$. It follows that
$
           expval(p^2/(2m)) & = expval(T) = (hbar omega)/4 = expval(H)/2 \
  expval((m omega^2 x^2)/2) & = expval(V) = (hbar omega)/4 = expval(H)/2
$
and $expval(x)=expval(p)=0$---which is true for any state. So
$
  expval((Delta x)^2) = expval(x^2)",  " expval((Delta p)^2) = expval(p^2)
$
so we have minimal uncertainty,
$
  expval((Delta x)^2)expval((Delta p)^2) = hbar^2/4
$
for the excited states
$
  expval((Delta x)^2) expval((Delta p)^2) = (n + 1/2)^2 hbar^2
$
which is easy to show.

=== Time-evolution of SHO
Now we work in Heisenberg picture. The equations of motion are
$
  dv(p, t) = - m omega^2 x",   " dv(x, t) = p/m
$
these are equivalent to
$
  dv(a, t) = - i omega a",   " dv(a^dagger, t) = i omega a^dagger
$
whose solutions are just
$
  a(t) = a(0) exp(-i omega t)",   " a^dagger (t) = a^dagger (0) exp(i omega t)
$
or
$
  x(t) + (i p(t))/(m omega) &= x(0) exp(- i omega t) + i (p(0))/(m omega) exp(- i omega t) \
  x(t) - (i p(t))/(m omega) &= x(0) exp(i omega t) - i (p(0))/(m omega) exp(i omega t)
$
equation Hermitian and anti-Hermitian parts
$
  x(t) & = x(0) cos omega t + (p(0))/(m omega) sin omega t \
  p(t) & = - m omega x(0) sin omega t + p(0) cos omega t
$
so the operators oscillate like their classical counterparts.

Another way to derive this would be using
$
  x(t) = exp((i H t)/hbar) x(0) exp((-i H t)/hbar)
$
here we need the Baker-Hausdorff lemma,
$
  exp(i G lambda) A exp(-i G lambda) &= A + i lambda [G,A] + ((i^2 lambda^2)/2!) [G,[G,A]] + dots \
  & dots + ((i^n lambda^n)/n!) [G,[G,[G,dots,[G,A]]] dots] + dots
$

Coherent states\*

#pagebreak()
== The wave-equation
We want to study the time evolution of $ket(alpha","t_0";"t)$ in the $x$-representation---in the Schrödinger picture, or we want to study how $psi(bold(x)', t) = braket(bold(x)', alpha","t_0";"t)$ behaves.

We take the Hamiltonian to be
$
  H = bold(p)^2/(2 m) + V(bold(x))
$
note that
$
  braket(bold(x)'', V(bold(x)), bold(x)') = V(bold(x)') delta^3 (bold(x)'-bold(x)'')
$
from the Schrödinger equation for at state we have
$
  i hbar pdv(, t) braket(bold(x)', alpha","t_0";"t) = braket(bold(x)', H, alpha","t_0";"t)
$
which we can do since in the Schrödinger picture $bra(bold(x)')$ is constant in time. Now we can use
$
  braket(bold(x)', bold(p)^2/(2m), alpha","t_0";"t) = - hbar^2/(2m) nabla'^2 braket(bold(x)', alpha","t_0";"t)
$
and $bra(bold(x)') V(bold(x)) = bra(bold(x)') V(bold(x)')$, to obtain
$
  i hbar pdv(, t) braket(bold(x)', alpha","t_0";"t) = - hbar^2/(2m) nabla'^2 braket(bold(x)', alpha","t_0";"t) + V(bold(x)') braket(bold(x)', alpha","t_0";"t)
$
or in familiar notation
$
  i hbar pdv(, t) psi(bold(x)', t) = - hbar^2/(2m) nabla'^2 psi(bold(x)', t)+V(bold(x)')psi(bold(x)', t)
$
this is the starting point of wave-mechanics.
=== Time-independent wave-equation
We now derive the partial differential equation satisfied by energy eigenfunctions. Recall that the time-dependence of a stationary state is given by
$
  exp((-i E_a' t)/hbar)
$
this lets us write
$
  braket(bold(x)', a'","t_0";"t) = braket(bold(x)', a') exp((- i E_a' t)/hbar)
$
where the system is in a simultaneous eigenstate of $A$ and $H$. Plugging this into the Schrödinger wave-equation we find
$
  - hbar^2/(2m) nabla'^2 braket(bold(x)', a') + V(bold(x)') braket(bold(x)', a') = E_a' braket(bold(x)', a')
$
this equation is satisfied by energy eigenfunctions $braket(bold(x)', a')$ with eigenvalue $E_a'$.

If we pick $A$ to be the function of $bold(x)$ and $bold(p)$ that coincides with $H$ then we can omit $a'$ and just write $ - hbar^2/(2m) nabla'^2 u_E (bold(x)') + V(bold(x)') u_E (bold(x)') = E u_E (bold(x)') $
this is Schrödinger's time-independent wave-equation.

As with any differential equation we need boundary conditions before we can solve it. Take $ E < lim_(abs(bold(x)') arrow oo) V(bold(x)') $
the proper boundary condition in this case is
$
  u_E (bold(x)') arrow 0 "as" abs(bold(x)') arrow oo
$
physically this means our particle is bound. We know from pde's that this boundary condition only yields discrete values of $E-->$ quantization. Likewise if the condition is not satisfied, then we get scattering states with continuous values of $E$.

=== Interpretation
We define the probability density
$
  rho(bold(x)', t) = abs(psi(bold(x)', t))^2 = abs(braket(bold(x)', alpha","t_0";"t))^2
$
from Schrödinger's wave-equation one can then find
$
  pdv(rho, t) + nabla dot bold(j) = 0
$
i.e. the continuity equation, here $bold(j) (bold(x)',t)$ is the probability flux defined by
$
  bold(j) (bold(x)',t) = - (i hbar)/(2 m) (psi^* nabla psi - (nabla psi)^* psi) = hbar/m Im (psi^* nabla psi)
$
$V$ being Hermitian is required for this, so a complex potential would lead to the disappearance of a particle. We can also obtain
$
  integral dd(x, 3) bold(j) (bold(x)',t) = expval(bold(p))_t/m
$
the nature of the continuity equation lead Born to interpret $abs(psi)^2$ as probability.

We can write
$
  psi(bold(x)', t) = sqrt(rho(bold(x)', t)) exp((i S(bold(x)',t))/hbar)
$
with $S$ real, giving
$
  bold(j) = (rho nabla S)/m
$
so the gradient of the phase $S$ characterizes $bold(j)$.

=== Classical limit
We substitute $psi$ written as before into the time-dependent wave-equation
$
  - hbar^2/(2m) & [nabla^2 sqrt(rho) + (2i)/hbar (nabla sqrt(rho)) dot (nabla S) - 1/hbar^2 sqrt(rho) abs(nabla S)^2 + i/hbar sqrt(rho) nabla^2 S] + sqrt(rho) V \
  &= i hbar [ pdv(sqrt(rho), t) + i/hbar sqrt(rho) pdv(S, t)]
$
we assume
$
  hbar abs(nabla^2 S) << abs(nabla S)^2
$
we essentially say that $hbar$ is really small. So we get
$
  1/(2m) abs(nabla S)^2 + V + pdv(S, t) = 0
$
this is just the Hamilton-Jacobi equation, with $S$ being Hamilton's principal function. So in the $hbar arrow 0$ limit we recover classical mechanics.
#pagebreak()
== Solutions to the wave-equation
We'll go through solutions for specific $V(bold(x))$.
=== Free particle
We start with $V(bold(x)) = 0$. The time-independent Schrödinger equation becomes
$
  nabla^2 u_E (bold(x)) = - (2 m E)/hbar^2 u_E (bold(x))
$
we define
$
  bold(k)^2 = k_x^2 +k_y^2 + k_z^2 equiv (2 m E)/hbar^2 = bold(p)^2/hbar^2
$
this is easily solved by $u_E (bold(x)) = u_x (x) u_y (y) u_z (z)$. Giving
$
  (1/u_x dv(u_x, x, 2) + k_x^2) + (1/u_y dv(u_y, y, 2) + k_y^2) + (1/u_z dv(u_z, z, 2) + k_z^2) = 0
$
this has solutions
$
  u_w (w) = c_w e^(i k_w w)" for "w={x,y,z}
$
so we obtain
$
  u_E (bold(x)) = c_x c_y c_z e^(i k_x x + i k_y z + i k_z z) = C e^(i bold(k) dot bold(x))
$
this cannot be normalized in the usual way. Instead we use big-box normalization, where we say all space is within a cube of side $L$---with periodic boundaries so
$
  u_x (x + L) = u_x (x) => k_x L = 2 pi n_x => k_x = (2 pi)/L n_x
$
and similarly for $y$ and $z$. Normalization gives
$
  1 = integral_0^L dd(x) integral_0^L dd(y) integral_0^L dd(z) u_E^* (bold(x)) u_E (bold(x)) = L^3 abs(C)^2 => C = 1/L^(3\/2)
$
so
$
  u_E (bold(x)) = 1/(L^(3\/2)) e^(i bold(k) dot bold(x))
$
with energies
$
  E = bold(p)^2/(2m) = hbar^2/(2m) ((2pi)/L)^2 (n_x^2 + n_y^2 + n_z^2)
$
we can find the density of states $dd(N)\/dd(E)$ by considering a shell in $bold(k)$ space with radius $abs(bold(k)) = 2pi abs(bold(n)) \/L$ and thickness $dd(abs(bold(k))) = 2 pi dd(abs(bold(n))) \/L$. All states in this shell have $E = hbar^2 bold(k)^2 \/2 m$. The number of states $dd(N)$ within the shell is $4 pi bold(n)^2 dd(abs(bold(n)))$---volume of shell. So
$
  dv(N, E) &= (4 pi bold(n)^2 dd(abs(bold(n))))/(hbar^2 abs(bold(k)) dd(abs(bold(k)))\/m) = (4 pi m)/hbar^2 (bold(n)^2)/(2 pi abs(bold(n))\/L) (dd(abs(bold(n))) )/(2 pi dd(abs(bold(n))) \/ L) \
  &= (4 pi m)/hbar^2 (L/(2 pi))^2 abs(bold(n)) \
  &= (4 pi m)/hbar^2 (L/(2 pi))^3 abs(bold(k)) = (4 pi m)/hbar^2 (L/(2 pi))^3 sqrt(2 m E)/hbar \
  &= (m^(3\/2) E^(1\/2) L^3)/(sqrt(2) pi^2 hbar^3)
$
In a more accurate representation the normalization would cancel and give the correct result.

=== SHO
We have $V(x) = m omega^2 x^2 \/2$, giving
$
  - hbar^2/(2m) dv(, x, 2) u_E (x) + 1/2 m omega^2 x^2 u_E (x) = E u_E (x)
$
we define $y equiv x\/x_0$ with $x_0 equiv sqrt(hbar\/m omega)$, and a diemensionless energy $epsilon equiv 2 E \/hbar omega$ giving
$
  dv(, y, 2) u(y) + (epsilon - y^2) u (y) = 0
$
the equation $w''(y) - y^2 w(y) = 0$ has solutions $w(y) prop exp(plus.minus y^2\/2)$, we need the minus sign, else normalizing our solution would be impossible. So we write
$
  u(y) = h(y) e^(- y^2 \/2)
$
where $h(y)$ satisfies
$
  dv(h, y, 2) - 2 y dv(h, y) + (epsilon -1) h(y) = 0
$
Now we introduce generating functions. Consider
$
  g(x,t) & equiv e^(-t^2 + 2 t x) \
         & equiv sum_(n=0)^oo H_n (x) t^n/n!
$
with $H_n (x)$ being the Hermite polynomials---notice that $H_0 (x) = 1$. And
$
  g(0,t) = e^(-t^2) = sum_(n=0)^oo (-1)^n/n! t^(2n)
$
so $H_n (0) = 0$ for odd $n$. For even $n$,
$
  g(0,t) = e^(-t^2) = sum_(n=0)^oo (-1)^(n\/2)/(n\/2)! t^n = sum_(n=0)^oo (-1)^(n\/2)/(n\/2)! n!/n! t^n => H_n (0) = ((-1)^(n\/2) n!)/(n\/2)!
$
and $H_n (-x) = (-1)^n H_n (x)$. By taking two different derivatives we can find
$
  H'_n (x) = 2 n H_(n-1) (x)
$
with this we can build all $H_n$.

This is relevant since by taking time-derivatives we find
$
  H_(n+1) (x) = 2 x H_n (x) - 2n H_(n-1) (x)
$
or using the previous recursion relation
$
  H''_n (x) = 2 x H'_n (x) - 2 n H_n (x)
$
which is equivalent to the transformed Schrödinger equation with $epsilon -1 = 2n$. So
$
  u_n (x) = c_n H_n (x sqrt((m omega)/hbar)) e^(- m omega x^2 \/2 hbar)
$
with $c_n$ being found by
$
  integral_(-oo)^oo H_n (x) H_m (x) e^(-x^2) dd(x) = pi^(1\/2) 2^n n! delta_(n m)
$

=== Linear potential
The linear potential is $V(x) = k abs(x)$. This potential has a classical turning point at some $x = a$ with $E = k a$.

The Schrödinger equation becomes
$
  - hbar^2/(2m) dv(u_E, x, 2) + k abs(x) u_E (x) = E u_E (x)
$
we can just treat $x >= 0$ since $V(-x) = V(x)$. We have two types of solutions $u_E (-x) = plus.minus u_E (x)$, in both cases $u_E (x) arrow 0$ as $x arrow oo$. If $u_E (-x) = - u_E (x)$ then $u_E (0) = 0$. If $u_E (-x) = u_E (x)$ then $u'_E (0) = 0$, since $u_E (epsilon) - u_E (-epsilon) equiv 0$ for $epsilon arrow 0$---these are referred to as even or odd parity.

We define
$
  x_0 = ((hbar^2)/(m k))^(1\/3) "and" E_0 = k x_0 = ((hbar^2 k^2)/m)^(1\/3)
$
giving dimensionless $y equiv x\/x_0$ and $epsilon equiv E\/E_0$, and we obtain
$
  dv(u_E, y, 2) - 2 (y-epsilon) u_E (y) = 0 "for" y>=0
$
note $y = epsilon$ when $x = E\/k$---the classical turning point $a$. We can define $z equiv 2^(1\/3) (y-epsilon)$ giving
$
  dv(u_E, z, 2) - z u_E (z) = 0
$
which is the Airy equation, with the solution being the Airy function $"Ai"(z)$. The boundary conditions becomes zeroes for $"Ai"'(z)$ and $"Ai" (z)$ with $z = - 2^(1\/3) epsilon$. These determine the quantized energies.

This potential actually corresponds to a quark-antiquark bound system with $k tilde r$. It also corresponds to the quantum bouncing ball with $k = m g$---this is only for $x >= 0$ as there is an infinite potential barrier at $x = 0$ causing the bounce---this means that only odd parity solutions are allowed.

=== WKB
The WKB (Wentzel, Kramers and Brillouin) approximation is a useful technique which uses the linear potential to join solutions near turning points.

We can write
$
  dv(u_E, x, 2) + (2 m)/hbar^2 (E - V(x)) u_E (x) = 0
$
we define
$
  k(x) &equiv [(2 m)/hbar^2 (E-V(x))]^(1\/2) "for" E > V(x) \
  k(x) equiv - i kappa(x) &equiv - i [(2 m)/hbar^2 (V(x)-E)]^(1\/2) "for" E < V(x)
$
so
$
  dv(u_E, x, 2) + k(x)^2 u_E (x) = 0
$
we assume $V(x)$ varies slowly and try a solution of the form
$
  u_E (x) equiv exp(i W(x) \/hbar)
$
giving
$
  i hbar dv(W, x, 2) - (dv(W, x))^2 + hbar^2 k(x)^2 = 0
$
varying slowly is quantified by the condition
$
  hbar abs(dv(W, x, 2)) << abs(dv(W, x))^2
$
this gives a lowest-order approximation for $W(x)$
$
  W'_0 (x) = plus.minus hbar k(x)
$
a first-order approximation is then obtained
$
  (dv(W_1, x))^2 & = hbar^2 k(x)^2 + i hbar W''_0 (x) \
                 & = hbar^2 k(x)^2 plus.minus i hbar^2 k' (x)
$
so
$
  W(x) approx W_1 (x) &= plus.minus hbar integral^x dd(x)' [k^2 (x') plus.minus i k' (x')]^(1\/2) \
  &approx plus.minus hbar integral^x dd(x)' k(x') [1 plus.minus i/2 (k'(x'))/(k^2 (x'))] \
  &= plus.minus hbar integral^x dd(x') k(x') + i/2 hbar ln(k(x))
$
the WKB approximation is then
$
  u_E (x) approx exp[i W(x) \/hbar] = 1/sqrt(k(x)) exp(plus.minus i integral^x dd(x)' k(x'))
$
this specifies the solutions for $E > V$ and $E < V$. We don't care about the joining procedure.

Instead consider a potential well with turning points $x_1$ and $x_2$ creating three regions. In the middle region the wave function behaves like our approximation with the first $k(x)$ and in the two outer regions with the second $k(x)$. In the neighborhood of the turning points the solutions are given by Airy function, since we assume a linear approximation in those regions. This leads to a consistency check
$
  integral_(x_1)^(x_2) dd(x) sqrt(2 m (E-V(x))) = (n+1/2) pi hbar
$
this gives approximate expressions for the energy levels. Again consider the bounding ball with
$
  V = cases(m g x "  for" x > 0, oo "     for" x < 0)
$
where $x$ is the height from the surface. We could use $x_1 = 0$ and $x_2 = E\/m g$, this is the classical turning points, but our wave function leaks into $x < x_1$ region, even though we require it must vanish. For this we use tha odd-parity solutions which vanish at $x = 0$.  So
$
  V(x) = m g abs(x)
$
with turning points $x_1 = - E\/m g$ and $x_2 = E\/m g$. Then
$
  integral_(- E\/m g)^(E\/m g) dd(x) sqrt(2m (E-m g abs(x))) = (n_"odd" + 1/2) pi hbar
$
or
$
  integral_0^(E\/m g) dd(x) sqrt(2 m (E- m g x)) = (n - 1/4) pi hbar
$
giving
$
  E_n = {[3(n - 1\/4) pi]^(2\/3)/2} (m g^2 hbar^2)^(1\/3)
$
\* interpretation of WKB limit.

The WKB limit is equivalent to
$
  lambda = hbar/sqrt(2 m [E-V(x)]) << (2 (E-V(x)))/abs(dd(V)\/dd(x))
$

== Propagators & the Path Integral
=== Propagators
We know that we can write
$
  ket(alpha","t_0";"t) = sum_a' ket(a') braket(a', alpha","t_0) exp[(-i E_a' (t-t_0))/hbar]
$
in terms of the wavefunction
$
  psi(bold(x)', t) = sum_a' c_a' (t_0) u_a' (bold(x)') exp[(-i E_a' (t-t_0))/hbar]
$
with $u_a' (bold(x)') = braket(bold(x)', a')$ and
$
  c_a' (t_0) = braket(a', alpha","t_0) = integral dd(x', 3) u_a'^* (bold(x)') psi(bold(x)', t_0)
$
we can write
$
  psi(bold(x)'', t) = integral dd(x', 3) K(bold(x)'',t;bold(x)',t_0) psi(bold(x)', t_0)
$
where the propagator $K$ is
$
  K(bold(x)'',t;bold(x)',t_0) = sum_a' braket(bold(x)'', a') braket(a', bold(x)') exp[(-i E_a' (t-t_0))/hbar]
$
note that if $bold(x)'$ and $t_0$ is fixed then the propagator satisfies the Schrödinger wave equation, and in the limit $t arrow t_0$ we have
$
  lim_(t arrow t_0) K(bold(x)'',t;bold(x)',t_0) = delta^3 (bold(x)''-bold(x)')
$
this looks a lot like a Green function. The propagator is exactly the Greens function for the Schrödinger wave equation, since it satisfies
$
  [- hbar^2/(2m) nabla''^2 + V(bold(x)'') - i hbar pdv(, t)] K(bold(x)'',t;bold(x)',t_0) = - i hbar delta^3 (bold(x)''-bold(x)') delta(t-t_0)
$
given $K = 0$ for $t < t_0$.

Let $t_0 = 0$ and $bold(x)''=bold(x)'$, then
$
  G(t) &equiv integral dd(x', 3) K(bold(x)',t;bold(x)',0) \
  &= integral dd(x', 3) sum_a' abs(braket(bold(x)', a'))^2 exp((-i E_a' t)/hbar) \
  &= sum_a' exp((-i E_a' t)/hbar)
$
this looks a lot like a partition function, we can write
$
  Z = sum_a' exp(- beta E_a') "with" beta=(i t)/hbar
$
we can also find
$
  tilde(G) (E) = sum_a' 1/(E-E_a')
$
where $tilde(G)$ is the Laplace-Fourier transform of $G$.

We can write the propagator as
$
  K(bold(x)'',t;bold(x)',t_0) &= sum_a' braket(bold(x)'', exp((-i H t)/hbar), a') braket(a', exp((i H t_0)/hbar), bold(x)') \
  &= braket(bold(x)'', cal(U)(t,t_0), bold(x)') \
  &= braket(bold(x)''","t, bold(x)'","t_0)
$
with $ket(bold(x)'","t_0)$ and $bra(bold(x)''","t)$ being an eigenket and an eigenbra of the position operator in the Heisenberg picture. This shows us that we can identify the propagator as the transition amplitude for a particle to go from $(bold(x)',t_0) -> (bold(x)'',t)$. In fact we can split our interval say $(t',t''')$ into $(t',t'')$ and $(t'',t''')$ to write
$
  braket(bold(x)'''","t''', bold(x)'","t') = integral dd(x'', 3) braket(bold(x)'''","t''', bold(x)''","t'') braket(bold(x)''","t'', bold(x)'","t')
$
since
$
  integral dd(x'', 3) braket(bold(x)''","t'', bold(x)''","t'') = 1
$
this can of course be done for more subdivisions---the composition property.

=== Path integrals in the book
We stay in one-dimension. We consider the transition amplitude for $(x_1,t_1) arrow (x_N, t_N)$, with $N-1$ equal intervals of length
$
  t_j - t_(j-1) = Delta t = (t_N-t_1)/(N-1)
$
then by the composition property
$
  braket(x_N","t_N, x_1","t_1) &= integral dd(x_(N-1)) dots integral dd(x_2) braket(x_N","t_N, x_"N-1"","t_"N_1") dots braket(x_2","t_2, x_1","t_1)
$
this means that for each $t_(n-1) arrow t_n$ we integrate over $x_2,x_3,dots, x_(N-1)$, so we essentially sum over all possible paths. In classical mechanics the path taken by a particle is unique and determined by the least action principle
$
  delta integral_(t_1)^(t_2) dd(t) L_"classical" (x,dot(x)) = 0
$
in quantum mechanics however we must consider every possible path---but we should still retrieve the classical path in the limit $hbar -> 0$. This is the problem Feynman tried to solve.

We define
$
  S(n,n-1) equiv integral_(t_(n-1))^t_n dd(t) L_"classical" (x,dot(x))
$
consider some small path segment $(x_(n-1),t_(n-1)) -> (x_n,t_n)$, following Dirac we associate $exp[(i S(n,n-1))\/hbar]$ to that segment. Then going along some path we get
$
  product_(n=2)^N exp[(i S(n,n-1))/hbar] = exp[i/hbar sum_(n=2)^N S(n,n-1)] = exp[(i S(N,1))/hbar]
$
this gives a contribution to $braket(x_N","t_N, x_1","t_1)$ following a particular path---so we can write
$
  braket(x_N","t_N, x_1","t_1) tilde sum_"all paths" exp[(i S(N,1))/hbar]
$
in the limit $hbar arrow 0$ then the exponential will oscillate a bunch, leading to many paths cancelling. But if a path satisfies $delta S(N,1) = 0$---i.e. it would be the classically correct path---then slight deformations don't change $exp(i S\/hbar)$, so in the $hbar -> 0$ limit this leads to constructive interference with the classical path being singled out and all others vanishing.

To make this precise we write
$
  braket(x_n","t_n, x_(n-1)","t_(n-1)) = 1/w(Delta t) exp[(i S(n,n-1))/hbar]
$
where we assume $t_n-t_(n-1)$ is an infinitesimal. We'd like to determine the weight factor $w(Delta t)$. Due to $Delta t -> 0$ we can treat the path between $t_(n-1) -> t_n$ as a straight line
$
  S(n,n-1) & = integral_(t_(n-1))^(t_n) dd(t) [(m dot(x)^2)/2 - V(x)] \
           & = Delta t { m/2 [(x_n-x_(n-1))/(Delta t)]^2 - V((x_n +x_(n-1))/2)}
$
since the weight doesn't depend on $V(x)$ we find it for the free particle. For $V = 0$ we have
$
  braket(x_n","t_n, x_(n-1)","t_(n-1)) = 1/w(Delta t) exp[(i m(x_n-x_(n-1))^2)/(2 hbar Delta t)]
$
in the limit $t_n = t_(n-1)$ this should be $delta(x_n - x_(n-1))$ giving
$
  1/w(Delta t) = sqrt(m/(2 pi i hbar Delta t))
$
so we have
$
  braket(x_n","t_n, x_(n-1)","t_(n-1)) = sqrt(m/(2 pi i hbar Delta t)) exp[(i S(n,n-1))/hbar]
$
so
$
  braket(x_N","t_N, x_(N-1)","t_(N-1)) &= lim_(N arrow oo) (m/(2 pi i hbar Delta t))^((N-1)\/2) integral dd(x_(N-1)) dots integral dd(x_2) product_(n=2)^N exp[(i S(n,n-1))/hbar]
$
by convention we define the operator
$
  integral_(x_1)^(x_N) cal(D) [x(t)] equiv lim_(N -> oo) (m/(2 pi i hbar Delta t))^((N-1)\/2) integral dd(x_(N-1)) dots integral dd(x_2)
$
so
$
  braket(x_N","t_N, x_(N-1)","t_(N-1)) &= integral_(x_1)^(x_N) cal(D) [x(t)] exp[i integral_(t_1)^(t_N) dd(t) (L_"classical" (x,dot(x)))/hbar] \
  &= integral_(x_1)^(x_N) cal(D) [x(t)] exp[(i S(N,1))/hbar]
$
this is the Feynman path integral.

The last thing we'll do is show that this formulation is equivalent to Schrödinger's wave mechanics---by showing that the path integral for $braket(x_N","t_N, x_1","t_1)$ satisfies the Schrödinger equation, and is just the propagator. \*

=== Path integrals in the notes
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

#pagebreak()
== Gauge Transformations
Consider adding a constant potential $V_0$ to some $V(bold(x))$. How is the state ket $ket(alpha","t_0";"t)$ for $V(bold(x))$ related to the state ket $tilde(ket(alpha","t_0";"t))$ for $tilde(V(bold(x)))$?---we assume they coincide at $t = t_0$, and obtain
$
  tilde(ket(alpha","t_0";"t)) &= exp[-i (bold(p)^2/(2m) + V(x) + V_0) (t-t_0)/hbar] ket(alpha) \
  &= exp[(-i V_0 (t-t_0))/hbar] ket(alpha","t_0";"t)
$
so $V_0$ just adds a phase. For a stationary state with time-dependence $ exp[(-i E(t-t_0))/hbar] =>^"equivalent" E -> E+V_0 $
this is an example of a gauge transformation. For $V(bold(x))+V_0 (t)$ we find
$
  ket(alpha","t_0";"t) -> exp[- i integral_(t_0)^(t) dd(t') (V_0 (t'))/hbar] ket(alpha","t_0";"t)
$
this has real measurable effects---which are purely quantum mechanical---see pg. 122-125.

=== Electromagnetism
We have
$
  bold(E) = - grad phi",  " bold(B) = curl bold(A)
$
classically
$
  H = 1/(2 m) (bold(p) - (e bold(A))/c)^2 + e phi
$
in QM we say $phi$ and $bold(A)$ are functions of the operator $bold(x)$---so $bold(p)$ and $bold(A)$ don't commute. So we write
$
  (bold(p) - (e bold(A))/c)^2 -> p^2 - (e/c) (bold(p) dot bold(A) + bold(A) dot bold(p)) + (e/c)^2 bold(A)^2
$
we can obtain---in the Heisenberg picture
$
  dv(x_i, t) = [x_i,H]/(i hbar) = (p_i - e A_i\/c)/m
$
we denote $bold(p)$ by the canonical momentum, and this new quantity by the kinematical momentum
$
  bold(Pi) equiv m dv(bold(x), t) = bold(p) - (e bold(A))/c
$
we have
$
  [Pi_i,Pi_j] = (i hbar e)/c epsilon_(i j k) B_k
$
we can also obtain the QM equivalent Lorentz law
$
  m dv(bold(x), t, 2) = dv(bold(Pi), t) = e [bold(E) + 1/(2 c) (dv(bold(x), t) times bold(B)-bold(B) times dv(bold(x), t))]
$
sandwiching $H$ between $bra(bold(x)')$ and $ket(alpha","t_0";"t)$ gives
$
  1/(2m) [- i hbar nabla' - (e bold(A) (bold(x)'))/c ] dot [- i hbar nabla' - (e bold(A) (bold(x)'))/c] braket(bold(x)', alpha","t_0";"t) + e phi (bold(x)') braket(bold(x)', alpha","t_0";"t) = i hbar pdv(, t) braket(bold(x)', alpha","t_0";"t)
$
this is minimal coupling and amounts to $bold(p) -> bold(Pi)$ in the Schrödinger equation.

The common gauge transformation is given by
$
  phi -> phi - 1/c pdv(Lambda, t)",  " bold(A) -> bold(A) + grad Lambda
$
which leaves $bold(E)$ and $bold(B)$ invariant, given
$
  bold(E) = - grad phi - 1/c pdv(bold(A), t)",  " bold(B) = curl bold(A)
$
but we don't concern our selves with time-dependent fields or potentials so the time derivative is zero. Classically $bold(p)$ is not gauge-invariant, but $bold(Pi)$ is.

In QM we demand that the expectation values behave similarly to the classical case, so $expval(bold(x))$ and $expval(bold(Pi))$ should be gauge-invariant while $expval(bold(p))$ is expected to change. Let $tilde(bold(A)) = bold(A) + grad Lambda$---then we require
$
  braket(alpha, bold(x), alpha) = braket(tilde(alpha), bold(x), tilde(alpha))
$
and
$
  expval((bold(p)-(e bold(A))/c), alpha) = expval((bold(p)-(e tilde(bold(A)))/c), tilde(alpha))
$
the norm should also be preserved $braket(alpha) = braket(tilde(alpha))$. We want to construct an operator $cal(G)$ defined by
$
  ket(tilde(alpha)) = cal(G) ket(alpha)
$
it must therefore satisfy
$
  cal(G)^dagger bold(x) cal(G) &= bold(x) \
  cal(G)^dagger (bold(p) - (e bold(A))/c - (e grad Lambda)/c) cal(G) &= bold(p) - (e bold(A))/c
$
we claim
$
  cal(G) = exp[(i e Lambda(bold(x)))/(hbar c)]
$
this is unitary so that's a good start, and note
$
  exp[(-i e Lambda)/(hbar c)] bold(p) exp[(i e Lambda)/(hbar c)] &= exp[(-i e Lambda)/(hbar c)] (bold(p), exp[(i e Lambda)/(hbar c)]) + bold(p) \
  &= - exp[(-i e Lambda)/(hbar c)] i hbar grad (exp[(i e Lambda)/(hbar c)]) + bold(p) \
  &= bold(p) + (e grad Lambda)/c
$
so a gauge transformation just corresponds to adding a phase factor $cal(G) = exp[i e Lambda (bold(x))\/hbar c]$---this can also be shown less elegantly using the Schrödinger equation directly. We find
$
  tilde(psi) (bold(x)',t) = exp[(i e Lambda(bold(x)'))/(hbar c)] psi(bold(x)', t)
$
notably the probability flux is gauge-invariant.

\*Aharanov-Bohm effect and magnetic monopoles.


