//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Formalism
The Stern-Gerlach experiment illustrates why we need quantum mechanics. And understanding two-state systems as in the Stern-Gerlach experiemnt is essential.

As should be familiar the Stern-Gerlach experiment uses $bold(mu) prop bold(S)$ and $F_z tilde mu_z partial_z B_z$. We shoot a neutral beam through a magnetic field and it deflects depending on the value of $mu_z$ or equivalently $S_z$. Classically we expect a continuous spread of $mu_z$ meaning the beam would smear. We instead observe splitting corresponding to spin up and spin down. Sequential Stern-Gerlach experiments show we can not measure $S_x$ and $S_z$ simultaneously etc.

== Dirac notation
All physical states are represented by state vectors in a complex vector space. Given the vector space is a complete inner product space we call it a Hilbert space $cal(H)$. We call these states kets and denote them by $ket(alpha)$. We postulate that $ket(alpha)$ and $c ket(alpha)$ correspond to the same state $ket(alpha) tilde c ket(alpha)$.

All observables are represented by operators acting on the vector space. For some operator $A$ we have $A ket(alpha) = ket(beta)$. We define eigenkets by the property
$
  A ket(a') = a' ket(a')
$
with $a'$ being the eigenvalue. We can expand any $ket(alpha)$ in terms of the eigenkets of some operator $A$
$
  ket(alpha) = sum_a' c_a' ket(a')
$
if the vector space is finite.

We define the bra space as the dual to the ket space. We postulate that all $ket(alpha)$ have a corresponding bra denoted by $bra(alpha)$. We have the dual correspondence
$
  c_alpha ket(alpha) <-->^"DC" c_alpha^* bra(alpha)
$
We define the innerproduct by $braket(beta, alpha)$. We postulate $braket(beta, alpha) = braket(alpha, beta)^*$ implying
$
  braket(alpha, alpha) >= 0
$

$ket(alpha)$ and $ket(beta)$ are orthogonal if
$
  braket(alpha, beta)=0
$
General operators $X dots$ act like $X ket(alpha)$ and  $bra(beta) X$. We only consider linear operators
$
  X (c_alpha ket(alpha) + c_beta ket(beta)) = c_alpha X ket(alpha) + c_beta X ket(beta)
$
and define the Hermitian adjoint as the dual operator $ X ket(alpha) <->^"DC" bra(alpha)X^dagger $ if $X = X^dagger$ we say $X$ is Hermitian. We can multiply operators and $(X Y)^dagger = Y^dagger X^dagger$ since
$
  X Y ket(alpha) = X(Y ket(alpha)) <->^"DC" (bra(alpha)Y^dagger)X^dagger = bra(alpha) Y^dagger X^dagger
$
We define the outer product by $ketbra(beta, alpha)$. This is easily seen to be an operator.

We assume the associative axiom of multiplication. Meaning kets, bras and operators are all associative under multiplication. This is used to show $X = ketbra(beta, alpha)$ is an operator. We have $X^dagger = ketbra(alpha, beta)$ and $braket(beta, X, alpha) = braket(alpha, X^dagger, beta)^*$.

== Hermiticity and representations
#thm[
  The eigenvalues of a Hermitian operator $A$ are real and the eigenkets of $A$ corresponding to different eigenvalues are orthogonal.
]
#proof[

  We have $A ket(a') = a' ket(a')$ and $bra(a'')A^dagger = bra(a'')A = a''^* bra(a'')$. Combining we obtain $ (a'-a''^*) braket(a'', a') = 0 $
  Let $a''=a'$ then $a'=a'^*$. Assuming $a'eq.not a''$ then $a'-a''^*=a'-a'' eq.not 0$ by assumption and we must have $braket(a'', a')=0$.
]

This theorem is very important since we want observables to have real eigenvalues. This is then satisfied if the corresponding operator is Hermitian!

We normalize $ket(a')$ such that $braket(a'', a') = delta_(a'' a')$. By construction ${ket(a')}$ is complete. Then for any $ket(alpha)$ we can write
$
  ket(alpha) = sum_(a') c_a' ket(a') = sum_(a') ket(a') braket(a', alpha)
$
where the second equality can be seen by $braket(a'', alpha)$. We find the completeness relation
$
  sum_(a') ketbra(a') = 1
$
which is very useful. As an example we can find
$
  braket(alpha, alpha) = sum_(a') underbracket(abs(braket(a', alpha))^2, abs(c_a')^2) =^"normalized" 1
$
which is nice! We define the projection operator along $ket(a')$ by
$
  Lambda_(a') equiv ketbra(a', a')
$
since it picks out the part along $ket(a')$. The completeness relation then becomes
$
  sum_(a') Lambda_(a') = 1
$

Consider
$
  X = sum_(a'') sum_(a') ket(a'') underbracket(braket(a'', X, a'), "matrix element") bra(a')
$
we write the elements in a matrix
$
  X eq^dot mat(braket(a^((1)), X, a^((1))), braket(a^((1)), X, a^((2))), dots; braket(a^((2)), X, a^((1))), braket(a^((2)), X, a^((2))), dots; dots.v, dots.v, dots.down)
$
this is called the matrix representation of $X$. Assuming $X = A$ is an observable with eigenkets $ket(a')$ then
$
  A = sum_(a'') sum_(a') ket(a'') braket(a'', A, a') bra(a') = sum_(a') a' Lambda_(a')
$
so it becomes diagonal!

As an example consider a spin $1/2$ system. We use $ket(S_z\;plus.minus) = ket(plus.minus)$ with eigenvalues $plus.minus hbar/2$. The completeness relation becomes
$
  bb(1) = ketbra(+, +) + ketbra(-, -)
$
and we can write $S_z$ as
$
  S_z = hbar/2 (ketbra(+, +)-ketbra(-, -))
$
We define the ladder operators
$
  S_+ equiv hbar ketbra(plus, minus)";  " S_- equiv hbar ketbra(-, +)
$
and their action is obvious. The matrix representation is now easily found
$
  ket(+) eq^dot vec(1, 0)",  " ket(-) eq^dot vec(0, 1)
$
$
  S_z eq^dot hbar/2 mat(1, 0; 0, -1)",  " S_+ eq^dot hbar mat(0, 1; 0, 0)",  " S_- eq^dot hbar mat(0, 0; 1, 0)
$

== Measurement and compatibility
Before a measurement of some observable $A$ we can write the state $ket(alpha)$ as
$
  ket(alpha) = sum_(a') c_(a') ket(a') = sum_(a') ket(a') braket(a', alpha)
$
when we make a measurement the system collapses into one of the eigenstates of $A$
$ ket(alpha) -->^"m.m" ket(a') $
typically this changes the state unless $ket(a') -->^"m.m" ket(a')$. As a result we measure $a'$ i.e. an eigenvalue of the observable!

We assume the probability of measuring $a'$ is $abs(braket(a', alpha))^2$ given $ket(alpha)$ is normalized. This is taken as an axiom of the theory. As a consequence repeated measurement gives the same value since $braket(a') = 1$. Likewise given $braket(a'', a')=0$ the probability of $ket(a') -->^"m.m" ket(a'')$ is zero. This interpretation makes sense since the probability is non-negative and sums to one.



We define the expectation value of $A$ with respect to $ket(alpha)$ as
$
  expval(A) equiv expval(A, alpha) #h(1em) (= expval(A)_alpha)
$
we can write
$
  expval(A) = sum_(a') sum_(a'') braket(alpha, a'') braket(a'', A, a') braket(a', alpha) = sum_(a') a' abs(braket(a', alpha))^2
$
matching what we expect from an expectation value!

Consider a spin $1/2$ system. By doing sequential Stern-Gerlach expreiments we can show
$
  ket(S_x\; plus.minus) = 1/sqrt(2) (ket(+) plus.minus ket(-))";  " ket(S_y\; plus.minus) = 1/sqrt(2) (ket(+) plus.minus i ket(-))
$
we can write $S_x$ and $S_y$ as
$
  S_x & = hbar/2 (ketbra(+, -) + ketbra(-, +))";  " S_y = hbar/2 (-i ketbra(+, -) + i ketbra(-, +))
$
and $S_plus.minus$ becomes
$
  S_plus.minus = S_x plus.minus i S_y
$
The $S_i$ satify the relations
$
  [S_i,S_j] = i epsilon_(i j k) hbar S_k";  " {S_i, S_j} = 1/2 hbar^2 delta_(i j)
$
We define $bold(S)^2 equiv bold(S) dot bold(S) = S_x^2 + S_y^2 + S_z^2$. By the anticommutation relation $ bold(S)^2 = 3/4 hbar^2 $
implying $[bold(S)^2, S_i]=0$.

We say $A$ and $B$ are compatible if $[A,B]=0$. We want a relation between the eigenkets of compatible operators.

#thm[
  Assume $[A,B]=0$ and non-degeneracy. Then all elements $braket(a'', B, a')$ are diagonal.
]
#proof[
  By definition
  $
    braket(a'', [A,B], a') = (a''-a') braket(a'', B, a') = 0
  $
  so $braket(a'', B, a') = 0$ unless $a'' = a'$.
]

So $A$ and $B$ are both diagonal with the same eigenkets. We write
$
  B = sum_(a'' a') ket(a'') braket(a'', B, a') bra(a') = sum_a'' ket(a'') braket(a'', B, a'') bra(a'')
$
then
$
  B ket(a') = sum_(a'') ket(a'') braket(a'', B, a'') braket(a'', a') = braket(a', B, a') ket(a')
$
so $ket(a')$ is a simultaneous eigenket of $A$ and $B$. We write this as $ket(a'b') equiv overbracket(ket(K'), "collective" #linebreak() "index")$. This generalizes to more mutually compatible observables and degeneracies.

Compatible observables are nice since they do not interfere upon measurement. This follows from
$
  A B ket(a' b') = a' b' ket(a' b') = B A ket(a' b')
$
which implies $[A, B] = 0$. We see that any incompatible observables interfere upon measurement. Then as the example shows we can not measure $S_i$ and $S_j$ simultaneously for $i eq.not j$. While we can find a simultaneous eigenket for $S_i$ and $bold(S)^2$.

== Uncertainty
We define the operator $Delta A equiv A-expval(A)$ giving the dispersion
$ expval((Delta A)^2) = expval(A^2)-expval(A)^2 $
this measures the _fuzziness_ of $A$.
#thm[
  Let $A$ and $B$ be observables. Then for any state
  $
    expval((Delta A)^2)expval((Delta B)^2) >= 1/4 abs(expval([A,B]))^2
  $
]
#proof[
  Define
  $
    ket(alpha) equiv Delta A ket(dot)";  " ket(beta) equiv Delta B ket(dot)
  $
  where $ket(dot)$ is any state.

  From Cauchy-Schwarz we have
  $
    expval((Delta A)^2)expval((Delta B)^2) >= abs(expval(Delta A Delta B))^2
  $
  here we use hermiticity of $Delta A$ and $Delta B$. We write  $ Delta A Delta B = 1/2 underbracket([Delta A, Delta B], "anti-Hermitian") + 1/2 underbracket({Delta A, Delta B}, "Hermitian") $
  with $[Delta A, Delta B] = [A, B]$. Taking $expval(dots)$ we find
  $
    expval(Delta A Delta B) = 1/2 underbracket(expval([A,B]), "imaginary") + 1/2 underbracket(expval({Delta A, Delta B}), "real")
  $
  Then
  $
    abs(expval(Delta A Delta B))^2 = 1/4 abs(expval([A,B]))^2 + underbracket(1/4 abs(expval({Delta A, Delta B}))^2, >= 0 " "-> "remove it")
  $
  and we are done.
]

This is the uncertainty principle and its importance should be obvious.

== Transformation operator
Consider $A$ and $B$ with $[A, B] eq.not 0$. We want a transformation between the two bases ${ket(a')}$ and ${ket(b')}$.

#thm[
  Let ${ket(a')}$ and ${ket(b')}$ be two bases satisfying orthonormality and completeness. Then a unitary operator $U$ with
  $
    ket(b^((1))) = U ket(a^((1))), ket(b^((2))) = U ket(a^((2))),dots
  $
  exists.
]

#proof[
  We claim
  $
    U = sum_k ketbra(b^((k)), a^((k)))
  $
  Then
  $
    U ket(a^((l)))=ket(b^((l)))
  $
  and
  $
    U^dagger U = sum_k sum_l ket(a^((l))) braket(b^((l)), b^((k))) bra(a^((k))) = sum_k ketbra(a^((k)))=1
  $
  so it is unitary.
]

The matrix elements are simple
$
  braket(a^((k)), U, a^((l))) = braket(a^((k)), b^((l)))
$
these are nice!

Consider an expansion
$
  ket(alpha) = sum_l ket(a^((l))) braket(a^((l)), alpha)
$
Then
$
  braket(b^((k)), alpha) = sum_l braket(b^((k)), a^((l))) braket(a^((l)), alpha) = sum_l braket(a^((k)), U^dagger, a^((l))) braket(a^((l)), alpha)
$
and
$
  braket(b^((k)), X, b^((l))) &= sum_m sum_n braket(b^((k)), a^((m)))braket(a^((m)), X, a^((n)))braket(a^((n)), b^((l))) \
  &= sum_m sum_n braket(a^((k)), U^dagger, a^((m))) braket(a^((m)), X, a^((n)))braket(a^((n)), U, a^((l)))
$
We define the trace by
$
  tr(X) & = sum_a' braket(a', X, a')
$
this has some nice properties worked out in the exercises e.g.
$
  tr(ketbra(b', a')) & equiv braket(a', b')
$

Consider
$
  B ket(b') = b' ket(b')
$
we write this as
$
  sum_a' braket(a'', B, a') braket(a', b') = b' braket(a'', b')
$
Then if $ket(b')$ is the $l$th eigenket of $B$ this becomes
$ sum_j B_(i j) C_j^((l)) = b^((l)) C_i^((l)) $
where $B_(i j) equiv braket(a^((i)), B, a^((j)))$ and $C_k^((l)) equiv braket(a^((k)), b^((l))) tilde U_(k l)$. Solutions are given by
$
  det(B- lambda bb(1)) = 0
$
with the roots being $b^((l))$.

#thm[
  Let ${ket(a')}$ and ${ket(b')}$ be connected by $U$. Then $A$ and $U A U^dagger$ are unitary equivalent.
]
#proof[
  We have
  $
    A ket(a^((l))) = a^((l)) ket(a^((l)))
  $
  implying
  $
    U A U^dagger U ket(a^((l))) = a^((l)) U ket(a^((l)))
  $
  so we can write
  $
    [U A U^dagger] ket(b^((l))) = a^((l)) ket(b^((l)))
  $
]

Meaning $ket(b')$ are eigenkets of $U A U^dagger$ with the same eigenvalues as $A$. We see that $B$ and $U A U^dagger$ have simultaneous eigenkets. Typically they are the same operator $B = U A U^dagger$.

== Continuous observables
Many important observables are continuous. The state space in these cases becomes infinite.

We write
$ xi ket(xi') = xi' ket(xi') $
as before and replace everything by their continuous equivalent. As an example we have
$
  braket(xi', xi'') & = delta(xi'-xi'') \
              bb(1) & =integral dd(xi') ketbra(xi')
$

We assume the position eigenkets ${ket(x')}$ are complete meaning we can expand any $ket(alpha)$ as
$
  ket(alpha) = integral_(-oo)^oo dd(x') ket(x') braket(x', alpha)
$
We consider a detector making a measurement when a particle is near $x'$. This makes the state collapse
$
  ket(alpha) = integral_(-oo)^oo dd(x'') ket(x'') braket(x'', alpha) ->^"mm." integral_(x'-dd(x')\/2)^(x' + dd(x')\/2) dd(x'')ket(x'') braket(x'', alpha)
$
We assume $braket(x'', alpha) tilde$ constant within $x' plus.minus dd(x')\/2$. Then the probability that we detect the particle is $abs(braket(x', alpha))^2 dd(x')$. The probability of finding the particle anywhere is then
$
  integral_(-oo)^oo dd(x') abs(braket(x', alpha))^2 =^"normalized" 1 = braket(alpha)
$
The above generalizes to three dimensions. We assume the position eigenkets ${ket(bold(x)')}$ are complete meaning we can expand any $ket(alpha)$ as
$
  ket(alpha) = integral dd(x', 3) ket(bold(x)') braket(bold(x)', alpha)
$
where $ket(bold(x)')$ is a simultaneous eigenket of $x, y "and" z$. We implicitly assume $[x_i,x_j] = 0$.

== Translation operator
We seek an operator with the action $bold(x)' arrow bold(x)' +dd(bold(x)')$. We define the translation operator
$
  cal(J)(dd(bold(x)')) ket(bold(x)') equiv ket(bold(x)'+dd(bold(x)'))
$
meaning $ket(bold(x)')$ is not an eigenket of $cal(J) (dd(bold(x)'))$. The action on $ket(alpha)$ is
$
  cal(J) (dd(bold(x)'))ket(alpha) &= integral dd(x', 3) ket(bold(x)'+dd(bold(x)')) braket(bold(x)', alpha) =^(bold(x)' -> bold(x)'-dd(bold(x)')) integral dd(x', 3) ket(bold(x)') braket(bold(x)'-dd(bold(x)'), alpha)
$

We want $cal(J) (dd(bold(x))')$ to have certain properties. We require $cal(J)^dagger cal(J) = 1$ since $cal(J) (dd(bold(x)')) ket(alpha)$ should be normalized if $ket(alpha)$ is normalized. We also want
$
  cal(J) (dd(bold(x)'')) cal(J) (dd(bold(x)')) &= cal(J) (dd(bold(x)') + dd(bold(x)'')) \
  cal(J) (- dd(bold(x)')) &= cal(J)^(-1) (dd(bold(x)')) \
  lim_(dd(bold(x)') arrow 0) cal(J) (dd(bold(x)'))&=bb(1)
$
with the reason being obvious.

We claim
$
  cal(J) (dd(bold(x)')) = 1 - i bold(K) dot dd(bold(x)')
$
where $K_i$ are Hermitian. Checking this satisfies all the above is simple.

Now we want a relation between $bold(K)$ and $bold(x)$. We have
$
  [bold(x),cal(J) (dd(bold(x)'))] ket(bold(x)') = dd(bold(x)') ket(bold(x)' + dd(bold(x)')) tilde.eq dd(bold(x)') ket(bold(x)')
$
meaning
$
  bold(x) bold(K) dot dd(bold(x)') - bold(K) dot dd(bold(x)') bold(x) = i dd(bold(x)')
$
Let $dd(bold(x)')$ be along $hat(bold(x))_j$. Then acting with $hat(bold(x))_i dot (dots)$ gives
$
  [x_i, K_j] = i delta_(i j)
$
We define
$
  bold(K) equiv bold(p)/hbar
$
where $hbar$ has units of action. So we define the momentum operator $bold(p)$ as the generator of translations. We obtain
$
  cal(J) (dd(bold(x)')) = 1 - (i bold(p) dot dd(bold(x)'))/hbar
$
and
$ [x_i, p_j] = i hbar delta_(i j) $
A finite translation is given by
$
  cal(J) (dd(x', d: Delta) hat(bold(x))) &= lim_(N arrow oo) (1 - (i p_x dd(x', d: Delta))/(N hbar))^N eq exp(- (i p_x dd(x', d: Delta))/hbar)
$
where as usual
$
  exp(X) equiv 1 + X + X^2/2! + dots
$
We require $[cal(J) (dd(x'_i, d: Delta) hat(bold(x))_i), cal(J) (dd(x'_j, d: Delta) hat(bold(x))_j) ] = 0$ which gives
$
  [p_i,p_j] = 0
$
meaning the translation group is Abelian. Consider the simultanoeus eigenket $ket(bold(p)')$. We have
$
  cal(J) (dd(bold(x)')) ket(bold(p)') &= (1 - (i bold(p) dot dd(bold(x)'))/hbar) ket(bold(p)') = (1 - (i bold(p)' dot dd(bold(x)'))/hbar) ket(bold(p)')
$
meaning it is an eigenket of $cal(J) (dd(bold(x)'))$. Which is obvious since $[bold(p),cal(J) (dd(bold(x)'))] = 0$.

The main result of the above analysis are the canonical commutation relations
$
  [x_i,x_j] = [p_i,p_j] = 0";  " [x_i,p_j] = i hbar delta_(i j)
$
which are very important!

== The wavefunction
Consider the one-dimensional case where
$
  x ket(x') = x' ket(x')";  " braket(x'', x') = delta(x''-x')
$
and any state can be expanded as
$
  ket(alpha) = integral dd(x') ket(x') underbracket(braket(x', alpha), psi_alpha (x'))
$
we call $psi_alpha (x') equiv braket(x', alpha)$ the wavefunction.

Then
$
  braket(beta, alpha) = integral dd(x') psi_beta^* (x') psi_alpha (x')
$
this represents the probability that $ket(alpha)$ is found in $ket(beta)$.

Consider a discrete observable $A$ then
$
  braket(x', alpha) = sum_a' underbracket(braket(x', a'), u_a' (x')) braket(a', alpha) = sum_a' c_a' u_a' (x')
$
where $u_a' (x')$ is the eigenfunction of $A$ in the $x$-representation.

An element can be expanded as
$
  braket(beta, A, alpha) &= integral dd(x') integral dd(x'') braket(beta, x') braket(x', A, x'') braket(x'', alpha)
$
with $A = f(x)$ we find
$
  braket(beta, f(x), alpha) = integral dd(x') psi_beta^* (x') f(x') psi_alpha (x')
$
The action of $p$ on $ket(alpha)$ is
$
  (1 - (i p dd(x', d: Delta))/hbar) ket(alpha) &= integral dd(x') ket(x') braket(x' - dd(x', d: Delta), alpha) \
  &tilde.eq^"Taylor" integral dd(x') ket(x') (braket(x', alpha) - dd(x', d: Delta) pdv(, x') braket(x', alpha) )
$
by comparison we obtain
$
  p ket(alpha) = integral dd(x') ket(x') (- i hbar pdv(, x') braket(x', alpha))
$
acting with $bra(x'') dot (dots)$ gives the action of $p$ in the $x$-representation. And we find
$
  braket(beta, p, alpha) = integral dd(x') psi_beta^* (x') (- i hbar pdv(, x')) psi_alpha (x')
$
which looks familiar!

Likewise we have $p ket(p') = p' ket(p')$ and $ ket(alpha) = integral dd(p') ket(p') underbracket(braket(p', alpha), phi_alpha (p')) $
with $braket(p', alpha) = phi_alpha (p')$ being the momentum-space wavefunction.

We want the transformation operator from $x -> p$. We define the transformation function by $braket(x', p')$. Considering $braket(x', p, p')$ gives
$
  braket(x', p') = N exp((i p' x')/hbar)
$
for fixed $p'$ this is just a plane wave! To determine $N$ consider
$
  braket(x', x'') & = integral dd(p') braket(x', p') braket(p', x'') \
                  & = abs(N)^2 integral dd(p') exp((i p'(x'-x''))/hbar) \
                  & = 2 pi hbar abs(N)^2 delta(x'-x'') \
                  & =^! delta(x'-x'')
$
so $N = (2 pi hbar)^(-1\/2)$.

We can relate $psi_alpha (x')$ and $phi_alpha (p')$ using
$
  braket(x', alpha) & = integral dd(p') braket(x', p') braket(p', alpha) \
  braket(p', alpha) & = integral dd(x') braket(p', x') braket(x', alpha)
$
giving
$
  psi_alpha (x') &= 1/sqrt(2pi hbar) integral dd(p') exp((i p'x')/hbar) phi_alpha (p') \
  phi_alpha (p') &= 1/sqrt(2pi hbar) integral dd(x') exp((-i p' x')/hbar) psi_alpha (x')
$
this is just a Fourier transform!

The generalization to three dimensions should be clear.
/*
$
  bold(x) ket(bold(x)') = bold(x)' ket(bold(x)')",  " braket(bold(x)', bold(x)'') = delta^3 (bold(x)'-bold(x)'')
$
with completeness
$
  integral dd(x', 3) ketbra(bold(x)') = 1 => ket(alpha) = integral dd(x', 3) ket(bold(x)') braket(bold(x)', alpha)
$
and $braket(bold(x)', alpha)$ is identified with the wave function $psi_alpha (bold(x)')$---everything is the same for momentum. The momentum operator becomes
$
  braket(beta, bold(p), alpha) = integral dd(x', 3) psi_beta^* (bold(x)')(-i hbar grad') psi_alpha (bold(x)')
$
and the transformation function
$
  braket(bold(x)', bold(p)') = 1/(2 pi hbar)^(3\/2) exp((i bold(p)'dot bold(x)')/hbar)
$
so we find
$
  psi_alpha (bold(x)') &= 1/(2 pi hbar)^(3\/2) integral dd(p', 3) exp((i bold(p)'dot bold(x)')/hbar) phi_alpha (bold(p)') \
  phi_alpha (bold(p)') &= 1/(2pi hbar)^(3\/2) integral dd(x', 3) exp((- i bold(p)' dot bold(x)')/hbar) psi_alpha (bold(x)')
$
*/


