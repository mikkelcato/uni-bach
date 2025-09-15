//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *quantum physics*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on quantum physics taken during the SDU course. Based primarily on _Modern Quantum Mechanics_ by Sakurai and other notes.
  ],
)

= Fundamentals
== Stern-Gerlach
The Stern-Gerlach experiment devised by O. Stern with the help of W. Gerlach illustrates that we need to let go of classical mechanics. One might even say that a two-state system similar to that of the Stern-Gerlach experiment is the least classical system---so understanding these types of systems are of great importance.

In the experiment we heat silver atoms in an oven, allowing small amounts of atoms to escape. This beam is focused and then subject to an inhomogeneous magnetic field---we want to find the effect of this field. For our purposes a silver atom consists of a nucleus and 47 electrons, with 46 electron forming a spherically symmetric cloud---thus having no net angular momentum. Ignoring nuclear spin, the atom as a whole then has angular momentum, which is solely due to spin---instrinsic angular momentum of the 47th electron. The atom is heavy so it has a magnetic moment equal to the spin magnetic moment of the 47th electron, i.e. $bold(mu) prop bold(S)$---we essentially probe the spin of the 47th electron with the magnetic field using the silver atoms which are neutral, therefore all motion will be due to spin.

The interaction energy of the magnetic moment is just $- bold(mu) dot bold(B)$, so the force experienced is
$
  F_z = pdv(, z) (bold(mu)dot bold(B)) tilde.eq mu_z pdv(B_z, z)
$
Since the atom is heavy the classical idea of a trajectory can be applied. So depending on $mu_z$ or equivalently $S_z$ the deflection direction will change---so our beam will split. Given the atoms are randomly oriented in our oven, there is no preferred direction for $bold(mu)$. If the electron was a classical spinning objet, then we'd expect a continuous bunch of beams with $mu_z$ values evenly between $plus.minus abs(mu)$. What is observed however are two spots corresponding to an _up_ and a _down_ orientation, meaning we have two possible values of $S_z$. We now know this is $S_z = plus.minus hbar\/2$, importantly from the experiment we find that the electrons spin is quantized. By changing orientation we find a similar result for $S_x$ and $S_y$.

Going a step further and doing sequential experiments with differently oriented fields shows that we can't measure $S_x$ and $S_z$ simultaneously---essentially performing a measurement erases some information. This is completely unlike classical mechanics where we can measure and specify e.g. $L_x$ and $L_z$ of a spinning top easily.

#pagebreak()
== Bras and Kets
We consider a complex vector space with dimensionality depending on the given system we are looking at. The vector space in question is a Hilbert space $cal(H)$, which is a complete inner product space. Any physical state in quantum mechanics is represented by a state vector in this complex vector space. These we call kets and are denoted by $ket(alpha)$ following Dirac who introduced them. Everything we'd like to know about a given state is contained in the corresponding ket. We can add them $ket(alpha)+ket(beta)=ket(gamma)$ and multiply them by some complex number $c$ to give new kets. As physicists we postulate that $ket(alpha)$ and $c ket(alpha)$ both correspond to the same physical state---so we have an equivalence relation $ket(alpha) tilde c ket(alpha)$.

All observables can be represented by operators that act on the vector space. In general for some operator $A$ we have $A dot (ket(alpha)) = A ket(alpha)$, with $A ket(alpha)$ also being a ket. We define eigenkets of an operator $ ket(a'), ket(a''), dots $ as those with the property
$
  A ket(a') = a' ket(a')",  " A ket(a'') = a'' ket(a''), dots
$
where $a', a'', dots$ are numbers, namely the eigenvalues of the given operator. As an example take the case of spin-$1\/2$ systems where we have
$
  S_z ket(S_z";"plus.minus) = plus.minus hbar/2 ket(S_z";"plus.minus)
$
and similarly for $ket(S_x";"plus.minus)$ and $ket(S_y";"plus.minus)$, which are eigenkets of the operators $S_x$ and $S_y$ respectively.

Usually we are interested in an $N$-dimensional vector space spanned by $N$ eigenkets of some observable $A$, in this case any ket $ket(alpha)$ can be written as $ket(alpha) = sum_(a') c_(a') ket(a')$.

Now we introduce the bra space---which is the dual to the ket space. We postulate that every ket $ket(alpha)$ has a corresponding bra $bra(alpha)$ in the bra space. So we have a dual correspondence between the two spaces. We postulate that the bra dual to $c ket(alpha)$ is $c^* bra(alpha)$ so,
$
  c_alpha ket(alpha) + c_beta ket(beta) <->^"DC" c_alpha^* bra(alpha) + c_beta^* bra(beta)
$

We define the inner product or bracket as $braket(beta, alpha) = bra(beta) dot ket(alpha)$, which always gives a complex number. We postulate that $braket(beta, alpha)=braket(alpha, beta)^*$, from this we immediately find that $braket(alpha, alpha) >=0$, i.e. it is positive definite, with equality holding for the null ket. Two kets are orthogonal if
$
  braket(alpha, beta)=0 <=> braket(beta, alpha)=0
$
note that we can always normalize kets
$
  ket(tilde(alpha)) = 1/(sqrt(braket(alpha, alpha))) ket(alpha)
$
with $braket(tilde(alpha), tilde(alpha))=1$. The root is known as the norm of $ket(alpha)$.

We also have more general operators $X, Y, Z, dots$. These act on kets from the left $X ket(alpha)$ and on bras from the right $bra(beta) X$. Two operators are equal $X=Y$ if $X ket(alpha)= Y ket(beta)$. An operator is the null operator if $X ket(alpha) = 0$. Operators can be added, and addition is commutative and associative. Most operators we consider are also linear
$
  X (c_alpha ket(alpha) + c_beta ket(beta)) = c_alpha X ket(alpha) + c_beta X ket(beta)
$
we define the Hermitian adjoint as the dual operator $X ket(alpha) <->^"DC" bra(alpha)X^dagger$, if $X = X^dagger$ then $X$ is Hermitian.

Two operators can be multiplied, this is non-commutative in general, but it is associative
$
  X(Y ket(alpha)) = (X Y)ket(alpha)=X Y ket(alpha)
$
similarly for bras. Notice also $(X Y)^dagger = Y^dagger X^dagger$, since
$
  X Y ket(alpha) = X(Y ket(alpha)) <->^"DC" (bra(alpha)Y^dagger)X^dagger = bra(alpha) Y^dagger X^dagger
$
We can also consider the outer product $ketbra(beta, alpha)$ this is an operator.

Multiplication amongst operators is associative, but we postulate that it holds among everything, so kets, bras and operators are all associative under multiplication---Dirac calls this the associative axiom of multiplication. So we can say that $(ketbra(beta, alpha)) dot ket(gamma) = ket(beta) (braket(alpha, gamma))$, where the inner product is just a number, so acting with the outer product gives us another ket---it is an operator. Notice that if $X = ketbra(beta, alpha)$ then $X^dagger = ketbra(alpha, beta)$. As a point of notation we have $(bra(beta))dot (X ket(alpha)) = (bra(beta)X) dot (ket(alpha)) = braket(beta, X, alpha)$ and
$
  braket(beta, X, alpha) = bra(beta) dot (X ket(alpha)) = {(bra(alpha)X^dagger)dot ket(beta)}^* = braket(alpha, X^dagger, beta)^*
$
#pagebreak()
== Base Kets
We consider the eigenkets and eigenvalues of some Hermitian operator $A$.
#thm[The eigenvalues of a Hermitian operator $A$ are real, and the eigenkets of $A$ corresponding to different eigenvalues are orthogonal.]
#proof[
  Recall $A ket(a') = a' ket(a')$ and since $A$ is Hermitian we also have $bra(a'')A = a''^* bra(a'')$. Combining these gives $(a'-a''^*) braket(a'', a') = 0$. Now we let $a''=a'$, in this case we require $a'=a'^*$, so they are real. Assuming they are different $a'eq.not a''$, then $a'-a''^*=a'-a'' eq.not 0$ by assumption. So we require $braket(a'', a')=0$ meaning they are orthogonal.
]

Physically we expect observables to have real eigenvalues---this theorem guarantees that eigenvalues are real if the operator is Hermitian. Usually we normalize $ket(a')$ so that ${ket(a')}$ form an orthonormal set $braket(a'', a')=delta_(a'' a')$. By construction this set is complete, since we asserted at the very beginning that the ket space is spanned by the eigenkets of $A$. Given some arbritrary $ket(alpha)$ we can write
$
  ket(alpha) = sum_(a') c_(a') ket(a')
$
multiplying on the left by $bra(a'')$ we find $c_(a') = braket(a', alpha)$ so
$
  ket(alpha) = sum_(a') ket(a') braket(a', alpha)
$
since $ket(alpha)$ is arbitrary we get the completeness relation
$
  sum_(a') ketbra(a', a') = 1
$
which is very useful. Consider
$
  braket(alpha, alpha) = bra(alpha) dot (sum_(a') ketbra(a', a')) dot ket(alpha) = sum_(a') abs(braket(a', alpha))^2
$
given $ket(alpha)$ is normalized then
$
  sum_(a') abs(c_(a'))^2 = sum_(a') abs(braket(a', alpha))^2 = 1
$
Consider again $ketbra(a', a')$ and let it act on $ket(alpha)$
$
  (ketbra(a', a'))ket(alpha) = c_(a') ket(a')
$
so it pick outs the part of $ket(alpha)$ parallel to $ket(a')$. For this reason we call it the projection operator along the base ket $ket(a')$
$
  Lambda_(a') equiv ketbra(a', a')
$
and the completeness relation becomes
$
  sum_(a') Lambda_(a') = 1
$

=== Matrices
We can represent operators by matrices; using the completeness relation twice we find
$
  X = sum_(a'') sum_(a') ket(a'') braket(a'', X, a') bra(a')
$
so there are $N^2$ numbers of the form $braket(a'', X, a')$. We can arrange these elements in a matrix and represent the operator as
$
  X eq^dot mat(braket(a^((1)), X, a^((1))), braket(a^((1)), X, a^((2))), dots; braket(a^((2)), X, a^((1))), braket(a^((2)), X, a^((2))), dots; dots.v, dots.v, dots.down)
$
we can write $braket(a'', X, a') = braket(a', X^dagger, a'')^*$, so the Hermitian conjugate has been related to the complex conjugate transpose.
/*
This is also in agreement with the usual way we do matrix multiplication $Z = X Y$ becomes
$
  braket(a'', Z, a') = braket(a'', X Y, a') = sum_(a''') braket(a'', X, a''')braket(a''', Y, a')
$
Now consider $ket(gamma)=X ket(alpha)$, then
$
  braket(a', gamma) = braket(a', X, alpha)=sum_(a'')braket(a', X, a'')braket(a'', alpha)
$
this can be seen as multiplying a square matrix by a column matrix, if we write
$
  ket(alpha) eq^dot vec(braket(a^((1)), alpha), dots.v)",  " ket(gamma) eq^dot vec(braket(a^((1)), gamma), dots.v)
$
similarly for $bra(gamma)=bra(alpha) X$ with $braket(gamma, a') = sum_(a'') braket(alpha, a'') braket(a'', X, a')$, so a bra becomes a row matrix
$
  bra(gamma) eq^dot vecrow(braket(gamma, a^((1))), dots) = vecrow(braket(a^((1)), gamma)^*, dots)
$
We can also write the inner product as
$
  braket(beta, alpha) = sum_(a') braket(beta, a')braket(a', alpha) = vecrow(braket(a^((1)), beta)^*, dots) vec(braket(a^((1)), alpha), dots.v)
$
And the outer product becomes
$
  ketbra(beta, alpha) eq^dot mat(braket(a^((1)), beta)braket(a^((1)), alpha)^*, dots; dots.v, dots.down)
$
*/

For an observable this representation becomes very simple---if the eigenkets are used as the base kets,
$
  A & = sum_(a'') sum_(a') ket(a'') braket(a'', A, a') bra(a') \
    & = sum_(a') a' ketbra(a', a') = sum_(a') a' Lambda_(a')
$
since $braket(a'', A, a') = braket(a', A, a') delta_(a'' a') = a' delta_(a'' a')$.

As an example consider a spin-$1\/2$ system. The base kets are $ket(S_z";"plus.minus) = ket(plus.minus)$. The simplest operator is the identity
$
  1 = ketbra(+, +) + ketbra(-, -)
$
per the previous we can write
$
  S_z = hbar/2 (ketbra(+, +)-ketbra(-, -))
$
immediately
$
  S_z ket(plus.minus) = plus.minus hbar/2 ket(plus.minus)
$
by orthonormality of $ket(plus.minus)$. We can also look at two other operators
$
  S_+ equiv hbar ketbra(plus, minus)",  " S_- equiv hbar ketbra(-, +)
$
these are ladder operators and their action is obvious. From these we can construct $S_x$ and $S_y$ by $S_(plus.minus) = S_x plus.minus i S_y$, which we show later. Now we can construct the matrix representation
$
  ket(+) eq^dot vec(1, 0)",  " ket(-) eq^dot vec(0, 1)
$
$
  S_z eq^dot hbar/2 mat(1, 0; 0, -1)",  " S_+ eq^dot hbar mat(0, 1; 0, 0)",  " S_- eq^dot hbar mat(0, 0; 1, 0)
$

#pagebreak()
== Measurement
Before any measurement of observable $A$, the system is assumed to be represented by some linear combination
$
  ket(alpha) = sum_(a') c_(a') ket(a') = sum_(a') ket(a') braket(a', alpha)
$
when a measurement is made the system collapses into one of the eigenstates of $A$
$ ket(alpha) -->^"m.m" ket(a') $
usually a measurement would therefore change the state, but if the system is already in some eigenstate then $ket(a') -->^"m.m" ket(a')$. When the measurement is made and $ket(alpha) --> ket(a')$, then $A$ is measured to be $a'$---so the result always yields an eigenvalue of the observable.

Which value will be taken is not known, but we postulate that the probability for $a'$ is $abs(braket(a', alpha))^2$, given $ket(alpha)$ is normalized. Empirically to determine the probability we need to take measurements of identically prepared systems---an ensemble. This probabilistic interpretation is taken as an axion. One consequence is that repeated measurement of the same observable will give the same value, since $braket(a', a') = 1$. Similarly there is no chance for some system to go from $ket(a') --> ket(a'')$ since $braket(a'', a')=0$---so orthogonal kets correspond to mutually exclusive options. Further probability must be nonnegative, and the total probability should be unity, both of these are satisfied.

We define the expectation value of $A$ with respect to $ket(alpha)$ as
$
  expval(A) equiv expval(A, alpha)
$
sometimes we write $expval(A)_alpha$, this definition agrees with what we'd expect
$
  expval(A) = sum_(a') sum_(a'') braket(alpha, a'') braket(a'', A, a') braket(a', alpha) = sum_(a') a' abs(braket(a', alpha))^2
$

=== Spin-$1\/2$
We visit the spin-$1\/2$ system again$dots$
/*
In the Stern-Gerlach sending a $S_x +$ beam into $"SG"hat(z)$ splits it into two equal parts. So
$
  abs(braket(+, S_x";"+))=abs(braket(-, S_x";"+))=1/sqrt(2)
$
so we can write
$
  ket(S_x";"+) = 1/sqrt(2) ket(+) + 1/sqrt(2) e^(i delta_1) ket(-)
$
where we use that the overall phase common to $ket(plus.minus)$ doesn't matter. This also gives us
$
  ket(S_x";"-) = 1/sqrt(2) ket(+) - 1/sqrt(2) e^(i delta_1) ket(-)
$
since it must be orthogonal. Then we construct the operator
$
  S_x & = hbar/2 [(ketbra(S_x";"+, S_x";"+))-(ketbra(S_x";"-, S_x";"-))] \
      & = hbar/2 [e^(- i delta_1) (ketbra(+, -)) + e^(i delta_1)(ketbra(-, +)) ]
$
similarly gives
$
  ket(S_y";"plus.minus) &= 1/sqrt(2) ket(+) plus.minus 1/sqrt(2) e^(i delta_2) ket(-) \
  S_y &= hbar/2 [e^(-i delta_2) (ketbra(+, -)) + e^(i delta_2)(ketbra(-, +))]
$
lastly we use $"SG"hat(x) arrow "SG"hat(y)$ to get
$
  abs(braket(S_y";"plus.minus, S_x";"+)) = abs(braket(S_y";"plus.minus, S_x";"-))=1/sqrt(2)
$
plugging everything in gives
$
  1/2 abs(1 plus.minus e^(i (delta_1-delta_2))) = 1/sqrt(2) => delta_2-delta_1 = plus.minus pi/2
$
it is simplest to pick $delta_1=0$ and then $delta_2 = pi\/2$ turns out to be correct.
*/

To summarize
$
  ket(S_x";"plus.minus) &= 1/sqrt(2) ket(+) plus.minus 1/sqrt(2) ket(-)",  " ket(S_y";"plus.minus) = 1/sqrt(2) ket(+) plus.minus i/sqrt(2) ket(-) \
$
$
  S_x = hbar/2 [(ketbra(+, -))+(ketbra(-, +))]",  " S_y = hbar/2 [-i (ketbra(+, -))+i(ketbra(-, +))]
$
we notice
$
  S_plus.minus = S_x plus.minus i S_y
$
these satisfy
$
  [S_i,S_j] = i epsilon_(i j k) hbar S_k
$
and
$
  {S_i, S_j} = 1/2 hbar^2 delta_(i j)
$
with the usual definitions
$
  [A,B] = A B - B A",  " {A, B} = A B + B A
$
we can also define $bold(S) dot bold(S) = bold(S)^2 equiv S_x^2 + S_y^2 + S_z^2$. It follows immediately from the anticommutation relation that
$
  bold(S)^2 = 3/4 hbar^2
$
so $[bold(S)^2, S_i]=0$.

=== Compatible and incompatible observables
Observables $A$ and $B$ are compatible if their operators commute $[A,B]=0$ and incompatible if they don't $[A,B]eq.not 0$.

We consider compatible $A$ and $B$ first, and we assume the ket space is spanned by the eigenkets of $A$. We want to know how the eigenkets of $A$ and of $B$ are related. We assume the eigenvalues of $A$ are non-degenerate, since degeneracy will pose a problem.

#thm[
  Suppose $A$ and $B$ are compatible, and the eigenvalues of $A$ are non-degenerate. Then all matrix elements $braket(a'', B, a')$ are diagonal.
]
#proof[
  By definition notice
  $
    braket(a'', [A,B], a') = (a''-a') braket(a'', B, a') = 0
  $
  so $braket(a'', B, a')$ must vanish unless $a'' = a'$.
]

We can write $braket(a'', B, a') = delta_(a'' a') braket(a', B, a')$, and both $A$ and $B$ can be represented by diagonal matrices with the same base kets. With this we can write
$
  B = sum_(a') ket(a'') braket(a'', B, a'') bra(a'')
$
and letting this act on some eigenket of $A$
$
  B ket(a') = sum_(a'') ket(a'') braket(a'', B, a'') braket(a'', a') = (braket(a', B, a')) ket(a')
$
so $ket(a')$ is a simultaneous eigenket of $A$ and $B$---we sometimes write $ket(a'b')$, or $ket(K')$ where $K'$ is a collective index. This is also valid for $n$-fold degeneracy
$
  A ket(a'^((i))) = a' ket(a'^((i)))"  for " i=1,2,dots,n
$
where $ket(a'^((i)))$ are $n$ mutually orthonormal eigenkets of $A$ with the same eigenvalue $a'$.

We can generalize this by considering more observables which are mutually compatible
$
  [A,B] = [B,C] = [A,C] = dots = 0
$
assume we have a maximal set of compatible observables. The eigenvalues of individual operators $A,B,C,dots$ may have degeneracies, but if we specify a combination $(a',b',c',dots)$, then the simultaneous eigenket of $A,B,C,dots$ is unique. The orthonormality condition states
$
  braket(K'', K') = delta_(K'K'') = delta_(a a') delta_(b b') dots
$
where $K'$ stands for $(a',b',c',dots)$. Similarly completeness can be written as
$
  sum_(K') ketbra(K', K') = sum_(a') sum_(b') sum_(c') dots ketbra(a' b' c' dots) = 1
$

The nice thing about compatible observables is that they don't interfere with each other upon measurement since they share eigenkets. Consider two compatible observables $A$ and $B$, after measuring $A$ giving some $a'$ our system collapses into some linear combination
$
  sum_i^n c_(a')^((i)) ket(a'b^((i)))
$
where $n$ is the degeneracy. Measuring $B$ would then pick one of the terms, but measuring $A$ again would still give $a'$. So they don't interfere, which is why we call them compatible.

Incompatible observables are more annoying. These don't have a complete set of simultaneous eigenkets. To show this we can proceed by contradiction, if it were true there would exist a set of simultaneous eigenkets satisfying
$
  A ket(a' b') = a' ket(a' b')",  " B ket(a' b') = b' ket(a' b')
$
so
$
  A B ket(a' b') = A b' ket(a' b') = a' b' ket(a' b')",  " B A ket(a' b') = B a' ket(a' b') = a' b' ket(a' b')
$
or $A B ket(a' b') = B A ket(a' b')$ implying that $[A, B]=0$ which is a contradiction. However, this can be true in a ket subspace, even if $A$ and $B$ are incompatible.
/*
Consider a measurement of $C$ taking $ket(a') arrow^"C" ket(c')$, and in the case where we perform a selective measurement first $ket(a') arrow^"B" ket(b') arrow^"C" ket(c')$. Since probabilites are multiplicative we have in the latter case probability $abs(braket(c', b'))^2 abs(braket(b', a'))^2$ to measure $ket(c')$. To get the total probability we sum over all $b'$ since this would seem to remove the filter
$
  sum_(b') abs(braket(c', b'))^2 abs(braket(b', a'))^2 = sum_(b') braket(c', b') braket(b', a') braket(a', b') braket(b', c')
$
the total probability in the first case is just
$
  abs(braket(c', a'))^2 = abs(sum_(b') braket(c', b')braket(b', a'))^2 = sum_(b') sum_(b'') braket(c', b')braket(b', a') braket(a', b'')braket(b'', c')
$
these are very different. Which is cool since in each case the $ket(a')$ can be regared as built from $ket(b')$, $ket(a') = sum_(b') ket(b') braket(b', a')$. This shows how taking a measurement of $B$ actually changes things---in the first we actually record which eigenvalues were realized, in the second we just imagine that $ket(a')$ is built from the $ket(b')$. So actually writing down the different $b'$ routes changes the result. These expressions are however equal for compatible observables, i.e. in the case where $[A,B] = [B,C] = 0$. To see this consider
$
  sum_(b') abs(braket(c', b'))^2 abs(braket(b', a'))^2
$
since both $A$ and $C$ are compatible with $B$ every $ket(a')$ and $ket(c')$ is one of the $ket(b')$, so $braket(c', b') = delta_(c' b')$ and $braket(b', a') = delta_(b' a')$
$
  sum_(b') abs(delta_(c' b'))^2 abs(delta_(b' a'))^2
$
since $sum_(b') delta_(c' b') delta_(b' a') = delta_(c' a')$ we find
$
  abs(delta_(c' a'))^2 = abs(braket(c', a'))^2
$
*/

=== Uncertainty
We define an operator $Delta A equiv A-expval(A)$, this is used to find the dispersion
$ expval((Delta A)^2) = expval(A^2)-expval(A)^2 $
this characterizes fuzziness of an observable---for eigenstates of $A$ this clearly vanishes.

#lemma[Schwarz inequality][
  $ braket(alpha, alpha) braket(beta, beta) >= abs(braket(alpha, beta))^2 $
]
#proof[
  Note that
  $
    (bra(alpha)+lambda^* bra(beta)) dot (ket(alpha)+lambda ket(beta)) >= 0
  $
  this holds for any $lambda in CC$. So it holds for
  $
    lambda = - (braket(beta, alpha))/braket(beta, beta) => braket(alpha, alpha)braket(beta, beta) - abs(braket(alpha, beta))^2 >= 0
  $
]

#lemma[
  The expectation value of a Hermitian operator is real.
]
#proof[
  $
    braket(a, A, a) = braket(a, A^dagger, a)^* = braket(a, A, a)^*
  $
]

#lemma[
  The expectation value of an anti-Hermitian operator $C = - C^dagger$ is purely imaginary.
]
#proof[
  $
    braket(c, C, c) = braket(c, C^dagger, c)^* = - braket(c, C, c)^*
  $
]

Now we'll prove the uncertainty principle
#thm[
  Let $A$ and $B$ be observables, then for any state
  $
    expval((Delta A)^2)expval((Delta B)^2) >= 1/4 abs(expval([A,B]))^2
  $
]
#proof[
  Let
  $
    ket(alpha) = Delta A ket(dot)",  " ket(beta) = Delta B ket(dot)
  $
  where $ket(dot)$ is any state. From Schwarz we get
  $
    expval((Delta A)^2)expval((Delta B)^2) >= abs(expval(Delta A Delta B))^2
  $
  where we use that $Delta A$ and $Delta B$ are Hermitian. For the right-hand side note
  $
    Delta A Delta B = 1/2 [Delta A, Delta B] + 1/2 {Delta A, Delta B}
  $
  with $[Delta A, Delta B] = [A, B]$ since expectation values are numbers, and it is anti-Hermitian
  $
    ([A,B])^dagger = B A - A B = -[A,B]
  $
  as opposed to the anticommutator ${Delta A, Delta B}$ which is Hermitian. So in
  $
    expval(Delta A Delta B) = 1/2 expval([A,B]) + 1/2 expval({Delta A, Delta B})
  $
  the first term is purely imaginary and the second term is purely real---from the lemmas. So the right-hand side becomes
  $
    abs(expval(Delta A Delta B))^2 = 1/4 abs(expval([A,B]))^2 + 1/4 abs(expval({Delta A, Delta B}))^2
  $
  and we are done since removing the second term can only make the inequality stronger.
]

== Changing basis
Given two incompatible observables $A$ and $B$ the ket space can be spanned by either ${ket(a')}$ or ${ket(b')}$. We want to figure out how these are related, and how we can change our basis or representation. We'll do this through a transformation operator.

#thm[
  Given two sets of base kets, both satisfying orthonormality and completeness, there exists a unitary operator $U$ such that
  $
    ket(b^((1))) = U ket(a^((1))), ket(b^((2))) = U ket(a^((2))),dots
  $
]

A unitary operator satisfies $U^dagger U = U U^dagger = 1$.

#proof[
  We proceed by construction. We claim
  $
    U = sum_k ketbra(b^((k)), a^((k)))
  $
  is good enough. Clearly
  $
    U ket(a^((l)))=ket(b^((l)))
  $
  since ${ket(a')}$ satisfies orthonormality. And it is unitary
  $
    U^dagger U = sum_k sum_l ket(a^((l))) braket(b^((l)), b^((k))) bra(a^((k))) = sum_k ketbra(a^((k)))=1
  $
  where we use orthonormality of ${ket(b')}$ and completeness of ${ket(a')}$.
]

The matrix elements---the transformation matrix---are easy to find
$
  braket(a^((k)), U, a^((l))) = braket(a^((k)), b^((l)))
$
so the matrix elements are built from the old bras and new kets.

Consider some expansion
$
  ket(alpha) = sum_l ket(a^((l))) braket(a^((l)), alpha)
$
then it's easy to get the new expansion coefficients,
$
  braket(b^((k)), alpha) = sum_l braket(b^((k)), a^((l))) braket(a^((l)), alpha) = sum_l braket(a^((k)), U^dagger, a^((l))) braket(a^((l)), alpha)
$
similarly
$
  braket(b^((k)), X, b^((l))) &= sum_m sum_n braket(b^((k)), a^((m)))braket(a^((m)), X, a^((n)))braket(a^((n)), b^((l))) \
  &= sum_m sum_n braket(a^((k)), U^dagger, a^((m))) braket(a^((m)), X, a^((n)))braket(a^((n)), U, a^((l)))
$
The trace of an operator is defined as
$
  tr(X) & = sum_a' braket(a', X, a')
$
/*
and as shown it is representation independent. Likewise we can prove other things
$
  tr(X Y) & = sum_a' braket(a', X Y, a') \
          & = sum_b' sum_a' braket(a', X, b') braket(b', Y, a') \
          & = sum_b' sum_a' braket(b', Y, a') braket(a', X, b') \
          & = sum_b' braket(b', Y X, b')
$
and
$
  tr(U^dagger X U) & = tr(U^dagger Y) \
                   & = tr(Y U^dagger) \
                   & = tr(X U U^dagger) \
                   & = tr(X)
$
and
$
  tr(ketbra(a', a'')) & = sum_b' braket(b', a')braket(a'', b') \
                      & = sum_b' braket(a'', b') braket(b', a') \
                      & = braket(a'', a') \
                      & = delta_(a' a'')
$
similarly
*/
and it is ndependent of representation (see exercises for more), a nice property is
$
  tr(ketbra(b', a')) & = braket(a', b')
$
We are interested in finding $b'$ and $ket(b')$ with
$
  B ket(b') = b' ket(b')
$
we can write this as
$
  sum_a' braket(a'', B, a') braket(a', b') = b' braket(a'', b')
$
if $ket(b')$ is the $l$'th eigenket of $B$ then this can be written as $ sum_j B_(i j) C_j^((l)) = b^((l)) C_i^((l)) $
where $B_(i j) =braket(a^((i)), B, a^((j)))$ and $C_k^((l)) = braket(a^((k)), b^((l)))$---this could be written as explicit matrices. Nontrivial solutions are given by
$
  det(B- lambda II) = 0
$
with the $N$ roots being various $b^((l))$. Knowing these we can find $C_k^((l))$ up to a constant---as should be clear $C_k^((l))$ are the elements of the transformation matrix.

#thm[
  Consider to orthonormal bases ${ket(a')}$ and ${ket(b')}$ connected by $U$. Then we can construct a unitary transform of $A$: $U A U^dagger$, then $A$ and $U A U^dagger$ are unitary equivalent observables.
]
#proof[
  We can write
  $
    A ket(a^((l))) = a^((l)) ket(a^((l)))
  $
  which implies
  $
    U A U^dagger U ket(a^((l))) = a^((l)) U ket(a^((l)))
  $
  which we can rewrite
  $
    (U A U^dagger) ket(b^((l))) = a^((l)) ket(b^((l)))
  $
]

So the eigenkets $ket(b')$ of $U A U^dagger$ have the same eigenvalues as $A$---so they have identical spectra. This also implies that $B$ and $U A U^dagger$ have simultaneous eigenkets---typically they are the same operator $B = U A U^dagger$.

#pagebreak()
== $x$, $p$ and translation.
So far we have looked at observables with discrete spectra, but many observables have continuous spectra. In these cases the ket space becomes infinite---many results carry over.

We have $xi ket(xi') = xi' ket(xi')$ as before and most other things get replaced by their continuous analog---sums become integrals, and Kronecker deltas become Dirac delta functions:
$
  braket(xi', xi'') & = delta(xi'-xi'') => braket(xi'', xi, xi') = xi' delta(xi''-xi') \
  1 & =integral dd(xi') ketbra(xi') => ket(alpha) =integral dd(xi') ket(xi') braket(xi', alpha) => braket(beta, alpha) = integral dd(xi') braket(beta, xi') braket(xi', alpha) \
  1 & = integral dd(xi') abs(braket(xi', alpha))^2
$

=== $x$ eigenkets
The eigenkets $ket(x')$ satisfy $x ket(x') = x' ket(x')$, and are postulated to form a complete set. So for any $ket(alpha)$
$
  ket(alpha) = integral_(-oo)^oo dd(x') ket(x') braket(x', alpha)
$
now suppose we use a detector that only clicks when the particle is near $x'$, then
$
  ket(alpha) = integral_(-oo)^oo dd(x'') ket(x'') braket(x'', alpha) ->^"mm." integral_(x'-Delta\/2)^(x' + Delta\/2) dd(x'')ket(x'') braket(x'', alpha)
$
assuming that $braket(x'', alpha)$ is unaffected within the interval, the probability that our detector clicks is $abs(braket(x', alpha))^2 dd(x')$, where we've let $Delta arrow dd(x')$. The probability of finding the particle somewhere is
$
  integral_(-oo)^oo dd(x') abs(braket(x', alpha))^2
$
which is unity given $ket(alpha)$ is normalized
$
  braket(alpha, alpha) = integral_(-oo)^oo dd(x') braket(alpha, x') braket(x', alpha) = 1
$
This generalizes to three dimensions---assuming that the position eigenkets $ket(bold(x)')$ are complete
$
  ket(alpha) = integral dd(x', 3) ket(bold(x)') braket(bold(x)', alpha)
$
where $ket(bold(x)')$ is a simultaneous eigenket of $x, y$ and $z$
$
  ket(bold(x)') equiv ket(x'","y'","z')" with" x ket(bold(x)') = x' ket(bold(x)')" etc."
$
doing this we assume that $x, y$ and $z$ are compatible:
$
  [x_i,x_j] = 0
$

=== Translation and $p$
We want an operator which can take a state from $bold(x)' arrow bold(x)' +dd(bold(x)')$. This is the translation operator and it's defined by
$
  cal(J)(dd(bold(x)')) ket(bold(x)') = ket(bold(x)'+dd(bold(x)'))
$
clearly $ket(bold(x)')$ is not an eigenket, now consider the action on some arbitrary $ket(alpha)$
$
  cal(J) (dd(bold(x)'))ket(alpha) &= integral dd(x', 3) ket(bold(x)'+dd(bold(x)')) braket(bold(x)', alpha) = integral dd(x', 3) ket(bold(x)') braket(bold(x)'-dd(bold(x)'), alpha)
$
where the shift is allowed since we integrate over all of space. So the translated state is obtained by subbing $bold(x)' -> bold(x)' - dd(bold(x)')$.

We want $cal(J) (dd(bold(x))')$ to have certain properties. It should be unitary, since if $ket(alpha)$ is normalized then $cal(J) (dd(bold(x)')) ket(alpha)$ should also be normalized,
$
  braket(alpha) = braket(alpha, cal(J)^dagger (dd(bold(x)')) cal(J) (dd(bold(x)')), alpha) => cal(J)^dagger (dd(bold(x)')) cal(J) (dd(bold(x)')) = 1
$
this is what being unitary means. Some obvious properties:
$
  cal(J) (dd(bold(x)'')) cal(J) (dd(bold(x)')) &= cal(J) (dd(bold(x)') + dd(bold(x)'')) \
  cal(J) (- dd(bold(x)')) &= cal(J)^(-1) (dd(bold(x)')) \
  lim_(dd(bold(x)') arrow 0) cal(J) (dd(bold(x)'))&=1
$
but what does $cal(J) (dd(bold(x)'))$ look like? If we define
$
  cal(J) (dd(bold(x)')) = 1 - i bold(K) dot dd(bold(x)')
$
where $K_i$ are Hermitian operators---then everything is satisfied. Which is easily checked.

Now we want a relation between $bold(K)$ and $bold(x)$---notice
$
  bold(x) cal(J) (dd(bold(x)')) ket(bold(x)') &= bold(x) ket(bold(x)' + dd(bold(x)')) = (bold(x)' + dd(bold(x)'))ket(bold(x)' + dd(bold(x)')) \
  cal(J) (dd(bold(x)')) bold(x) ket(bold(x)') &= bold(x)' cal(J) (dd(bold(x)')) ket(bold(x)') = bold(x)' ket(bold(x)' + dd(bold(x)'))
$
combining these
$
  [bold(x),cal(J) (dd(bold(x)'))] ket(bold(x)') = dd(bold(x)') ket(bold(x)' + dd(bold(x)')) tilde.eq dd(bold(x)') ket(bold(x)')
$
so we have found the identity $ [bold(x),cal(J) (dd(bold(x)'))] = dd(bold(x)') =>
-i bold(x) bold(K) dot dd(bold(x)') + i bold(K) dot dd(bold(x)') bold(x) = dd(bold(x)') $
here $dd(bold(x)')$ is a number multiplied by the identity operator in the ket space spanned by $ket(bold(x)')$. Letting $dd(bold(x)')$ be along $hat(bold(x))_j$ and forming the scalar product with $hat(bold(x))_i$ we obtain
$
  [x_i, K_j] = i delta_(i j)
$
but what is $bold(K)$ physically? Borrowing from classical generating functions it seems like $bold(K)$ might be related to $bold(p)$. For this to be the case we need
$
  bold(K) = bold(p)/"constant with unit of action"
$
The constant turns out to be $hbar$---consider de Broglie's relation $2 pi\/lambda = p\/hbar$---and the $K$ operator appears to correspond to the wave number $k$. So we find
$
  cal(J) (dd(bold(x)')) = 1 - i bold(p) dot dd(bold(x)')\/hbar
$
where $bold(p)$ is the momentum operator, and we obtain: $ [x_i, p_j] = i hbar delta_(i j) $

What happens for a finite translation $dd(x', d: Delta)$? We can compound $N -> oo$ translations with displacement $dd(x', d: Delta)\/N$ and obtain
$
  cal(J) (dd(x', d: Delta) hat(bold(x))) &= lim_(N arrow oo) (1 - (i p_x dd(x', d: Delta))/(N hbar))^N = exp(- (i p_x dd(x', d: Delta))/hbar)
$
where we define
$
  exp(X) equiv 1 + X + X^2/2! + dots
$
a fundamental property of translations is that successive translations in different directions commute $[cal(J) (dd(y', d: Delta) hat(bold(y))), cal(J) (dd(x', d: Delta) hat(bold(x))) ] = 0$. Using this we obtain
$
  [cal(J) (dd(y', d: Delta) hat(bold(y))), cal(J) (dd(x', d: Delta) hat(bold(x))) ] &= [(1-(i p_y dd(y', d: Delta))/hbar + dots), (1-(i p_x dd(x', d: Delta))/hbar + dots)] \
  &tilde.eq - ((dd(x', d: Delta)) (dd(y', d: Delta)) [p_y,p_x])/hbar^2
$
the displacements are arbitrary so:
$
  [p_x,p_y]=0 => [p_i,p_j] = 0
$
this is a direct result of the fact that translations in different directions commute---so the translation group in three dimensions is Abelian. This implies we can find a simultaneous eigenket of $p_x, p_y$ and $p_z$:
$
  ket(bold(p)') equiv ket(p'_x","p'_y","p'_z)",  " p_x ket(bold(p)') = p'_x ket(bold(p)') "etc."
$
the effect of the translation operator on this state is
$
  cal(J) (dd(bold(x)')) ket(bold(p)') &= (1 - (i bold(p) dot dd(bold(x)'))/hbar) ket(bold(p)') = (1 - (i bold(p)' dot dd(bold(x)'))/hbar) ket(bold(p)')
$
so $ket(bold(p)')$ is an eigenket of $cal(J) (dd(bold(x)'))$ which is obvious since $[bold(p),cal(J) (dd(bold(x)'))] = 0$.

What we have found are the canonical commutation relations:
$
  [x_i,x_j] = 0", " [p_i,p_j] = 0", " [x_i,p_j] = i hbar delta_(i j)
$
following Dirac one might notice that these are very similar to Poisson brackets
$
  [dot,dot]_"class" -> [dot,dot]/(i hbar)
$
but we won't care.

#pagebreak()
== Wave-functions
=== Position space
In one dimension we have
$
  x ket(x') = x' ket(x') " with " braket(x'', x') = delta(x''-x')
$
any state can be written as
$
  ket(alpha) = integral dd(x') ket(x') braket(x', alpha)
$
with $abs(braket(x', alpha))^2 dd(x')$ being the probability to find our particle near $x'$. We refer to $braket(x', alpha)$ as the wave function $psi_alpha (x')$.

Consider the inner product
$
  braket(beta, alpha) = integral dd(x') braket(beta, x') braket(x', alpha) = integral dd(x') psi_beta^* (x') psi_alpha (x')
$
so $braket(beta, alpha)$$<-->$overlap---in general $braket(beta, alpha)$ represents the probability for $ket(alpha)$ to be found in $ket(beta)$.

We can also look at
$
  ket(alpha) = sum_a' ket(a') braket(a', alpha) => braket(x', alpha) = sum_a' braket(x', a') braket(a', alpha)
$
this is typically written $psi_alpha (x') = sum_a' c_a' u_a' (x')$ where $u_a' (x') = braket(x', a')$ is an eigenfunction of an operator $A$---in the $x$-representation---with eigenvalue $a'$. Similarly
$
  braket(beta, A, alpha) &= integral dd(x') integral dd(x'') braket(beta, x') braket(x', A, x'') braket(x'', alpha)
$
if $A = f(x)$
$
  braket(beta, f(x), alpha) = integral dd(x') psi_beta^* (x') f(x') psi_alpha (x')
$

We want $p$ in this basis, starting with $p$ as the generator of translations
$
  (1 - (i p dd(x', d: Delta))/hbar) ket(alpha) &= integral dd(x') cal(J) (dd(x', d: Delta)) ket(x')) braket(x', alpha) = integral dd(x') ket(x') braket(x' - dd(x', d: Delta), alpha) \
  &= integral dd(x') ket(x') (braket(x', alpha) - dd(x', d: Delta) pdv(, x') braket(x', alpha) )
$
comparison gives
$
  p ket(alpha) = integral dd(x') ket(x') (- i hbar pdv(, x') braket(x', alpha)) => braket(x', p, alpha) = - i hbar pdv(, x') braket(x', alpha)
$
so this is the action of $p$ in the $x$-representation, /* the matrix element is
                                                        $
                                                          braket(x', p, x'') = -i hbar pdv(, x') delta(x''-x')
                                                        $
                                                        */
and
$
  braket(beta, p, alpha) &= integral dd(x') braket(beta, x') (- i hbar pdv(, x') braket(x', alpha)) = integral dd(x') psi_beta^* (x') (- i hbar pdv(, x')) psi_alpha (x')
$
/*
by repeated application we can find
$
  braket(x', p^n, alpha) &= (-i hbar)^n pdv(, x', [n]) braket(x', alpha) \
  braket(beta, p^n, alpha) &= integral dd(x') psi_beta^* (x') (-i hbar)^n pdv(, x', [n]) psi_alpha (x')
$
*/

=== Momentum-space
Similarly $p ket(p') = p' ket(p')", "braket(p', p'') = delta(p'-p'')$ and
$
  ket(alpha) = integral dd(p') ket(p') braket(p', alpha)
$
we call $braket(p', alpha) = phi_alpha (p')$ the momentum-space wave function, and the analogous probability interpretation holds.

Recall that given two bases we can construct a transformation matrix. In this case $braket(x', p')$ is called the transformation function from $x arrow p$, this is what we want to find. Consider
$
  braket(x', p, p') & = -i hbar pdv(, x') braket(x', p') = p' braket(x', p') => braket(x', p') = N exp((i p' x')/hbar)
$
this is the momentum eigenstate $ket(p')$ in the $x$-representation---if we treat it as a function of $x'$ with fixed $p'$ then it is just a plain wave. To determine $N$ consider
$
  braket(x', x'') & = integral dd(p') braket(x', p') braket(p', x'') = abs(N)^2 integral dd(p') exp((i p'(x'-x''))/hbar) \
  & = 2 pi hbar abs(N)^2 delta(x'-x'') =^! delta(x'-x'') => N = 1/sqrt((2 pi hbar))
$
where we let $N in RR$ by convention. Now we can relate $psi_alpha (x')$ and $phi_alpha (p')$,
$
  braket(x', alpha) &= integral dd(p') braket(x', p') braket(p', alpha)", " braket(p', alpha) = integral dd(x') braket(p', x') braket(x', alpha)
$
giving
$
  psi_alpha (x') &= 1/sqrt(2pi hbar) integral dd(p') exp((i p'x')/hbar) phi_alpha (p') \
  phi_alpha (p') &= 1/sqrt(2pi hbar) integral dd(x') exp((-i p' x')/hbar) psi_alpha (x')
$
this is just a Fourier transform! The generalization to three dimensions is obvious.
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

#pagebreak()
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
== SHO
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
