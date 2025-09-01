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

Light analogy leads us to consider a complex vector space.\*

#pagebreak()
== Bras and Kets
We consider a complex vector space with dimensionality depending on the given system we are looking at---in the Stern-Gerlach experiment the dimensionality would be two. We consider the case of infinite dimensionality later, which is appropriate for something like position or momentum. The vector space in question is known as a Hilbert space $cal(H)$, which is a complete inner product space---the reason for this will become clear since it is complete and has an inner product. Any physical state in quantum mechanics is represented by a state vector in this complex vector space. These we call kets and are denoted by $ket(alpha)$ following Dirac who introduced them. Everything we'd like to know about a given state is contained in the corresponding ket. Since these are members of a vector space we can add them $ket(alpha)+ket(beta)=ket(gamma)$ and multiply them by some complex number $c$ to give new kets---note that $c ket(alpha)=ket(alpha)c$, if $c=0$ we get the null ket. As physicists we postulate that $ket(alpha)$ and $c ket(alpha)$ both correspond to the same physical state.

All observables can be represented by operators that act on the vector space. In general for some operator $A$ we have $A dot (ket(alpha)) = A ket(alpha)$, with $A ket(alpha)$ also being a ket. We define eigenkets of an operator $ ket(a'), ket(a''), dots $ as those with the property
$
  A ket(a') = a' ket(a')",  " A ket(a'') = a'' ket(a''), dots
$
where $a', a'', dots$ are numbers, namely the eigenvalues of the given operator---typically denoted as ${a'}$. The physical state corresponding to the eigenket we call the eigenstate. As an example take the case of spin-$1\/2$ systems where we have
$
  S_z ket(S_z";"plus.minus) = plus.minus hbar/2 ket(S_z";"plus.minus)
$
and similarly for $ket(S_x";"plus.minus)$ and $ket(S_y";"plus.minus)$, which are eigenkets of the operators $S_x$ and $S_y$ respectively. Usually we are interested in an $N$-dimensional vector space spanned by $N$ eigenkets of some observable $A$, in this case any ket $ket(alpha)$ can be written as $ket(alpha) = sum_(a') c_(a') ket(a')$.

So far we've been looking at the ket space, now we introduce the bra space---which is the dual to the ket space. We postulate that every ket $ket(alpha)$ has a corresponding bra $bra(alpha)$ in the bra space---which is spanned by eigenbras ${bra(a')}$ corresponding to the eigenkets ${ket(a')}$. Essentially we have a dual correspondence between the two spaces. We postulate that the bra dual to $c ket(alpha)$ is $c^* bra(alpha)$, so we have
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
with $braket(tilde(alpha), tilde(alpha))=1$. The root is known as the norm of $alpha$. Since $ket(alpha)$ and $c ket(alpha)$ both represent the same physical state we typically normalize everything.

Beside operators representing observables we also have more general operators, denote these by $X, Y, Z, dots$. These act on kets from the left $X ket(alpha)$ and on bras from the right $bra(beta) X$. Two operators are equal $X=Y$ if $X ket(alpha)= Y ket(beta)$. An operator is the null operator if $X ket(alpha) = 0$. Operators can be added and addition is commutative and associative, further most of the operators we consider are linear,
$
  X (c_alpha ket(alpha) + c_beta ket(beta)) = c_alpha X ket(alpha) + c_beta X ket(beta)
$
we define the Hermitian adjoint as the dual operator $X ket(alpha) <->^"DC" bra(alpha)X^dagger$, if $X = X^dagger$ then $X$ is Hermitian.

Two operators can be multiplied, but this is non-commutative in general, but it is associative and
$
  X(Y ket(alpha)) = (X Y)ket(alpha)=X Y ket(alpha)
$
similarly for bras. Notice also $(X Y)^dagger = Y^dagger X^dagger$, since
$
  X Y ket(alpha) = X(Y ket(alpha)) <->^"DC" (bra(alpha)Y^dagger)X^dagger = bra(alpha) Y^dagger X^dagger
$
Another product we can consider is the outer product $ketbra(beta, alpha)$, this is an operator---it is also worth noting that there are many nonsensical product, e.g. $ket(alpha)ket(beta)$ when both kets are from the same space.

Multiplication amongst operators is associative, but we postulate that it holds among everything, so kets, bras and operators are all associative under multiplication---Dirac calls this the associative axiom of multiplication. So we can say that $(ketbra(beta, alpha)) dot ket(gamma) = ket(beta) (braket(alpha, gamma))$, where the inner product is just a number, so acting with the outer product gives us another ket, so the outer product is indeed an operator. Notice that if $X = ketbra(beta, alpha)$ then $X^dagger = ketbra(alpha, beta)$. As a point of notation we have $(bra(beta))dot (X ket(alpha)) = (bra(beta)X) dot (ket(alpha)) = braket(beta, X, alpha)$ and
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

On physical grounds we expect observables to have real eigenvalues---this theorem guarantees that eigenvalues are real if the operator is Hermitian. This is why we care about Hermitian observables. Usually we normalize $ket(a')$ so that ${ket(a')}$ form an orthonormal set $braket(a'', a')=delta_(a'' a')$. By construction this set is complete, since we asserted at the very beginning that the ket space is spanned by the eigenkets of $A$---then by definition we can use the eigenkets as basekets. Given some arbritrary $ket(alpha)$ we can write
$
  ket(alpha) = sum_(a') c_(a') ket(a')
$
multiplying on the left by $bra(a'')$ we find $c_(a') = braket(a', alpha)$ so
$
  ket(alpha) = sum_(a') ket(a') braket(a', alpha)
$
by associativity we can $ket(a') braket(a', alpha)$ as the operator $ketbra(a', a')$ acting on $ket(alpha)$. Since $ket(alpha)$ is arbitrary we get the completeness relation
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
Consider again the operator $ketbra(a', a')$ and let it act on $ket(alpha)$
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
we can write $braket(a'', X, a') = braket(a', X^dagger, a'')^*$, so the Hermitian conjugate has been related to the complex conjugate transpose. This is also in agreement with the usual way we do matrix multiplication $Z = X Y$ becomes
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
For an observable this representation becomes very simple if the eigenkets are used as the base kets. We have
$
  A & = sum_(a'') sum_(a') ket(a'') braket(a'', A, a') bra(a') \
    & = sum_(a') a' ketbra(a', a') = sum_(a') a' Lambda_(a')
$
since
$
  braket(a'', A, a') = braket(a', A, a') delta_(a'' a') = a' delta_(a'' a')
$

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
these are ladder operators and their effet is obvious. From these we can construct $S_x$ and $S_y$ by $S_(plus.minus) = S_x plus.minus i S_y$, which we show later. Now we can construct the matrix representation
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
sometimes we write $expval(A)_alpha$, this definition agree with what we'd expect
$
  expval(A) = sum_(a') sum_(a'') braket(alpha, a'') braket(a'', A, a') braket(a', alpha) = sum_(a') a' abs(braket(a', alpha))^2
$
so it's the measured value times the probability of measuring it.

=== Spin-$1\/2$
We visit the spin-$1\/2$ system again. In the Stern-Gerlach sending a $S_x +$ beam into $"SG"hat(z)$ splits it into two equal parts. So
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
it is simplest to pick $delta_1=0$ and then $delta_2 = pi\/2$ turns out to be correct. So to summarize
$
  ket(S_x";"plus.minus) &= 1/sqrt(2) ket(+) plus.minus 1/sqrt(2) ket(-)",  " ket(S_y";"plus.minus) = 1/sqrt(2) ket(+) plus.minus i/sqrt(2) ket(-) \
$
$
  S_x = hbar/2 [(ketbra(+, -))+(ketbra(-, +))]",  " S_y = hbar/2 [-i (ketbra(+, -))+i(ketbra(-, +))]
$
and we see now
$
  S_plus.minus = S_x plus.minus i S_y
$
further these satisfy the commutation relations
$
  [S_i,S_j] = i epsilon_(i j k) hbar S_k
$
and the anticommutation relations
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
this is just an eigenvalue equation for $B$ with eigenvalue $braket(a', B, a')$. So $ket(a')$ is a simultaneous eigenket of $A$ and $B$---we sometimes write $ket(a'b')$, or $ket(K')$ where $K'$ is a collective index.

This was for the case of some $A$ with non-degenerate eigenkets, this is also valid for $n$-fold degeneracy
$
  A ket(a'^((i))) = a' ket(a'^((i)))"  for " i=1,2,dots,n
$
where $ket(a'^((i)))$ are $n$ mutually orthonormal eigenkets of $A$ with the same eigenvalue $a'$---to see this we just need to construct appropriate linear combinations of the eigenkets that diagonalize the $B$ operator.

We can generalize this easily by considering more observables which are mutually compatible
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
where $n$ is the degeneracy. Measureing $B$ would then pick one of the terms, but measuring $A$ again would still give $a'$. So they don't interfere, which is why we call them compatible.

Incompatible observables are more annoying. These don't have a complete set of simultaneous eigenkets. To show this we can proceed by contradiction, if it were true there would exist a set of simultaneous eigenkets satisfying
$
  A ket(a' b') = a' ket(a' b')",  " B ket(a' b') = b' ket(a' b')
$
so
$
  A B ket(a' b') = A b' ket(a' b') = a' b' ket(a' b')",  " B A ket(a' b') = B a' ket(a' b') = a' b' ket(a' b')
$
or $A B ket(a' b') = B A ket(a' b')$ implying that $[A, B]=0$ which is a contradiction. However, this can be true in a ket subspace, even if $A$ and $B$ are incompatible.

Consider a measurement of $C$ taking $ket(a') arrow^"C" ket(c')$, and in the case where we perform a selective measurement first $ket(a') arrow^"B" ket(b') arrow^"C" ket(c')$. Since probabilites are multiplicative we have in the latter case probability $abs(braket(c', b'))^2 abs(braket(b', a'))^2$ to measure $ket(c')$. To get the total probability we sum over all $b'$ since this would seem to remove the filter
$
  sum_(b') abs(braket(c', b'))^2 abs(braket(b', a'))^2 = sum_(b') braket(c', b') braket(b', a') braket(a', b') braket(b', c')
$
the total probability in the first case is just
$
  abs(braket(c', a'))^2 = abs(sum_(b') braket(c', b')braket(b', a'))^2 = sum_(b') sum_(b'') braket(c', b')braket(b', a') braket(a', b'')braket(b'', c')
$
these are very different. Which is cool since in each case the $ket(a')$ can be regared as built from $ket(b')$, $ket(a') = sum_(b') ket(b') braket(b', a')$. This shows how taking a measurement of $B$ actually changes things---in the first we actually record which eigenvalues were realized, in the second we just imagine that $ket(a')$ is built from the $ket(b')$. So actually writing down the different $b'$ routes changes the result. These expressions are however equal for compatible observables as we would expect, i.e. in the case where $[A,B] = [B,C] = 0$. To see this consider
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

=== Uncertainty

#pagebreak()
= Quantum Dynamics

#pagebreak()
= Angular Momentum
