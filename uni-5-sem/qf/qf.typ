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

$$

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
which is just the similarity transformation $X' = U^dagger X U$. The trace of an operator is defined as
$
  tr(X) &= sum_a' braket(a', X, a') \
  &= sum_a' sum_b' sum_b'' braket(a', b') braket(b', X, b'')braket(b'', a') \
  &= sum_b' sum_b'' (sum_a' braket(b'', a') braket(a', b')) braket(b', X, b'') \
  &= sum_b' sum_b'' braket(b'', b') braket(b', X, b'') \
  &= sum_b' braket(b', X, b')
$
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
$
  tr(ketbra(b', a')) & = braket(a', b')
$
The problem we just solved is equivalent to finding a unitary matrix which diagonalizes $B$. We are interested in finding $b'$ and $ket(b')$ with
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
So far we have looked at observables with discrete spectra, but many observables---position, momentum, etc.---have continuous spectra. In these cases the ket space becomes infinite, luckily many results carry over, but others don't.

We have the eigenvalue equation $xi ket(xi') = xi' ket(xi')$ as before, most other things get replaced by their continuous analog---sums become integrals, and Kronecker deltas become Dirac delta functions:
$
      braket(xi', xi'') & = delta(xi'-xi'') \
                      1 & =integral dd(xi') ketbra(xi') \
             ket(alpha) & =integral dd(xi') ket(xi') braket(xi', alpha) \
                      1 & = integral dd(xi') abs(braket(xi', alpha))^2 \
    braket(beta, alpha) & = integral dd(xi') braket(beta, xi') braket(xi', alpha) \
  braket(xi'', xi, xi') & = xi' delta(xi''-xi')
$
where we use the completeness relation.

=== $x$ eigenkets
The eigenkets $ket(x')$ satisfy $x ket(x') = x' ket(x')$, and are postulated to form a complete set. So we can expand any state ket $ket(alpha)$ as
$
  ket(alpha) = integral_(-oo)^oo dd(x') ket(x') braket(x', alpha)
$
now suppose we take a measurement using a detector that only clicks when the particle is very close $x'$, then
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
This generalizes to three dimensions, where we assume that the position eigenkets $ket(bold(x)')$ are complete. Then any state can be written as---if we ignore internal degrees of freedom,
$
  ket(alpha) = integral dd(x', 3) ket(bold(x)') braket(bold(x)', alpha)
$
here $ket(bold(x)')$ is a simultaneous eigenket of $x, y$ and $z$
$
  ket(bold(x)') equiv ket(x'","y'","z')" with" x ket(bold(x)') = x' ket(bold(x)')" etc."
$
doing this we assume that $x, y$ and $z$ are compatible
$
  [x_i,x_j] = 0
$

=== Translation and $p$
We want an operator which can take a state from $bold(x)' arrow bold(x)' +dd(bold(x)')$ without affecting anything else. This is the translation operator, and its action is defined by
$
  cal(J)(dd(bold(x)')) ket(bold(x)') = ket(bold(x)'+dd(bold(x)'))
$
clearly $ket(bold(x)')$ is not an eigenket, now consider the action on some arbitrary $ket(alpha)$
$
  cal(J) (dd(bold(x)'))ket(alpha) &= integral dd(x', 3) ket(bold(x)'+dd(bold(x)')) braket(bold(x)', alpha) \
  &= integral dd(x', 3) ket(bold(x)') braket(bold(x)'-dd(bold(x)'), alpha)
$
where we the shift is allowed since we integrate over all of space. So the translated state is obtained by subbing $bold(x)' -> bold(x)' - dd(bold(x)')$.

As all other operators this one has some nice properties. We require that it is unitary, since if $ket(alpha)$ is normalized then $cal(J) (dd(bold(x)')) ket(alpha)$ should also be normalized,
$
  braket(alpha) = braket(alpha, cal(J)^dagger (dd(bold(x)')) cal(J) (dd(bold(x)')), alpha) => cal(J)^dagger (dd(bold(x)')) cal(J) (dd(bold(x)')) = 1
$
in general it is clear that the norm is preserved under unitary transformation---by definition. Some obvious properties:
$
  cal(J) (dd(bold(x)'')) cal(J) (dd(bold(x)')) &= cal(J) (dd(bold(x)') + dd(bold(x)'')) \
  cal(J) (- dd(bold(x)')) &= cal(J)^(-1) (dd(bold(x)')) \
  lim_(dd(bold(x)') arrow 0) cal(J) (dd(bold(x)'))&=1
$
now we want to find $cal(J) (dd(bold(x)'))$.

It turns out if we define
$
  cal(J) (dd(bold(x)')) = 1 - i bold(K) dot dd(bold(x)')
$
where $K_i$ are Hermitian operators, then everything is satisfied. First we check unitarity
$
  cal(J)^dagger (dd(bold(x)')) cal(J) (dd(bold(x)')) &= (1 + i bold(K)^dagger dd(bold(x)'))(1-i bold(K) dot dd(bold(x)')) \
  &= 1 - i(bold(K)-bold(K)^dagger) dot dd(bold(x)') + cal(O) [(dd(bold(x)'))^2] \
  &tilde.eq 1
$
and the additive property
$
  cal(J) (dd(bold(x)'')) cal(J) (dd(bold(x)')) &= (1 - i bold(K) dot dd(bold(x)''))(1 - i bold(K) dot dd(bold(x)')) \
  &tilde.eq 1 - i bold(K) dot (dd(bold(x)') + dd(bold(x)'')) \
  &=cal(J) (dd(bold(x)') + dd(bold(x)''))
$
the other two are immediately satisfied by the definition. Now we want to find a relation between $bold(K)$ and $bold(x)$. Notice
$
  bold(x) cal(J) (dd(bold(x)')) ket(bold(x)') = bold(x) ket(bold(x)' + dd(bold(x)')) = (bold(x)' + dd(bold(x)'))ket(bold(x)' + dd(bold(x)'))
$
and
$
  cal(J) (dd(bold(x)')) bold(x) ket(bold(x)') = bold(x)' cal(J) (dd(bold(x)')) ket(bold(x)') = bold(x)' ket(bold(x)' + dd(bold(x)'))
$
combining these gives
$
  [bold(x),cal(J) (dd(bold(x)'))] ket(bold(x)') = dd(bold(x)') ket(bold(x)' + dd(bold(x)')) tilde.eq dd(bold(x)') ket(bold(x)')
$
so we have found the identity $ [bold(x),cal(J) (dd(bold(x)'))] = dd(bold(x)') $
which implies
$
  -i bold(x) bold(K) dot dd(bold(x)') + i bold(K) dot dd(bold(x)') bold(x) = dd(bold(x)')
$
here $dd(bold(x)')$ is a number multiplied by the identity operator in the ket space spanned by $ket(bold(x)')$. Letting $dd(bold(x)')$ be along $hat(bold(x))_j$ and forming the scalar product with $hat(bold(x))_i$ we obtain
$
  [x_i, K_j] = i delta_(i j)
$
but what is $bold(K)$ physically?

In classical mechanics a similar translation can be seen as a canonical transformation
$
  bold(x)_"new" equiv bold(X) = bold(x) + dd(bold(x))",  " bold(p)_"new" equiv bold(P) = bold(p)
$
with the generating function
$
  F(bold(x),bold(P)) = bold(x) dot bold(P) + bold(p) dot dd(bold(x))
$
this looks like the translation operator we defined, especially since $bold(x) dot bold(P)$ is the generating function for the identity transformation---so maybe $bold(K)$ is related to momentum? For this to be the case we need to fix units since $bold(K)$ has units of $1\/"length"$,
$
  bold(K) = bold(p)/"constant with unit of action"
$
so what constant do we choose? Consider de Broglie's relation
$
  (2pi)/lambda = p/hbar
$
the universal constant we are after turns out to be $hbar$, and the $K$ operator is the operator corresponding to the wave number, typically denoted by $k$. With this identification we have found
$
  cal(J) (dd(bold(x)')) = 1 - i bold(p) dot dd(bold(x)')\/hbar
$
where $bold(p)$ is the momentum operator, and the commutation relation becomes
$
  [x_i, p_j] = i hbar delta_(i j)
$
which as we know is very powerful.

Now we see what happens for a finite translation $dd(x', d: Delta)$? We can compound $N -> oo$ translations with displacement $dd(x', d: Delta)\/N$ and we obtain
$
  cal(J) (dd(x', d: Delta) hat(bold(x))) &= lim_(N arrow oo) (1 - (i p_x dd(x', d: Delta))/(N hbar))^N \
  &= exp(- (i p_x dd(x', d: Delta))/hbar)
$
where we define
$
  exp(X) equiv 1 + X + X^2/2! + dots
$
a fundamental property of translations is that successive translations in different directions commute
$
  [cal(J) (dd(y', d: Delta) hat(bold(y))), cal(J) (dd(x', d: Delta) hat(bold(x))) ] = 0
$
using this we obtain
$
  [cal(J) (dd(y', d: Delta) hat(bold(y))), cal(J) (dd(x', d: Delta) hat(bold(x))) ] &= [(1-(i p_y dd(y', d: Delta))/hbar + dots), (1-(i p_x dd(x', d: Delta))/hbar + dots)] \
  &tilde.eq - ((dd(x', d: Delta)) (dd(y', d: Delta)) [p_y,p_x])/hbar^2
$
or since the displacements are arbitrary we find
$
  [p_x,p_y]=0 => [p_i,p_j] = 0
$
this commutation relation is a direct result of the fact that translations in different directions commute. When the generators of transformations commute the corresponding group is Abelian---so the translation group in three dimensions is Abelian.

This implies we can find a simultaneous eigenket of $p_x, p_y$ and $p_z$
$
  ket(bold(p)') equiv ket(p'_x","p'_y","p'_z)",  " p_x ket(bold(p)') = p'_x ket(bold(p)') "etc."
$
we can find the effect of the translation operator on this state
$
  cal(J) (dd(bold(x)')) ket(bold(p)') &= (1 - (i bold(p) dot dd(bold(x)'))/hbar) ket(bold(p)') \
  &= (1 - (i bold(p)' dot dd(bold(x)'))/hbar) ket(bold(p)')
$
so $ket(bold(p)')$ is an eigenket of $cal(J) (dd(bold(x)'))$, which is obvious since
$
  [bold(p),cal(J) (dd(bold(x)'))] = 0
$
but the eigenvalues are complex, since $cal(J) (dd(bold(x)'))$ is unitary.

What we have just found are the canonical commutation relations:
$
  [x_i,x_j] = 0", " [p_i,p_j] = 0", " [x_i,p_j] = i hbar delta_(i j)
$
following Dirac one might notice that these are very similar to Poisson brackets
$
  [dot,dot]_"class" -> [dot,dot]/(i hbar)
$
this is plausible since the commutator and Poisson bracket share very similar algebraic propeties---these could have been used to immediately find the canonical commutation relations, but using the translation operator is much nicer.

#pagebreak()
== Wave-functions
=== Position space
In one dimension we have
$
  x ket(x') = x' ket(x')
$
with
$
  braket(x'', x') = delta(x''-x')
$
any state can then be written
$
  ket(alpha) = integral dd(x') ket(x') braket(x', alpha)
$
with $abs(braket(x', alpha))^2 dd(x')$ being the probability to find our particle near $x'$. In our formalism we refer to $braket(x', alpha)$ as the wave function $psi_alpha (x')$ for the state $ket(alpha)$---so the wavefunction is just an expansion coefficient.

Consider the inner product
$
  braket(beta, alpha) = integral dd(x') braket(beta, x') braket(x', alpha) = integral dd(x') psi_beta^* (x') psi_alpha (x')
$
so $braket(beta, alpha)$ characterizes overlap between our wavefunctions---note that in general $braket(beta, alpha)$ is independent of representations, and represents the probability for $ket(alpha)$ to be found in $ket(beta)$.

We can also look at
$
  ket(alpha) = sum_a' ket(a') braket(a', alpha)
$
and obtain
$
  braket(x', alpha) = sum_a' braket(x', a') braket(a', alpha)
$
which is typically written $psi_alpha (x') = sum_a' c_a' u_a' (x')$ where $u_a' (x') = braket(x', a')$ is an eigenfunction of an operator $A$, in the $x$-representation, with eigenvalue $a'$. Similarly
$
  braket(beta, A, alpha) &= integral dd(x') integral dd(x'') braket(beta, x') braket(x', A, x'') braket(x'', alpha) \
  &= integral dd(x') integral dd(x'') psi_beta^* (x') braket(x', A, x'') psi_alpha (x'')
$
if $A = f(x)$ we get a delta-function eating one of the integrals and
$
  braket(beta, f(x), alpha) = integral dd(x') psi_beta^* (x') f(x') psi_alpha (x')
$

It would be nice to know what $p$ looks like in this basis. We start with momentum as the generator of translations
$
  (1 - (i p dd(x', d: Delta))/hbar) ket(alpha) &= integral dd(x') cal(J) (dd(x', d: Delta)) ket(x')) braket(x', alpha) \
  &= integral dd(x') ket(x') braket(x' - dd(x', d: Delta), alpha) \
  &= integral dd(x') ket(x') (braket(x', alpha) - dd(x', d: Delta) pdv(, x') braket(x', alpha) )
$
comparison gives
$
  p ket(alpha) = integral dd(x') ket(x') (- i hbar pdv(, x') braket(x', alpha)) => braket(x', p, alpha) = - i hbar pdv(, x') braket(x', alpha)
$
so this is the action of $p$ in the $x$-representation, the matrix element is
$
  braket(x', p, x'') = -i hbar pdv(, x') delta(x''-x')
$
we can also find
$
  braket(beta, p, alpha) &= integral dd(x') braket(beta, x') (- i hbar pdv(, x') braket(x', alpha)) \
  &= integral dd(x') psi_beta^* (x') (- i hbar pdv(, x')) psi_alpha (x')
$
by repeated application we can find
$
  braket(x', p^n, alpha) &= (-i hbar)^n pdv(, x', [n]) braket(x', alpha) \
  braket(beta, p^n, alpha) &= integral dd(x') psi_beta^* (x') (-i hbar)^n pdv(, x', [n]) psi_alpha (x')
$

=== Momentum-space
Similarly we have
$
  p ket(p') = p' ket(p')",  " braket(p', p'') = delta(p'-p'')
$
these like ${ket(x')}$ span the space so for arbitrary $ket(alpha)$ we can write
$
  ket(alpha) = integral dd(p') ket(p') braket(p', alpha)
$
we call $braket(p', alpha) = phi_alpha (p')$ the momentum-space wave function, and the analogous probability interpretation holds.

Recall that given two bases we can construct a transformation matrix. In this case $braket(x', p')$ is called the transformation function from $x arrow p$, this is what we want to find. Consider
$
  braket(x', p, p') & = -i hbar pdv(, x') braket(x', p') \
                    & = p' braket(x', p')
$
the solution is
$
  braket(x', p') = N exp((i p' x')/hbar)
$
before we normalize, it is worth mentioning that this is just the momentum eigenstate $ket(p')$ in the $x$-representation---if we treat it as a function of $x'$ with fixed $p'$---so the wave function of a momentum eigenstate is just a plain wave. To find $N$ consider
$
  braket(x', x'') & = integral dd(p') braket(x', p') braket(p', x'') \
  & = abs(N)^2 integral dd(p') exp((i p'(x'-x''))/hbar) \
  & = 2 pi hbar abs(N)^2 delta(x'-x'') =^! delta(x'-x'') => N = 1/sqrt((2 pi hbar))
$
where we let $N in RR$ by convention. So we have just found
$
  braket(x', p') = 1/(sqrt(2 pi hbar)) exp((i p' x')/hbar)
$
now we can relate $psi_alpha (x')$ and $phi_alpha (p')$,
$
  braket(x', alpha) &= integral dd(p') braket(x', p') braket(p', alpha)", " braket(p', alpha) = integral dd(x') braket(p', x') braket(x', alpha)
$
give
$
  psi_alpha (x') &= 1/sqrt(2pi hbar) integral dd(p') exp((i p'x')/hbar) phi_alpha (p') \
  phi_alpha (p') &= 1/sqrt(2pi hbar) integral dd(x') exp((-i p' x')/hbar) psi_alpha (x')
$
this is just a Fourier transform!

The generalization to three dimensions is obvious
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

#pagebreak()
= Quantum Dynamics
