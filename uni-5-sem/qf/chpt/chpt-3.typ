//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Theory of Angular Momentum
== Commutators
=== Classical rotation
As should be obvious some rotations commute while others don't---rotations about the same axis always commute, we want to know how rotations about different axes fail. Rotations can be represented by 3 by 3 orthogonal matrices $R$---ensuring length is conserved. An infinitesimal rotation is of the form
$
  R_z (epsilon) = mat(1-epsilon^2/2, -epsilon, 0; epsilon, 1-epsilon^2/2, 0; 0, 0, 1)
$
and similarly for $x$ and $y$. If we stick with order $epsilon^2$ then we can find
$
  R_x (epsilon) R_y (epsilon) - R_y (epsilon) R_x (epsilon) = R_z (epsilon^2) - R_"any" (0)
$
with $R_"any" (0) = 1$, we use this later.
=== Rotation operator
In QM we associate an operator $cal(D) (R)$ to each rotation matrix such that
$
  ket(alpha)_R = cal(D) (R) ket(alpha)
$
with $ket(alpha)_R$ being the state in the rotated system---we'd like to construct this operator. For translations and time-evolution we had operators of the form
$
  U_epsilon = 1 - i G epsilon
$
with $G$ being either $p$ or $H$ giving $epsilon -> {dd(x)', dd(t)}$. Classically angular momentum is the generator of rotations, so we define the angular momentum operator $J_k$ such that
$
  G -> J_k/hbar => epsilon -> dd(phi.alt)
$
with $J_k$ Hermitian the operator is guaranteed to be unitary and reduces to the identity in the limit $dd(phi.alt) -> 0$. So we obtain
$
  cal(D) (hat(bold(n)),dd(phi.alt)) = 1 - i ((bold(J)dot hat(bold(n)))/hbar) dd(phi.alt)
$
for a rotation about $hat(bold(n))$ by an angle $dd(phi.alt)$. As with the other infinitesimal operators we can get a finite rotation by compounding infinitesimal rotations---take a rotation about the $z$-axis
$
  cal(D)_z (phi.alt) &= lim_(N -> oo) [1 - i (J_z/hbar) (phi.alt/N)]^N \
  &= exp((-i J_z phi.alt)/hbar) \
  &= 1 - (i J_z phi.alt)/hbar - (J_z^2 phi.alt^2)/(2 hbar^2) + dots
$
To get the commutation relations we postulate that $cal(D) (R)$ has the same group properties as $R$---e.g.
$
  R_1 R_2 = R_3 & => cal(D) (R_1) cal(D) (R_2) = cal(D) (R_3) \
   R R^(-1) = 1 & => cal(D) (R) cal(D)^(-1) (R) = 1
$
and associativity.

=== Commutation relations
We take the classical commutator but with $R -> cal(D)(R)$---this is a mess, but equating terms of order $epsilon^2$ we obtain
$
  [J_i, J_j] = i hbar epsilon_(i j k) J_k
$
which is the fundamental commutation relations of angular momentum---so the group of rotations is non-Abelian.

To see the operator actually rotates our state consider
$
  expval(J_x) -> braket(alpha, cal(D)^dagger_z (phi.alt) J_x cal(D)_z (phi.alt), alpha)
$
and we have
$
  exp((i J_z phi.alt)/hbar) J_x exp((-i J_z phi.alt)/hbar) &= dots \
  &= J_x [1 - phi.alt^2/2! + dots] - J_y [phi.alt - phi.alt^3/3! + dots]\
  &= J_x cos phi.alt - J_y sin phi.alt
$

== Spin $1\/2$ systems
=== Rotation operator
We have seen that
$
  S_x & = hbar/2 {ketbra(+, -) + ketbra(-, +)} \
  S_y & = (i hbar)/2 {-ketbra(+, -) + ketbra(-, +)} \
  S_z & = hbar/2 {ketbra(+, +) - ketbra(-, -)}
$
satisfy the commutation relations with $J_k -> S_k$---so the relations are realized for $N=2$ dimensions. We also just showed that
$
  expval(S_x) -> expval(S_x) cos phi.alt - expval(S_y) sin phi.alt
$
similarly
$
  expval(S_y) -> expval(S_y) cos phi.alt + expval(S_x) sin phi.alt
$
and since $[S_z, cal(D)_z (phi.alt)] = 0$
$
  expval(S_z) -> expval(S_z)
$
this shows that the expectation values behave like a classical vector under rotation---in general
$
  expval(J_k) -> sum_l R_(k l) expval(J_l)
$
Now consider a general ket
$
  ket(alpha) = ket(+) braket(+, alpha) + ket(-) braket(-, alpha)
$
we see that
$
  exp((-i S_z phi.alt)/hbar) ket(alpha) = e^(-i phi.alt \/2) ket(+) braket(+, alpha) + e^(i phi.alt\/2) ket(-) braket(-, alpha)
$
for a rotation of $2 pi$ we see that $ket(alpha) -> - ket(alpha)$, so to get the same ket back we need a rotation of $4 pi$ or $720 degree$---note that this minus sign cancels when we compute expectation values.

=== Spin precession
We saw that the time-evolution operator was
$
  cal(U) (t,0) = exp((-i S_z omega t)/hbar)
$
but this is just the rotation operator with $phi.alt = omega t$! So we immediately obtain
$
  expval(S_x)_t = expval(S_x)_(t = 0) cos omega t - expval(S_y)_(t=0) sin omega t
$
etc. so we clearly see precession---the ket itself will be given by
$
  ket(alpha","t_0 = 0";"t) = e^(-i omega t\/2) ket(+) braket(+, alpha) + e^(i omega t\/2) ket(-) braket(-, alpha)
$
it's worth noting that
$
  tau_"precession" = (2 pi)/omega" and" tau_"state" = (4 pi)/omega
$
the _extra_ minus sign can be shown by experiment to be a real thing.

=== Pauli Formalism
Following Pauli we denote
$
  ket(+) eq^dot vec(1, 0) equiv chi_+" and " bra(+) eq^dot vecrow(1, 0) equiv chi_+^dagger
$
with obvious representation for the spin-down state. Then
$
  ket(alpha) eq^dot vec(braket(+, alpha), braket(-, alpha)) = chi equiv (c_+,c_-)
$
with $chi$ being the two-component spinor. The matrix elements of $S_k$ are written in terms of the Pauli matrices
$
  braket(plus.minus, S_k, +) equiv hbar/2 (sigma_k)_(plus.minus, plus)" and " braket(plus.minus, S_k, -) equiv hbar/2 (sigma_k)_(plus.minus,-)
$
so we can write
$
  expval(S_k) = braket(alpha, S_k, alpha) &= sum_(a',a'' = plus,minus) braket(alpha, a') braket(a', S_k, a'') braket(a'', alpha) \
  &= hbar/2 chi^dagger sigma_k chi
$
explicitly we can find
$
  sigma_1 = mat(0, 1; 1, 0)",  " sigma_2 = mat(0, -i; i, 0)",  " sigma_3 = mat(1, 0; 0, -1)
$
notably we have the relations
$
  {sigma_i, sigma_j} & = 2 delta_(i j) \
  [sigma_i, sigma_j] & = 2 i epsilon_(i j k) sigma_k
$
these have many other nice properties, e.g. they are Hermitian, have determinant $-1$ and are traceless.

Consider $bold(sigma) dot bold(a)$
$
  bold(sigma) dot bold(a) equiv sum_k a_k sigma_k = mat(a_3, a_1 -i a_2; a_1+i a_2, - a_3)
$
we have the identity
$
  (bold(sigma) dot bold(a))(bold(sigma) dot bold(b)) = bold(a) dot bold(b) + i bold(sigma) dot (bold(a) times bold(b))
$
#proof[
  $
    sum_j sigma_j a_j sum_k sigma_k b_k &= sum_(j,k) (1/2 {sigma_j, sigma_k} + 1/2 [sigma_j,sigma_k]) a_j b_k \
    &= sum_(j, k) (delta_(j k) + i epsilon_(j k l) sigma_l) a_j b_k \
    &= bold(a) dot bold(b) + i bold(sigma) dot (bold(a) times bold(b))
  $
]
given $a_k in RR$ we have $(bold(sigma) dot bold(a))^2 = abs(bold(a))^2$.

==== Rotation operator
We have
$
  cal(D) (hat(bold(n)),phi.alt) = exp((-i bold(S) dot hat(bold(n)) phi.alt)/hbar) eq^dot exp((-i bold(sigma) dot hat(bold(n))phi.alt)/2)
$
with
$
  (bold(sigma) dot hat(bold(n)))^n = cases(
    1 #h(40pt) & "for" n "even",
    bold(sigma) dot hat(bold(n)) & "for" n "odd"
  )
$
so we can write
$
  exp((-i bold(sigma) dot hat(bold(n)) phi.alt)/2) &= dots \
  &= II cos phi.alt/2 - i bold(sigma) dot hat(bold(n)) sin phi.alt/2 \
  &= mat(cos phi.alt/2 - i n_z sin phi.alt/2, (- i n_x - n_y) sin phi.alt/2; (-i n_x + n_y) sin phi.alt/2, cos phi.alt/2 + i n_z sin phi.alt/2)
$
this acts on the two-component spinor $chi$
$
  chi -> exp((-i bold(sigma)dot hat(bold(n))phi.alt)/2) chi
$
note that $sigma_k$ remain unchanged under rotations---but
$
  chi^dagger sigma_k chi -> sum_l R_(k l) chi^dagger sigma_l chi
$
we also have
$
  evaluated(exp((-i bold(sigma) dot hat(bold(n))phi.alt)/2))_(phi.alt = 2 pi) = -1 "for any" hat(bold(n))
$
Now consider
$
  bold(sigma) dot hat(bold(n)) chi = chi
$
so we want the eigenspinor with eigenvalue $1$. So we want
$
  bold(S) dot hat(bold(n)) ket(bold(S) dot hat(bold(n))";" +) = hbar/2 ket(bold(S) dot hat(bold(n))";"+)
$
we let the polar and azimuthal angles characterizing $hat(bold(n))$ be $beta$ and $alpha$. We start with $vec(1, 0)$ representing spin-up and then we first rotate about $y$ by $beta$, and then about $z$ by $alpha$. This is equivalent to
$
  chi &= (cos alpha/2 - i sigma_3 sin alpha/2) (cos beta/2 - i sigma_2 sin beta/2) vec(1, 0) \
  &= vec(cos beta/2 e^(-i alpha\/2), sin beta/2 e^(i alpha\/2))
$
which is what we previously found.

== Rotations and Groups
Rotations can be represented by $3 times 3$ orthorgonal matrices $R$---$R R^T = 1$. All orthogonal matrices form a group under multiplication---all group axioms are trivially satisfied. This group is typically denoted $"SO"(3)$ for special orthogonal group in three dimensions---we require $det R = 1$ since we don't care about inversions.

We can also represent rotations by unitary $2 times 2$ matrices acting on two-component spinors---this is what we saw in the previous section with
$
  cal(D) =^dot exp((-i bold(sigma) dot hat(bold(n)) phi.alt)/2)
$
acting on $chi$, unitarity ensures that $abs(c_+)^2 + abs(c_-)^2 = 1$ is left invariant. It's also easy to verify that this matrix has determinant $1$. Consider a general such matrix
$
  U (a,b) = mat(a, b; -b^*, a^*)
$
with $abs(a)^2+abs(b)^2 = 1$---trivially then $U^dagger (a,b) U(a,b) = 1$. By comparison we can find
$
  Re a & = cos phi.alt/2",  " Im a = -n_z sin phi.alt/2 \
  Re b & = -n_y sin phi.alt/2",  " Im b = -n_x sin phi.alt/2
$
the $a$ and $b$ are known as Cayley-Klein parameters. It's also trivial to show that these form a group namely $"SU"(2)$ for special unitary group in two dimensions. A general member of $"U"(2)$ can be written as
$
  U = e^(i gamma) mat(a, b; -b^*, a^*)
$
with $abs(a)^2+abs(b)^2=1$ and $gamma = gamma^*$---this is why the previously discussed gauge symmetry corresponds to $"U"(1)=e^(i gamma)$. All these groups are of course subgroups of $"GL"_n (CC)$. Note that even though both $"SO"(3)$ and $"SU"(2)$ can represent rotations then they are not isomorphic---by counterexample consider a $2 pi$ and $4 pi$ rotation both of these would be the identity in $"SO"(3)$ but in $"SU"(2)$ they correspond to the identity and minus the identity---in general $U(a, b)$ and $U(-a,-b)$ correspond to the same member in $"SO"(3)$.

Classically any rotation can be accomplished using three rotations and Euler angles by
$
  R(alpha,beta,gamma) = R_z (alpha) R_y (beta) R_z (gamma)
$
applying this to spin-$1\/2$ we find
$
  cal(D) (alpha,beta,gamma) = cal(D)_z (alpha) cal(D)_y (beta) cal(D)_z (gamma)
$
which has the representation
$
  cal(D) (alpha,beta,gamma) =^dot mat(e^(-i(alpha+gamma)\/2) cos beta/2, -e^(-i(alpha-gamma)\/2) sin beta/2; e^(i(alpha-gamma)\/2) sin beta/2, e^(i(alpha+gamma)\/2) cos beta/2)
$
this has the form of $U(a,b)$---and is called the $j = 1\/2$ irreducible representation of $cal(D) (alpha,beta,gamma)$ with matrix elements
$
  cal(D)_(m'm)^(1\/2) (alpha,beta,gamma) = braket(j=1/2","m', exp((-i J_z alpha)/hbar) exp((-i J_y beta)/hbar) exp((-i J_z gamma)/hbar), j =1/2","m)
$

#pagebreak()
== Eigenvalues and -states
=== Ladder operators
We define
$
  bold(J)^2 equiv J_x J_x + J_y J_y + J_z J_z
$
this commutes with all $J_k$
$
  [bold(J)^2, J_k] = 0
$

#proof[
  Consider just $k = z$
  $
    [J_x J_x + J_y J_y + J_z J_z, J_z] &= J_x [J_x, J_z] + [J_x, J_z] J_x + J_y [J_y, J_z] + [J_y,J_z] J_y \
    &= J_x (- i hbar J_y) + (-i hbar J_y) J_x + J_y (i hbar J_x) + (i hbar J_x) J_y \
    &= 0
  $
  by permutation it holds for $k = x,y$ as well.
]

This is nice since we can find simultaneous eigenstates for $bold(J)^2$ and $J_z$ (by convention). We denote
$
  bold(J)^2 ket(a","b) = a ket(a","b)",  " J_z ket(a","b) = b ket(a","b)
$
to determine $a, b$ we define the ladder operators
$
  J_plus.minus equiv J_x plus.minus i J_y
$
these satisfy
$
  [J_+ ,J_-] = 2 hbar J_z",  " [J_z,J_plus.minus] = plus.minus hbar J_z
$
which can be easily proven---note also $[bold(J)^2, J_plus.minus] = 0$. Now consider
$
  J_z (J_plus.minus ket(a","b)) &= ([J_z, J_plus.minus]+J_plus.minus J_z) ket(a","b) \
  &= (b plus.minus hbar) (J_plus.minus ket(a","b))
$
so applying $J_plus.minus$ to a $J_z$ eigenket gives a new eigenket of $J_z$ with eigenvalue $b plus.minus hbar$---note that however
$
  bold(J)^2 (J_plus.minus ket(a","b)) & = J_plus.minus bold(J)^2 ket(a","b) \
                                      & = a (J_plus.minus ket(a","b))
$
since they commute. So $J_plus.minus ket(a","b)$ are simultaneous eigenkets of $bold(J)^2$ and $J_z$ with eigenvalues $a$ and $b plus.minus hbar$---this also lets us write
$
  J_plus.minus ket(a","b) = c_plus.minus ket(a"," b plus.minus hbar)
$

=== Eigenvalues
Consider applying $J_+$ $n$ times to a simultaneous eigenket. Then we obtain an eigenket with the $J_z$ eigenvalue increased by $n hbar$ amd the $bold(J)^2$ eigenvalue unchanged. This process has an upper limit by
$
  a >= b^2
$

#proof[
  To see this note
  $
    bold(J)^2 - J_z^2 & = 1/2 (J_+ J_- + J_- J_+) \
                      & = 1/2 (J_+ J_+^dagger + J_+^dagger J_+)
  $
  the $J_+ J_+^dagger$ and $J_+^dagger J_+$ have non-negative expectation values so
  $
    braket(a","b, bold(J)^2-J_z^2, a","b) >= 0
  $
  or
  $
    braket(a","b, a - b^2, a","b) >= 0 => a >= b^2
  $
]

It follows that there is some $b_"max"$ with
$
  J_+ ket(a","b_"max") = 0 => J_- J_+ ket(a","b_"max") = 0
$
but
$
  J_- J_+ & = J_x^2 + J_y^2 - i (J_y J_x - J_x J_y) \
          & = bold(J)^2 - J_z^2 - hbar J_z
$

so
$
  (bold(J)^2 - J_z^2 - hbar J_z) ket(a","b_"max") = 0
$
this is only possible given
$
  a - b_"max"^2 - b_"max" hbar = 0 => a = b_"max" (b_"max" + hbar)
$
analogously there is some
$
  J_- ket(a","b_"min") = 0 => J_+ J_- ket(a","b_"min") = 0
$
with
$
  J_+ J_- = bold(J)^2 - J_z^2 + hbar J_z => a = b_"min" (b_"min" - hbar)
$
by comparison $b_"max" = - b_"min" > 0$. So the allowed values of $b$ lie in
$
  -b_"max" <= b <= b_"max"
$
we should be able to reach $b_"max"$ from $b_"min"$ by applying $J_+$ so
$
  b_"max" = b_"min" + n hbar => b_"max" = (n hbar)/2
$
or defining $j equiv b_"max"\/hbar$ we find
$
  j = n/2 => b_"max" = j hbar
$
with $j$ an integer or half-integer. We find
$
  a = hbar^2 j (j + 1)
$
and define
$
  b equiv m hbar
$
with allowed $m$ being
$
  m = -j, -j+1, dots, j-1, j
$
to summarize we have found
$
  bold(J)^2 ket(j","m) & = j(j+1)hbar^2 ket(j","m) \
        J_z ket(j","m) & = m hbar ket(j","m)
$
with $j$ being either integer or half-integer---notably this is a direct result of the commutation relations for $J_k$, which in turn follow from the properties of rotations, and the definition of $J_k$ as the generator of rotations.

=== Matrix elements
We assume $ket(j","m)$ is normalized. Then
$
  braket(j'","m', bold(J)^2, j","m) & = j(j+1)hbar^2 delta_(j'j) delta_(m'm) \
        braket(j'","m', J_z, j","m) & = m hbar delta_(j'j) delta_(m'm)
$
consider
$
  braket(j","m, J_+^dagger J_+, j","m) & = hbar^2 [j(j+1)-m^2-m]
$
we have seen
$
  J_+ ket(j","m) = c^+_(j m) ket(j","m+1)
$
so by comparison
$
  abs(c^+_(j m))^2 & = hbar^2 [j(j+1)-m(m+1)] \
                   & = hbar^2 (j-m)(j+m+1)
$
we obtain
$
  J_+ ket(j","m) & = sqrt((j-m)(j+m+1)) hbar ket(j","m+1) \
  J_- ket(j","m) & = sqrt((j+m)(j-m+1)) hbar ket(j","m-1)
$
so
$
  braket(j'","m', J_plus.minus, j","m) = sqrt((j minus.plus m)(j plus.minus m + 1)) hbar delta_(j'j) delta_(m', m plus.minus 1)
$
Now we can determine the matrix elements of $cal(D)_(m'm)^(j) (R)$ given by
$
  cal(D)_(m'm)^(j) (R) = braket(j","m', exp((-i bold(J) dot hat(bold(n))phi.alt)/hbar), j","m)
$
note we don't care about matrix elements of $cal(D)$ between states with different $j$-values since these vanish trivially. This is since
$
  bold(J)^2 cal(D) ket(j","m) & = cal(D) (R) bold(J)^2 ket(j","m) \
                              & = j(j+1) hbar^2 (cal(D) (R) ket(j","m))
$
so rotations can't change $j$.

The rotations characterized by some $j$ form a group. The identity is just $phi.alt = 0$---corresponding to a $(2j+1) times (2j+1)$ identity matrix. The inverse is given by $phi.alt -> - phi.alt$, and the product of any two members is also a member
$
  sum_m' cal(D)_(m''m')^(j) (R_1) cal(D)_(m' m)^(j) (R_2) = cal(D)_(m''m)^(j) (R_1 R_2)
$
The rotation matrix is unitary since the rotation operator is unitary, explicitly
$
  cal(D)_(m'm) (R^(-1)) = cal(D)_(m m')^* (R)
$
so the conjugate-transpose give the inverse.

Consider
$
  ket(j","m) -> cal(D) (R) ket(j","m)
$
we can expand the rotated state
$
  cal(D) (R) ket(j","m) & = sum_m' ket(j","m') braket(j","m', cal(D)(R), j","m) \
                        & = sum_m' ket(j","m') cal(D)_(m'm)^(j) (R)
$
consider now
$
  cal(D)_(m'm)^(j) (alpha,beta,gamma) &= braket(j","m', exp((-i J_z alpha)/hbar) exp((-i J_y beta)/hbar) exp((-i J_z gamma)/hbar), j","m) \
  &= e^(-i (m' alpha + m gamma)) braket(j","m', exp((-i J_y beta)/hbar), j","m) \
  &= e^(-i (m' alpha + m gamma)) d_(m' m)^(j) (beta)
$
where we define
$
  d_(m'm)^(j) (beta) equiv braket(j","m', exp((-i J_y beta)/hbar), j","m)
$
where we can use
$
  J_y = (J_+ - J_-)/(2 i)
$
consider the case $j = 1$, we have $m,m' = {1,0,-1}$ we can obtain
$
  J_y^(j=1) = hbar/2 mat(0, - sqrt(2) i, 0; sqrt(2) i, 0, - sqrt(2) i; 0, sqrt(2)i, 0)
$
one can show
$
  (J_y^(j=1)/hbar)^3 = J_y^(j=1)/hbar
$
giving
$
  exp((-i J_y beta)/hbar) = 1 - (J_y/hbar)^2 (1 - cos beta) - i (J_y/hbar) sin beta
$
explicitly then
$
  d^1 (beta) = mat(1/2 (1 + cos beta), - 1/sqrt(2) sin beta, 1/2 (1-cos beta); 1/sqrt(2) sin beta, cos beta, -1/sqrt(2) sin beta; 1/2 (1-cos beta), 1/sqrt(2) sin beta, 1/2 (1+cos beta))
$
this quickly becomes hell.

== Orbital Angular Momentum
=== Generator of Rotation
So far we've defined angular momentum as the generator of rotations. Instead if we ignore spin angular momentum, then we can define the angular momentum $bold(J)$ as the orbital angular momentum
$
  bold(L) equiv bold(x) times bold(p)
$
note that this definition implies
$
  [L_i, L_j] = i epsilon_(i j k) hbar L_k
$
as can be seen by just plugging in the definition. Now let
$
  1 - i (dd(phi.alt, d: delta)/hbar) L_z = 1 - i (dd(phi.alt, d: delta)/hbar) (x p_y - y p_x)
$
act on some position eigenket $ket(x'","y'","z')$. Since $x_i$ and $p_j$ commute for $i eq.not j$ and using that momentum is the generator of translations we find
$
  [1 - i (dd(phi.alt, d: delta)/hbar)L_z] ket(x'","y'","z') &= [1 - i (p_y/hbar) (dd(phi.alt, d: delta) x') + i (p_x/hbar) (dd(phi.alt, d: delta) y')] ket(x'","y'"."z') \
  &= ket(x'-y' dd(phi.alt, d: delta)"," y' + x' dd(phi.alt, d: delta)","z')
$
this is what we'd expect for an infinitesimal rotation about $hat(bold(z))$---so if $bold(p)$ generates translation then $bold(L)$ generates rotation. For a wavefunction:
$
  braket(x'","y'","z', [1 - i (dd(phi.alt, d: delta)/hbar) L_z], alpha) = braket(x'+y' dd(phi.alt, d: delta)","y'-x' dd(phi.alt, d: delta)","z', alpha)
$
doing $braket(x'","y'","z', alpha) -> braket(r"," theta"," phi.alt, alpha)$ we find
$
  braket(r","theta","phi.alt, [1- i (dd(phi.alt, d: delta)/hbar) L_z], alpha) &= braket(r","theta","phi.alt - dd(phi.alt, d: delta), alpha) \
  &= braket(r","theta","phi.alt, alpha) - dd(phi.alt, d: delta) pdv(, phi.alt) braket(r","theta","phi.alt, alpha) \
  & => braket(bold(x)', L_z, alpha) = - i hbar pdv(, phi.alt) braket(bold(x)', alpha)
$
as we recognize from familiar wave mechanics. We have a similar thing for rotations about the other axes:
$
  braket(bold(x)', L_x, alpha) &= - i hbar (- sin phi.alt pdv(, theta) - cot theta cos phi.alt pdv(, phi.alt)) braket(bold(x)', alpha) \
  braket(bold(x)', L_y, alpha) &= -i hbar ( cos phi.alt pdv(, theta)- cot theta sin phi.alt pdv(, phi.alt)) braket(bold(x)', alpha)
$
and the ladder operators
$
  braket(bold(x)', L_plus.minus, alpha) = - i hbar e^(plus.minus i phi.alt) (plus.minus i pdv(, theta) - cot theta pdv(, phi.alt)) braket(bold(x)', alpha)
$
then
$
  bold(L)^2 &= L_z^2 + (1/2) (L_+ L_- + L_- L_+) \
  &=> braket(bold(x)', bold(L)^2, alpha) = -hbar^2 [1/(sin^2 theta) pdv(, phi.alt, 2) + 1/(sin theta) pdv(, theta) (sin theta pdv(, theta))] braket(bold(x)', alpha)
$
these are ugly, however we recognize that $bold(L)^2$ is just the angular part of the Laplacian.

=== Spherical Harmonics
The eigenfunctions of a spinless particle in a spherically symmetric potential can be written in the form
$
  braket(bold(x)', n","l","m) = R_(n l) (r) Y_l^m (theta,phi.alt)
$
given the potential only depends on $r$ then the angular dependence is common. We consider
$
  braket(hat(bold(n)), l","m) = Y_l^m (theta, phi.alt) = Y_l^m (hat(bold(n)))
$
also if the Hamiltonian is spherically symmetric then it commutes with $L_z$ and $bold(L)^2$ meaning there exists simultaneous eigenkets---since these satisfy the angular-momentum commutation relation we immediately know the eigenvalues: ${l(l+1)hbar^2, m hbar}$. Consider
$
  L_z ket(l","m) &= m hbar ket(l","m) \
  - i hbar pdv(, phi.alt) braket(hat(bold(n)), l","m) &= m hbar braket(hat(bold(n)), l","m) \
  - i hbar pdv(, phi.alt) Y_l^m (theta,phi.alt) &= m hbar Y_l^m (theta, phi.alt)
$
so the $phi.alt$ dependence is of the form $tilde exp(i m phi.alt)$. Likewise
$
  bold(L)^2 ket(l","m) &= l(l+1)hbar^2 ket(l","m) \
  &=> 0 = [1/(sin theta) pdv(, theta) (sin theta pdv(, theta)) + 1/(sin^2 theta) pdv(, phi.alt, 2) + l(l+1)] Y_l^m
$
this is simply the equation that $Y_l^m$ satisfies. We also have
$
  braket(l'","m', l","m) &= delta_(l l') delta_(m m') \
  &=> integral_0^(2 pi) dd(phi.alt) integral_(-1)^1 dd((cos theta)) Y_(l')^(m'^*) (theta,phi.alt) Y_l^m (theta,phi.alt) ) delta_(l l') delta_(m m')
$
where we used
$
  integral dd(Omega_hat(bold(n))) ketbra(hat(bold(n))) = 1
$
now by considering $L_+ ket(l","l) = 0$ and then succesively using $braket(hat(bold(n)), l","m-1) tilde braket(hat(bold(n)), L_-, l","m)$ one can find
$
  Y_l^m (theta,phi.alt) = (-1)^l/(2^l l!) sqrt((2l+1)/(4 pi) ((1+m)!)/((1-m)!)) e^(i m phi.alt) 1/(sin^m theta) d^(l-m)/(dd((cos theta)^(l-m))) (sin theta)^(2 l)
$
for $m >= 0$ and
$
  Y_l^(-m) (theta, phi.alt) = (-1)^m [Y_l^m (theta, phi.alt)]^*
$
importantly here we don't allow half-integer $l$ and therefore half-integer $m$, if this was allowed we'd get
$
  e^(i m (2 pi)) = -1
$
so after a $2 pi$ rotation we get a $-$ sign. This would make the wavefunction doubly-valued which we don't allow---at least if we define $bold(L) = bold(x) times bold(p)$ to be the generator of rotation---there are also other mathematically motivated reasons.

\* as rotation matrices


== Central Potentials
=== Radial equation
We consider problems of the form
$
  H = bold(p)^2/(2 m) + V (r)",  " r^2 = bold(x)^2
$
since these are all over physics, notably it is spherically symmetric---meaning we have
$
  [bold(L),bold(p)^2] = [bold(L),bold(x)^2] = 0 => [bold(L),H]=[bold(L)^2,H] = 0
$
these problems are referred to as central potentials for the obvious reason.

We have for the energy eigenstates $ket(alpha) = ket(E l m)$:
$
          H ket(E l m) & = E ket(E l m) \
  bold(L)^2 ket(E l m) & = l(l+1) hbar^2 ket(E l m) \
        L_z ket(E l m) & = m hbar ket(E l m)
$
acting with $bra(bold(x)')$ on the first and using the second with
$
  braket(bold(x)', E l m) &= R_(E l) Y_l^m \
  1/(2m) braket(bold(x)', bold(p)^2, alpha) &= - hbar^2/(2m) (pdv(, r, 2) braket(bold(x)', alpha) + 2/r pdv(, r) braket(bold(x)', alpha) - 1/(hbar^2 r^2) braket(bold(x)', bold(L)^2, alpha))
$
gives the radial equation
$
  [- hbar^2/(2 m r^2) dv(, r) (r^2 dv(, r)) + (l(l+1)hbar^2)/(2m r^2) + V(r) ] R_(E l) (r) = E R_(E l) (r)
$
with the full eigenfunction being $R_(E l) Y_l^m$. Now defining $u_(E l) (r) = r R_(E l) (r)$ we find
$
  - hbar^2/(2 m) dv(u_(E l), r, 2) + [(l(l+1)hbar^2)/(2 m r^2) + V (r)] u_(E l) (r) = E u_(E l) (r)
$
with
$
  1 = integral dd(r) r^2 R_(E l)^* (r) R_(E l) (r) = integral dd(r) u_(E l)^* (r) u_(E l) (r)
$
notably our problem has now (as we'd expect) reduced to a one-dimensional problem with the effective potential
$
  V_"eff" (r) = V(r) + (l(l+1)hbar^2)/(2 m r^2)
$

=== Examples
==== Free particle & spherical well
Defining
$
  E equiv (hbar^2 k^2)/(2 m)",  " rho equiv k r
$
we obtain
$
  dv(R, rho, 2) + 2/rho dv(R, rho) + [1 - (l(l+1))/rho^2] R = 0
$
which is famous---with the solutions being the spherical Bessel functions $j_l (rho)$ and $n_l (rho)$, given by
$
  j_l (rho) & = (- rho)^l [1/rho dv(, rho)]^l ((sin rho)/rho) \
  n_l (rho) & = - (- rho)^l [1/rho dv(, rho)]^l ((cos rho)/rho)
$
in the limit $rho -> 0$ we have $j_l (rho) -> rho^l$ and $n_l (rho) -> rho^(-l-1)$, so $n_l (rho)$ diverges at the origin and we discard it here. We can write
$
  j_l (z) = 1/(2 i^l) integral_(-1)^1 dd(s) e^(i z s) P_l (s)
$
this can of course immediately be applied to the spherical well, since we now have $V(r) = 0$ within some $r < a$, so we require the wavefunction be zero at $r = a$. So we require $j_l (k a) = 0$, for $l = 0 => k a = {pi, 2pi, 3 pi, dots}$, so
$
  E_(l=0) = hbar^2/(2 m a^2) {pi^2, (2 pi)^2, (3 pi)^2, dots}
$

==== Isotropic harmonic oscillator
We have
$
  H = bold(p)^2/(2 m) + 1/2 m omega^2 r^2
$
defining
$
  E equiv 1/2 hbar omega lambda",  " r = [hbar/(m omega)]^(1\/2) rho
$
we obtain
$
  dv(u, rho, 2) - (l(l+1))/rho^2 u(rho) + (lambda - rho^2) u(rho) = 0
$
we write $u (rho) = rho^(l+1) e^(- rho^2\/2) f(rho)$ (this removes large and small $rho$ behavior). This gives
$
  rho dv(f, rho, 2) + 2[(l+1)-rho^2)] dv(f, rho) + [lambda - (2 l +3)] rho f (rho) = 0
$
we make the anzats
$
  f (rho) = sum_(n=0)^oo a_n rho^n
$
for order $rho^0$ the only surviving term is $2(l+1)a_1$ so $a_1 = 0$. The $rho^1$ term relates $a_2$ and $a_0$, we get
$
  sum_(n=2)^oo {(n+2)(n+1)a_(n+2) + 2(l+1)(n+2)a_(n+2) - 2 n a_n + [lambda - (2 l +3)]a_n} rho^(n+1) = 0
$
giving
$
  a_(n+2) = (2 n + 2 l + 3 - lambda)/((n+2)(n + 2l + 3)) a_n
$
so $a_1 = 0 => a_"odd" = 0$. For $n -> oo$ we find
$
  a_(n+2)/a_n -> 2/n = 1/q
$
giving
$
  f(rho) -> K times sum_q 1/q! (rho^2)^q prop e^(rho^2)
$
so it seems like this grows forever, this is not normalizable, so we require the series truncates, this requires
$
  2 n + 2l + 3 - lambda = 0
$
for some even $n = 2 q$ giving
$
  E_(q l) = (2 q + l + 3/2) hbar omega equiv (N + 3/2) hbar omega
$
with $N equiv 2 q + l in NN$, giving clear degeneracy. Another way to see this is by writing the Hamiltonian as
$
  H = H_x + H_y + H_z
$
where each $H_i = a_i^dagger a_i + 1\/2$ is an independent harmonic oscillator. Labeling the eigenstates $ket(n_x","n_y","n_z)$ we immediately find
$
  E = (n_x + 1/2 + n_y + 1/2 + n_z + 1/2) hbar omega = (N + 3/2) hbar omega
$
with $N equiv n_x+n_y+n_z$.

==== Coulomb potential
This potential is obviously a big deal,
$
  V (bold(x)) = - (Z e^2)/r
$
so this represents a one-electron atom with atomic number $Z$.

Consider $V(r) -> 0$ for large $r$, then for $r -> oo$ we can write
$
  dv(u_E, r, 2) = kappa^2 u",  " kappa^2 equiv - (2 m E)/hbar^2 > 0
$
given $E < 0$ for bound states. The solution is just $u_E (r) prop e^(- kappa r)$. Similarly considering a $V(r)$ where $lim_(r -> 0) r^2 V(r) = 0$ then for small $r$ we can write
$
  dv(u_(E l), r, 2) = (l(l+1))/r^2 u_(E l) (r)
$
with solution
$
  u_(E l) (r) = A r^(l + 1) + B/r^l
$
we set $B = 0$, since we require normalizability (for $l = 0$ it is more involved but this also leads to issues), so
$
  R_(E l) -> r^l",  " r -> 0
$
now defining $rho equiv kappa r$, and removing both these limits we write
$
  u_(E l) (rho) = rho^(l+1) e^(- rho) w (rho)
$
with $w (rho)$ being nice
$
  dv(w, rho, 2) + 2 ((l+1)/rho -1) dv(w, rho) + [V/E - (2(l+1))/rho] w = 0
$
then one solves for $w ( rho)$ using $V (r = rho\/kappa)$.

Everything above is valid for the Coulomb potential. We define
$
  rho_0 equiv [(2 m)/(-E)]^(1\/2) (Z e^2)/hbar = [(2 m c^2)/(-E)]^(1\/2) Z alpha
$
then
$
  rho dv(w, rho, 2) + 2(l+1-rho) dv(w, rho) + [rho_0 - 2(l+1)] w (rho) = 0
$
this can be written as Kummer's equation
$
  x dv(F, x, 2) + (c-x) dv(F, x) - a F = 0
$
with $x = 2 rho$, $c = 2 (l+1)$ and $2 a = 2(l+1) - rho_0$. The solution is a hypergeometric function written as
$
  F (a;c;x) = 1 + a/c x/1! + (a(a+1))/(c(c+1)) x^2/2! + dots
$
so
$
  w ( rho) = F(l+1-rho_0/2 ; 2(l+1) ; 2 rho)
$
again we can show for large $rho$ that $w (rho) tilde e^rho$ meaning we need this series to terminate. This requires $a + N = 0$ for some $N in NN$. Giving
$
  rho_0 = 2 (N + l + 1)
$
we define $n equiv N + l + 1$ so the energy is
$
  rho_0 = 2 n = [(2 m c^2)/(-E)]^(1\/2) Z alpha => E = - 1/2 m c^2 (Z^2 alpha^2)/n^2 = - 13.6 "eV" Z^2/n^2
$
and
$
  1/kappa = hbar/(m c alpha) n/Z equiv a_0 n/Z
$
giving the relative length scale---with $a_0$ being the Bohr radius. Since the energies only depend on $n$ we can find the degeneracy easily
$
  d = sum_(l=0)^(n-1) (2 l + 1) = n^2
$
we could also write the wavefunction $psi_(n l m) = R_(n l) Y_l^m$ but this is a mess.

== Addition of Angular Momentum
=== Examples
So far we've either ignored internal degrees of freedom (spin) or ignored external degrees of freedom (position, momentum, etc.). We'd like to know how we can account for both. We say the base ket for a spin $1\/2$ particle is in the direct-product space spanned by ${ket(bold(x)')}$ and ${ket(+),ket(-)}$,
$
  ket(bold(x)'","plus.minus) = ket(bold(x)') times.circle ket(plus.minus) in cal(H)_"orb" times.circle cal(H)_"spin"
$
with anything in ${ket(bold(x)')}$ commuting with ${ket(+),ket(-)}$. The rotation operator is the same but with $bold(J) = bold(L) + bold(S) = bold(L) times.circle 1 + 1 times.circle bold(S)$,
$
  cal(D) (R) = cal(D)^(("orb")) (R) times.circle cal(D)^(("spin")) (R) = exp((-i bold(L) dot hat(bold(n))phi.alt)/hbar) times.circle exp((-i bold(S) dot hat(bold(n)) phi.alt)/hbar)
$
since they commute, now by $(A times.circle B) (x times.circle y) = (A x) times.circle (B y)$ it acts like
$
  cal(D) (R) (ket(psi)_"orb" times.circle ket(chi)_"spin") = (cal(D)^(("orb")) (R) ket(psi)_"orb") times.circle (cal(D)^(("spin")) (R) ket(chi)_"spin")
$
The wave function is just
$
  braket(bold(x)'","plus.minus, alpha) = psi_plus.minus (bold(x)')
$
these are often written in a vector.

We can also consider a system of two spin $1\/2$ particles, then the total spin operator is $ bold(S) = bold(S)_1 + bold(S)_2 = bold(S)_1 times.circle 1 + 1 times.circle bold(S)_2 $
again the different spin operators commute, but within each space we still have
$
  [S_(1 x), S_(1 y)] = i hbar S_(1 z)",  " [S_(2 x), S_(2 y)] = i hbar S_(2 z)
$
giving
$
  [S_x,S_y] = i hbar S_z
$
for the total spin operator. We denote the different eigenvalues by
$
  bold(S)^2 & : s (s+1) hbar^2 \
        S_z & : m hbar \
    S_(1 z) & : m_1 hbar \
    S_(2 z) & : m_2 hbar
$
now we can expand the ket of a spin-state in terms of either the eigenkets of $bold(S)^2$ and $S_z$ or of $S_(1 z)$ and $S_(2 z)$.

1. For the second, the ${m_1,m_2}$ representation, we have eigenkets,
$
  {ket(+ +), ket(+ -), ket(- +), ket(- -)}
$
with $m_i = plus.minus 1\/2$.

2. For the first, the ${s,m}$ representation, we have eigenkets,
$
  {ket(s=1","m=plus.minus 1","0), ket(s=0","m=0)}
$
with $s=1$ being the triplet and $s=0$ being the singlet.

In each we have four basekets, they are related as
$
   ket(s=1","m=1) & = ket(+ +) \
   ket(s=1","m=0) & = 1/sqrt(2) (ket(+ -) + ket(- +)) \
  ket(s=1","m=-1) & = ket(- -) \
   ket(s=0","m=0) & = 1/sqrt(2) (ket(+-) - ket(- +))
$
by the RHS of the first both must be spin-up meaning $s = 1, m=1$. Then the second can be obtained by applying $S_-$,
$
  S_- & equiv S_(1 -) + S_(2 -) \
      & = (S_(1 x) - i S_(1 y)) + (S_(2 x) - i S_(2 y))
$
then
$
  S_- ket(s=1"," m=1) &= (S_(1 -) +S_(2 -)) ket(+ +) \
  sqrt((1+1)(1-1+1)) ket(s=1","m=0) &= sqrt((1/2+1/2)(1/2-1/2+1)) (ket(- +) + ket(+ -)) \
$
immediately simplifying to the second. The third can then be obtained by applying $S_-$ to this. The fourth is then obtained by forcing it to be orthogonal to the others. The coefficients of the RHS are called Clebsch-Gordan coefficients and are the elements of the transformation matrix for ${m_1,m_2} -> {s,m}$.

=== Formal theory & Clebsch-Gordan
We consider two angular-momentum operators $bold(J)_1$ and $bold(J)_2$ satisfying the usual commutation relations, and $[J_(1 k), J_(2 l)] = 0$. Then the rotation operator affecting both subspaces is written as
$
  (1 - (i bold(J)_1 dot hat(bold(n)) dd(phi.alt, d: delta))/hbar) times.circle ((1 - (i bold(J)_2 dot hat(bold(n)) dd(phi.alt, d: delta))/hbar)) = 1 - (i (bold(J)_1 times.circle 1 + 1 times.circle bold(J)_2) dot hat(bold(n)) dd(phi.alt, d: delta))/hbar
$
we define the total angular momentum
$
  bold(J) equiv bold(J)_1 times.circle 1 + 1 times.circle bold(J)_2 = bold(J)_1 + bold(J)_2
$
the finite angle rotation is as we've seen
$
  cal(D)_1 (R) times.circle cal(D)_2 (R) = exp((-i bold(J)_1 dot hat(bold(n))phi.alt)/hbar) times.circle exp((-i bold(J)_2 dot hat(bold(n))phi.alt)/hbar)
$
note also that the total angular-momentum satisfies $[J_i,J_j] = i hbar epsilon_(i j k) J_k$, so $bold(J)$ is also an angular-momentum and can be treated as the generator of the entire system. This gives us two options for the choice of basekets:

1. Simultaneous eigenkets of $bold(J)_1^2, bold(J)_2^2, J_(1 z)$, and $J_(2 z)$: $ket(j_1 j_2";" m_1 m_2)$, with defining equations:
$
  bold(J)_1^2 ket(dots) & = j_1(j_1+1) hbar^2 ket(dots) \
      J_(1 z) ket(dots) & = m_1 hbar ket(dots) \
  bold(J)_2^2 ket(dots) & = j_2(j_2+1) hbar^2 ket(dots) \
      J_(2 z) ket(dots) & = m_2 hbar ket(dots)
$

2. Simultaneous eigenkets of $bold(J)^2, bold(J)_1^2, bold(J)_2^2$, and $J_z$. These mutually commute, e.g. $[bold(J)^2, bold(J)_1^2] = 0$ which follows from
$
  bold(J)^2 = bold(J)_1^2 + bold(J)_2^2 + 2 J_(1 z) J_(2 z) + J_(1 +) J_(2 -) + J_(1-) J_(2 +)
$
in this case we write $ket(j_1 j_2 ";" j m)$ with
$
  bold(J)_1^2 ket(dots) & = j_1 (j_1+1) hbar^2 ket(dots) \
  bold(J)_2^2 ket(dots) & = j_2 (j_2+1) hbar^2 ket(dots) \
    bold(J)^2 ket(dots) & = j(j+1) hbar^2 ket(dots) \
          J_z ket(dots) & = m hbar ket(dots)
$
typically this is just written as $ket(j m)$.

Note that importantly $[bold(J)^2, J_(1 z)] eq.not 0 eq.not [bold(J)^2,J_(2 z)]$, which is why we can't add $bold(J)^2$ to the first set of operators---the sets we have constructed are the two maximal sets of mutually compatible observables we can construct.

We can consider a unitary transformation connecting these by
$
  ket(j_1 j_2 ";" j m) = sum_(m_1, m_2) ket(j_1 j_2 ";" m_1 m_2) underbracket(braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m), "Clebsch-Gordan")
$
where we just insert an identity. The matrix elements of the RHS are the Clebsch-Gordan coefficients. These are very important, now we'll go through their properties.

==== Basic properties
If we don't have $m = m_1 + m_2$ then they vanish.
#proof[
  Note that
  $
    (J_z - J_(1 z) - J_(2 z)) ket(j_1 j_2";" j m) = 0
  $
  then acting with $bra(j_1 j_2 ";" m_1 m_2)$ we find
  $
    (m - m_1 - m_2) braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m) = 0
  $
  so unless the Clebsch-Gordan coefficients vanish we have $m = m_1 + m_2$.
]

The coefficients don't vanish unless $abs(j_1 - j_2) <= j <= j_1 + j_2$.

#proof[(partial)][
  We show that if it holds then the dimensionality of the space spanned by ${ket(j_1 j_2 ";" m_1 m_2)}$ is the same as that spanned by ${ket(j_1 j_2 ";" j m)}$. For $(m_1, m_2)$ we have $2 j_1 + 1$ values of $m_1$ and similarly for $j_2$ so $N = (2 j_1 + 1)(2 j_2 + 1)$. For $(j,m)$ each $j$ has $2 j + 1$ states, and by assumption $j_1-j_2 <= j <= j_1 + j_2$ so
  $
    N = sum_(j_1- j_2)^(j_1+j_2) (2 j +1) = (2j_1 + 1)(2 j_2 + 1)
  $
  this shows it is atleast consistent.
]
The Clebsch-Gordan coefficients form a unitary matrix, and we take them to be real by convention. It follows that the inverse coefficients are just the coefficients themselves. Real unitary matrices are orthogonal, so we obtain
$
  sum_(j,m) braket(j_1 j_2 ";" m_1 m_2, j_1 j_2";" j m) braket(j_1 j_2 ";" m'_1 m'_2, j_1 j_2 ";" j m) = delta_(m_1 m'_1) delta_(m_2 m'_2)
$
likewise
$
  sum_(m_1 m_2) braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m) braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j'm') = delta_(j j') delta_(m m')
$
giving also
$
  sum_(m_1, m_2 = m-m_1) abs(braket(j_1 j_2";" m_1 m_2, j_1 j_2 ";" j m))^2 = 1
$
note that there are a million different ways to denote the Clebsch-Gordan coefficients.

==== Recursion relations
We let $j_1, j_2$ and $j$ be fixed. Then consider
$
  J_plus.minus ket(j_1 j_2 ";" j m) &= (J_(1 plus.minus) + J_(2 plus.minus)) sum_(m'_1,m'_2) ket(j_1 j_2 ";" m'_1 m'_2) braket(j_1 j_2 ";" m'_1 m'_2, j_1 j_2 ";" j m) \
  sqrt((j minus.plus m)(j plus.minus m + 1)) ket(j_1 j_2 ";" j","m plus.minus 1) &= sum_(m'_1 m'_2) {sqrt((j_1 minus.plus m'_1)(j_1 plus.minus m'_2 + 1)) ket(j_1 j_2 ";" m'_1 plus.minus 1"," m'_2) \
    &+ sqrt((j_2 minus.plus m'_2)(j_2 plus.minus m'_2 +1)) ket(j_1 j_2 ";" m'_1 "," m'_2 plus.minus 1) } \
  & times braket(j_1 j_2";" m'_1 m'_2, j_1 j_2 ";" j m)
$
now we act with $bra(j_1 j_2";" m_1 m_2)$. By orthonormality the first term vanishes unless $ m_1 = m'_1 plus.minus 1, m_2 = m'_2 $ and the second term vanishes unless $ m_1 = m'_1, m_2 =m'_2 plus.minus 1 $ so we obtain
$
  sqrt((j minus.plus m)(j plus.minus m +1)) & braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j "," m plus.minus 1) \ &= sqrt((j_1 minus.plus m_1 + 1)(j_1 plus.minus m_1)) braket(j_1 j_2 ";" m_1 minus.plus 1 "," m_2, j_1 j_2 ";" j m) \
  &+ sqrt((j_2 minus.plus m_2 + 1)(j_2 plus.minus m_2)) braket(j_1 j_2";" m_1"," m_2 minus.plus 1, j_1 j_2 ";" j m)
$
note the following
$
  abs(m_1) <= j_1",  " abs(m_2) <= j_2",  " -j <= m_1 + m_2 <= j
$
as an example consider adding the orbital and spin-angular momentum of a single spin $1\/2$ particle. We have
$
  j_1 = l " "("integer")",  " m_1 = m_l \
  j_2 = s = 1/2",  " m_2 = m_s = plus.minus 1/2
$
the allowed $j$ are given by $j = l plus.minus 1/2$ for $l > 0$ and $j = 1/2$ for $l = 0$. In $m$-space we'd get two rows, one for each $m_s = plus.minus 1/2$, using $J_-$ recursion we can stay in the upper while $m_l$ changes,
$
  sqrt((l+1/2+m+1)(l+1/2-m)) &braket(m-1/2"," 1/2, l+1/2","m) \
  &= sqrt((l+1+1/2)(l-m+1/2)) braket(m+1/2","1/2, l+1/2","m+1)
$
where we don't write $j_1 = l, j_2 = 1/2$ and used $m_1 = m_l = m-1/2, m_2 = m_s = 1/2$. So we obtain
$
  braket(m-1/2","1/2, l+1/2","m) = sqrt((l+m+1/2)/(l+m+3/2)) braket(m+1/2","1/2, l+1/2","m+1)
$
which is a mess. This can be done again
$
  braket(m-1/2","1/2, l+1/2","m) &= sqrt((l+m+1/2)/(l+m+3/2)) sqrt((l+m+3/2)/(l+m+5/2)) braket(m+3/2","1/2, l+1/2","m+2) \
  &= dots \
  &= sqrt((l+m+1/2)/(2 l + 1)) braket(l","1/2, l + 1/2","l+1/2)
$
consider the maximal case $m_l = l$ and $m_s = 1/2$. then $m = l + 1/2$ which is only possible for $j = l + 1/2$. So $ket(m_l = l","m_s = 1/2)$ is equal to $ket(j=l+1/2","m=l+1/2)$, up to a phase. This gives
$
  braket(m-1/2","1/2, l+1/2","m) = sqrt((l+m+1/2)/(2 l + 1))
$
this gives one matrix element of the transformation matrix $(m_l, m_s) -> (j,m)$, which has the form
$
  mat(cos alpha, sin alpha; - sin alpha, cos alpha)
$
(due to orthogonality) by comparison what we just found is $cos alpha$, then
$
  sin^2 alpha = 1 - cos^2 alpha = (l-m + 1/2)/(2 l +1)
$
up to a sign, we claim it is positive since all states are reachable by applying $J_-$ succesively, and the matrix elements of $J_-$ are positive by convention. So we find
$
  mat(sqrt((l+m+1/2)/(2 l +1)), sqrt((l-m+1/2)/(2 l +1)); - sqrt((l-m+1/2)/(2 l +1)), sqrt((l+m+1/2)/(2 l + 1)))
$
We define spin-angular functions as follows
$
  cal(Y)_l^(j=l plus.minus 1\/2, m) = 1/sqrt(2 l +1) vec(plus.minus sqrt(l plus.minus m + 1/2) Y_l^(m-1\/2) (theta,phi.alt), sqrt(l minus.plus m + 1/2) Y_l^(m+1\/2) (theta,phi.alt))
$
these are by construction simultaneous eigenfunctions of $bold(L)^2, bold(S)^2, bold(J)^2$, and $J_z$.

=== Rotation matrices
Let $cal(D)^((j_1)) (R)$ be the rotation operator in the ket space spanned by the eigenket with eigenvalue $j_1$. Then we can write the Clebsch-Gordan series
$
  cal(D)^((j_1))_(m_1 m'_1) (R) cal(D)^((j_2))_(m_2 m'_2) (R) &= sum_(j,m,m') braket(j_1 j_2 ";" m_1 m_2, j_1 j_2";" j m) \
  &times braket(j_1 j_2";"m'_1 m'_2, j_1 j_2 ";" j m') cal(D)_(m m')^((j)) (R)
$
with $j: abs(j_1 -j_2) -> j_1 + j_2$. This comes from $cal(D)^((j_1)) times.circle cal(D)^((j_2))$ which is reducible.
#proof[
  The LHS is equivalent to
  $
    braket(j_1 j_2 ";" m_1 m_2, cal(D) (R), j_1 j_2 ";" m'_1 m'_2) &= braket(j_1 m_1, cal(D) (R), j_1 m'_1) braket(j_2 m_2, cal(D) (R), j_2 m'_2) \
    &= cal(D)^((j_1))_(m_1 m'_1) (R) cal(D)^((j_2))_(m_2 m'_2) (R)
  $
  but we can also insert a couple of completion relations
  $
    braket(j_1 j_2 ";" m_1 m_2, cal(D) (R), j_1 j_2 ";" m'_1 m'_2) &= sum_(j m,j'm') braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m) braket(j_1 j_2 ";" j m, cal(D) (R), j_1 j_2 ";" j' m') \
    & times braket(j_1 j_2";"j' m', j_1 j_2 ";" m'_1 m'_2) \
    &= sum_(j m, j' m') braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m) cal(D)^((j))_(m m') (R) delta_(j j') \ & times braket(j_1 j_2 ";" m'_1 m'_2, j_1 j_2 ";" j' m')
  $
  which is the RHS.
]

== Schwinger's model
We consider two simple harmonic oscillators, we denote these by the plus type and the minus type. We have accompanying annihilation and creation operators $a_+$ and $a_+^dagger$ etc. Similarly we define
$
  N_+ equiv a_+^dagger a_+",  " N_- equiv a_-^dagger a_-
$
we assume the usual commutation relations hold for each, e.g.
$
   [a_+,a_+^dagger] & = 1 \
          [N_+,a_+] & = - a_+ \
  [N_+, a_+^dagger] & = a_+^dagger
$
and similarly for $a_-$, but we assume they commute with eachother; $[a_+,a_-^dagger]=[a_-,a_+^dagger]=0$ etc. We call the oscillators decoupled, or independent. $N_-$ and $N_+$ commute, so there exists a simultaneous eigenket with eigenvalues $n_+$ and $n_-$. So we have
$
  N_+ ket(n_+","n_-) = n_+ ket(n_+","n_-)",   " N_- ket(n_+","n_-) = n_- ket(n_+","n_-)
$
likewise we have
$
  a_+^dagger ket(n_+","n_-) = sqrt(n_++1) ket(n_+ +1","n_-)",  " a_+ ket(n_+","n_-) = sqrt(n_+) ket(n_+-1","n_-)
$
and similarly for $a_-$ and $a_-^dagger$. From
$
  a_+ ket(0","0) = 0",  " a_- ket(0","0) = 0
$
we can obtain, by succesive application
$
  ket(n_+","n_-) = ((a_+^dagger)^(n_+) (a_-^dagger)^(n_-))/(sqrt(n_+ !) sqrt(n_- !)) ket(0","0)
$
now we define
$
  J_+ equiv hbar a_+^dagger a_-",  " J_- equiv hbar a_-^dagger a_+
$
and
$
  J_z equiv hbar/2 (a_+^dagger a_+ - a_-^dagger a_-) = hbar/2 (N_+ - N_-)
$
these satisfy
$
  [J_z,J_plus.minus] & = plus.minus hbar J_plus.minus \
           [J_+,J_-] & = 2 hbar J_z
$
#proof[
  For the second we have
  $
    hbar^2 [a_+^dagger a_-, a_-^dagger a_+] &= hbar^2 a_+^dagger a_- a_-^dagger a_+ - hbar^2 a_-^dagger a_+ a_+^dagger a_- \
    &= hbar^2 a_+^dagger (a_-^dagger a_- + 1) a_+ - hbar^2 a_-^dagger (a_+^dagger a_+ + 1) a_- \
    &= hbar^2 (a_+^dagger a_+ - a_-^dagger a_-) = 2 hbar J_z
  $
]
defining $N equiv N_+ + N_-$ we can write
$
  bold(J)^2 & equiv J_z^2 + 1/2 (J_+ J_- + J_- J_+) \
            & = (hbar^2/2) N (N/2 + 1)
$

We associate spin up with one quantum unit of the plus-type, and spin down with one quantum unit of the minus-type. $n_+$ and $n_-$ are then the amount of spin ups and downs respectively. Then $J_+$ destroys one unit of spin down with $z$-component $-hbar\/2$ and creates one unit of spin up with $z$-component $hbar\/2$, so it increases by $hbar$. A similar interpretation holds for $J_-$, and $J_z$ evidently just counts $hbar\/2$ times the difference of $n_+$ and $n_-$.

We can find
$
  J_+ ket(n_+","n_-) & = sqrt(n_- (n_++1)) hbar ket(n_+ + 1"," n_- -1) \
  J_- ket(n_+","n_-) & = sqrt(n_+ (n_- + 1)) hbar ket(n_+ -1"," n_- + 1) \
  J_z ket(n_+","n_-) & = 1/2 (n_+-n_-) hbar ket(n_+","n_-)
$
note that $n_+ + n_-$ is unchanged, and if we substitute $n_+ -> j + m$ and $n_- -> j - m$ then these reduce to the relations we are used to. And the eigenvalue of $bold(J)^2$ becomes $ hbar^2 j(j+1) $
this is of course expected since we constructed $J_plus.minus$ and $J_z$ such that they satisfy the same commutation relations. We can then write
$
  j equiv (n_+ + n_-)/2",  " m equiv (n_+ - n_-)/2
$
then the action of $J_+$ changes $n_+ -> n_+ + 1$ and $n_- -> n_- - 1$ leaving $j$ unchanged and $m -> m + 1$, similarly $J_-$ leaves $j$ unchanged and $m -> m - 1$, as we'd expect. We can also write the general eigenket as
$
  ket(j","m) = ((a_+^dagger)^(j+m) (a_-^dagger)^(j-m))/sqrt((j+m)!(j-m)!) ket(0)
$
letting $j = m$ we find the special case
$
  ket(j","j) = ((a_+^dagger)^(2 j))/sqrt((2j)!) ket(0)
$
which corresponds to the physical case where the eigenvalue of $J_z$ is as large as possible for a given $j$. Evidently we can treat this state like $2 j$ spin-$1\/2$ particles with all spins pointing up. In general some object with high $j$ can be treated as being made up from $j+m$ spin up particles and $j-m$ spin down particles. This is convenient but in no way physical, just that the object transforms in this way.

=== Rotation matrices
We apply $cal(D) (R)$ to $ket(j","m)$, the only non-trivial rotation is the second one (euler notation) so we just care about
$
  cal(D) (R) = cal(D) (alpha,beta,gamma)_(alpha = gamma = 0) = exp((-i J_y beta)/hbar)
$
then
$
  cal(D) (R) ket(j","m) = ([cal(D) (R) a_+^dagger cal(D)^(-1) (R)]^(j+m) [cal(D) (R) a_-^dagger cal(D)^(-1) (R)]^(j-m))/sqrt((j+m)! (j-m)!) cal(D) (R) ket(0)
$
now $cal(D) (R) ket(0) = ket(0)$ and
$
  cal(D) (R) a_plus.minus^dagger cal(D)^(-1) (R) = exp((-i J_y beta)/hbar) a_plus.minus^dagger exp((i J_y beta)/hbar)
$
then by Baker-Hausdorff
$
  cal(D) (R) a_plus.minus^dagger cal(D)^(-1) (R) = a_plus.minus^dagger cos beta/2 plus.minus a_minus.plus^dagger sin beta/2
$
then after using the binomial theorem we finally obtain
$
  cal(D) (alpha=0,beta, gamma=0) ket(j","m) &= sum_(k,l) ((j+m)!(j-m)!)/((j+m-k)! k! (j-m-l)! l !) \
  &times ([a_plus^dagger cos beta\/2]^(j+m-k) [a_-^dagger sin beta\/2]^k)/sqrt((j+m)!(j-m)!) \
  &times [-a_plus^dagger sin beta\/2]^(j-m-l) [a_-^dagger cos beta\/2]^l ket(0)
$
we can compare this with our previous result for the same thing
$
  cal(D) (alpha=0,beta,gamma=0) ket(j","m) &= sum_m' ket(j","m') d_(m'm)^(j) (beta) \
  &= sum_m' d_(m'm)^j (beta) ((a_+^dagger)^(j+m') (a_-^dagger)^(j-m'))/sqrt((j+m')! (j-m')!) ket(0)
$
we can find an expression for the Wigner functions $d_(m' m)^j (beta)$ by equating powers of $a_+^dagger$. We want to compare $j+m'$ and $2j-k-l$, giving $l = j-k-m'$. The $k$ and $l$ sums are not independent, so we eliminate $l$ in favor of $k$---note that the powers of $a_-^dagger$ match with the same constraint. Lastly we look at the exponents of $cos beta\/2$ and $sin beta\/2$ giving
$
  j+m-k+l & = 2j- 2k +m - m' \
  k+j-m-l & = 2k -m +m' \
    j-m-l & = k-m+m'
$
so we obtain
$
  d^j_(m'm) (beta) &= sum_k (-1)^(k-m+m') sqrt((j+m)! (j-m)! (j+m')! (j-m')!)/((j+m-k)! k! (j-k-m')! (k-m+m')!) \
  &times (cos beta/2)^(2j-2k+m-m') (sin beta/2)^(2k-m+m')
$
we sum over $k$ when none of the factorials in the denominator become negative.

=== A detour
We'll briefly discuss Landau levels. The Hamiltonian in the presence of a magnetic field is given by
$
  H = 1/(2 m) (bold(p)- e/c bold(A))^2 "with" bold(B) = nabla times bold(A)
$
and the equation of motion is
$
  m dot.double(bold(x)) = e/c dot(bold(x)) times bold(B)
$
in the case of $bold(B) = B_0 hat(bold(z))$. We use the symmetric gauge
$
  bold(A) = - 1/2 bold(x) times bold(B) = B_0/2 vec(-y, x, 0)
$
then
$
  H &= 1/(2 m) vec(p_x + (e B_0)/(2 c) y, p_y - (e B_0)/(2 c) x)^2 \
  &= underbrace(1/(2 m) (p_x^2 + p_y^2) + 1/2 m (omega_L/2)^2 (x^2+y^2), "2D SHO with" omega= omega_L\/2) - omega_L/2 underbrace((x p_y - y p_x), L_z)
$
by Schwinger
$
  x_i = sqrt(hbar/(m omega_L)) (a_i + a_i^dagger)
$
similarly for $p_i$. Then we define
$
  a_plus.minus = 1/sqrt(2) (a_1 minus.plus i a_2);"  " [a_+,a_-] = 0 dots
$
and $[a_+,a_+^dagger] = 1 dots$, so these are uncoupled, and
$
  L_z = hbar(a_+^dagger a_+ - a_-^dagger a_-)
$
so we can write
$
  H &= hbar (omega_L)/2 (N_+ + N_- + 1) - hbar omega_L/2 (a_+^dagger a_+ - a_-^dagger a_-) \
  &= hbar omega_L (a_-^dagger a_- + 1/2) -> E_n = hbar omega_L (n_- + 1/2)
$
and $[H,L_z] = 0$ with eigenket $ket(n_+","n_-)$ so
$
  L_z = hbar (n_+-n_-)
$
the found energies are the Landau levels. These have infinite degeneracy since to each $n_-$ we have infinite $n_+$ in principle, and thereby different $L_z$.


#pagebreak()
== Bell's inequality
We consider a two-electron system in a spin-singlet state, so with total spin zero. We can write this state as
$
  ket("singlet") = 1/sqrt(2) (ket(hat(bold(z))+";"hat(bold(z))-)- ket(hat(bold(z))-";"hat(bold(z))+))
$
with the obvious interpretation. In this case the two electrons are correlated in the sense that if we measure the spin of one electron and collapse to e.g. $ket(hat(bold(z))+";"hat(bold(z))-)$ then the second electron is forced to take some value. This correlation happens even at distance, a relevant example is the decay $eta -> mu^+ + mu^-$ where $eta$ has spin zero and each $mu^(plus.minus)$ has spin half. We could imagine measuring the spin of $mu^+$ at some detector A and measureing the spin of $mu^-$ at some detector B. Then the measurements would be correlated in the obvious manner, even at macroscopic distance.

Recall for this type of system
$
  ket(hat(bold(x)) plus.minus) = 1/sqrt(2) (ket(hat(bold(z))+) plus.minus ket(hat(bold(z))-))
$
and similarly with $hat(bold(z))<--> hat(bold(x))$. Then one can show by substitution that
$
  ket("singlet") = 1/sqrt(2) (ket(hat(bold(x))-";"hat(bold(x))+)-ket(hat(bold(x))+";"hat(bold(x))-))
$
which is the same apart form an overall phase. Now we can consider measuring $S_x$ and/or $S_z$ at A and B. If we measure $S_z$ at A and $S_x$ at B then the measurements are randomly correlated. If we measure the same at both A and B then we find total anti-correlation. Lastly if we do no measurement at A then the measurement at B is random. It seems like the second particle knows which spin component of the first particle is being measured instantaneously. This seems to contradict Einstein's locality principle and the problem itself is known as the Einstein-Podolsky-Rosen paradox. Some ascribe the aforementioned measurement behaviour to be a result of hidden variables. We'd like to know if these hidden variable theories make predictions that can distinguish them from quantum mechanics---in this manner we can try to find some true theory. This is what Bell's inequality does.

We derive the Bell's inequality using a framework by Wigner. We consider an ensemble of spin half particles with zero total angular momentum, now we define different species of particles. So some have $(hat(bold(z))+,hat(bold(x))-)$ while other have $(hat(bold(z))-,hat(bold(x))-)$ etc. To ensure the angular momentum vanishes any particular pair must match so if particle 1 has $(hat(bold(z))+,hat(bold(x))-)$ then particle 2 must have $(hat(bold(z))-,hat(bold(x))+)$. For two particles we have four possible pairs with equal probability $1/4$, this leads to the same measurement correlations as the previous. Now we consider the arbitrary unit vectors ${hat(bold(a)), hat(bold(b)), hat(bold(c))}$ then we have e.g. the type $[hat(bold(a))-,hat(bold(b))+,hat(bold(c))+]$ etc. and again we must necessarily have matching pairs, now we just have 8 types with each type having a population $N_i$.

Now suppose A measures $bold(S)_1 dot hat(bold(a))$ to be $+$ and B measures $bold(S)_2 dot hat(bold(b))$ to be $+$ then this pair belongs to either $N_3$ or $N_4$. Note that
$
  N_3 + N_4 <= (N_2 +N_4) + (N_3 + N_7)
$
then we have
$
  P(hat(bold(a))+";"hat(bold(b))+) = (N_3 +N_4)/(sum N_i)
$
with the obvious meaning, similarly
$
  P(hat(bold(a))+";"hat(bold(c))+) & = (N_2 + N_4)/(sum N_i) \
  P(hat(bold(c))+";"hat(bold(b))+) & = (N_3+N_7)/(sum N_i)
$
so we immediately find the Bell's inequality
$
  P(hat(bold(a))+";"hat(bold(b))+) <= P(hat(bold(a))+";"hat(bold(c))+) + P(hat(bold(c))+";"hat(bold(b))+)
$
so for hidden variable theories this inequality is satisfied.

We now show that this inequality is violated in quantum mechanics. We characterize all singlet states by $ket("singlet")$ defined before, then we can calculate all terms in the Bell's inequality. We define $bold(hat(a)) dot bold(hat(b)) = cos theta_"ab"$ and take both to be in the $x y$-plane with $bold(hat(a))$ being $bold(hat(x))$ without loss of generality. Then
$
  ket(S_a";"plus.minus) equiv ket(S_x";"plus.minus) = 1/sqrt(2) (ketbra(+)-ketbra(-))
$
we rotate this by $theta_"ab"$ about the $z$-axis then
$
  ket(S_b";"+) &= exp((-i S_z theta_"ab")/hbar) ket(S_x";"+) \
  &= exp((-i theta_"ab")/2) 1/sqrt(2) ket(+) + exp((i theta_"ab")/2) 1/sqrt(2) ket(-)
$
then
$
  braket(S_a";"+, S_b";"+) & = cos theta_"ab"/2
$
so
$
  P(S_a +";" S_b +) = 1/2 cos^2 theta_"ab"/2
$
with the factor $1/2$ coming from $S_a$ being $+$. So we find
$
  sin^2 theta_"ab"/2 <= sin^2 theta_"ac"/2 + sin^2 theta_"cb"/2
$
this should hold for all angles for convenience pick $theta_"ab" = 2 theta$ and $theta_"ac" = theta_"bc" = theta$ then it is violated for $0 < theta < pi\/2$. This was verified experimentally by Aspect et al.
