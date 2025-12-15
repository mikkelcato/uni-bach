//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Theory of Angular Momentum
Rotations about the same axis commute and we want to know how rotations about different axes fail. Rotations can be represented by $3 times 3$ orthogonal matrices $R$. An infinitesimal rotation is of the form
$
  R_z (epsilon) = mat(1-epsilon^2/2, -epsilon, 0; epsilon, 1-epsilon^2/2, 0; 0, 0, 1)
$
and similarly for $x$ and $y$. To order $epsilon^2$ we can find
$
  R_x (epsilon) R_y (epsilon) - R_y (epsilon) R_x (epsilon) = R_z (epsilon^2) - R_"any" (0)
$
with $R_"any" (0) = 1$. We use this later.

== Rotation operator
We define the rotation operator $D (R)$ by
$
  ket(alpha)_R = D (R) ket(alpha)
$
with $ket(alpha)_R$ being the rotated state.

We propose an operator of the form
$
  U_epsilon = 1 - i G epsilon
$
and classically angular momentum generates rotations so we define the angular momentum operator $J_k$ such that
$
  G = J_k/hbar " makes " epsilon -> dd(phi.alt)
$
with $J_k$ Hermitian the operator is guaranteed to be unitary and reduces to $bb(1)$ in the limit $dd(phi.alt) -> 0$. We obtain
$
  D (hat(bold(n)),dd(phi.alt)) = 1 - i ((bold(J)dot hat(bold(n)))/hbar) dd(phi.alt)
$
for a rotation about $hat(bold(n))$ by $dd(phi.alt)$.

Let $hat(bold(n)) = hat(bold(z))$. Then a finite rotation is given by
$
  D_z (phi.alt) & = exp((-i J_z phi.alt)/hbar)
$
We assume $D (R)$ has the same group properties as $R$. Then the classical commutator outlined in the beginning eventually gives
$
  [J_i, J_j] = i hbar epsilon_(i j k) J_k
$
which is very important!

== Spin $1/2$ system
Recall
$
  S_x & = hbar/2 (ketbra(+, -) + ketbra(-, +)) \
  S_y & = (i hbar)/2 (-ketbra(+, -) + ketbra(-, +)) \
  S_z & = hbar/2 (ketbra(+, +) - ketbra(-, -))
$
these satisfy the commutation relations with $J_k -> S_k$.

By a short computation $expval(A) -> expval(D^dagger A D)$ one can show
$
  expval(S_x) & -> expval(S_x) cos phi.alt - expval(S_y) sin phi.alt \
  expval(S_y) & -> expval(S_y) cos phi.alt + expval(S_x) sin phi.alt
$
and since $[S_z, D_z (phi.alt)] = 0$
$
  expval(S_z) -> expval(S_z)
$
this shows that expectation values behave like vectors under rotation
$
  expval(J_k) -> sum_l R_(k l) expval(J_l)
$
Consider
$
  ket(alpha) = ket(+) braket(+, alpha) + ket(-) braket(-, alpha)
$
we find
$
  exp((-i S_z phi.alt)/hbar) ket(alpha) = e^(-i phi.alt \/2) ket(+) braket(+, alpha) + e^(i phi.alt\/2) ket(-) braket(-, alpha)
$
for $phi.alt = 2 pi$ we find $ket(alpha) -> - ket(alpha)$!

We previously found
$
  U (t,0) = exp((-i S_z omega t)/hbar)
$
this is just the rotation operator with $phi.alt = omega t$! So we immediately obtain
$
  expval(S_x)_t = expval(S_x)_(t = 0) cos omega t - expval(S_y)_(t=0) sin omega t
$
etc. so we see precession! The state is given by
$
  ket(alpha\,t_0 = 0\;t) = e^(-i omega t\/2) ket(+) braket(+, alpha) + e^(i omega t\/2) ket(-) braket(-, alpha)
$
and it is worth noting that
$
  tau_"precession" = (2 pi)/omega" and" tau_"state" = (4 pi)/omega
$
with the _extra_ minus sign being a real thing!

== Pauli formalism
We define
$
  ket(+) eq^dot vec(1, 0) equiv chi_+
$
and similarly for the spin down state. Then
$
  ket(alpha) eq^dot vec(braket(+, alpha), braket(-, alpha)) = chi equiv (c_+,c_-)
$
with $chi$ being the two-component spinor. We write $S_k$ in terms of the Pauli matrices $sigma_k$
$
  braket(plus.minus, S_k, +) equiv hbar/2 (sigma_k)_(plus.minus, plus)" and " braket(plus.minus, S_k, -) equiv hbar/2 (sigma_k)_(plus.minus,-)
$
so
$
  expval(S_k) = sum_(a',a'' = plus,minus) braket(alpha, a') braket(a', S_k, a'') braket(a'', alpha) = hbar/2 chi^dagger sigma_k chi
$
explicitly
$
  sigma_1 = mat(0, 1; 1, 0)";  " sigma_2 = mat(0, -i; i, 0)";  " sigma_3 = mat(1, 0; 0, -1)
$
We can find
$
  {sigma_i, sigma_j} & = 2 delta_(i j)";  " [sigma_i, sigma_j] = 2 i epsilon_(i j k) sigma_k
$
and many other nice properties!

Consider $bold(sigma) dot bold(a)$
$
  bold(sigma) dot bold(a) equiv sum_k a_k sigma_k = mat(a_3, a_1 -i a_2; a_1+i a_2, - a_3)
$
with the identity
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

Then we can write the rotation operator as
$
  D (hat(bold(n)),phi.alt) = exp((-i bold(S) dot hat(bold(n)) phi.alt)/hbar) eq^dot underbracket(exp((-i bold(sigma) dot hat(bold(n))phi.alt)/2), "spinor-representation")
$
the above identity gives
$
  (bold(sigma) dot hat(bold(n)))^n = cases(
    1 #h(40pt) & "for" n "even",
    bold(sigma) dot hat(bold(n)) & "for" n "odd"
  )
$
so we can write
$
  exp((-i bold(sigma) dot hat(bold(n)) phi.alt)/2) &= bb(1) cos phi.alt/2 - i bold(sigma) dot hat(bold(n)) sin phi.alt/2
$
This acts on the two-component spinor $chi$ by
$
  chi -> exp((-i bold(sigma)dot hat(bold(n))phi.alt)/2) chi
$

== Groups
We denote the group of $3 times 3$ orthogonal matrices with determinant $1$ by $"SO"(3)$.

We can also represent rotations by unitary $2 times 2$ matrices acting on two-component spinors. As we saw above with
$
  D =^dot exp((-i bold(sigma) dot hat(bold(n)) phi.alt)/2)
$
acting on $chi$. A general $M in "SU"(2)$ can be written as
$
  U (a,b) = mat(a, b; -b^*, a^*)
$
by comparison we can find
$
  Re a & = cos phi.alt/2",  " Im a = -n_z sin phi.alt/2 \
  Re b & = -n_y sin phi.alt/2",  " Im b = -n_x sin phi.alt/2
$
the $a$ and $b$ are known as Cayley-Klein parameters. A general member of $"U"(2)$ can be written as
$
  U = e^(i gamma) mat(a, b; -b^*, a^*)
$
with $abs(a)^2+abs(b)^2=1$ and $gamma = gamma^*$. All these groups are subgroups of $"GL"_n (CC)$.

Classically any rotation can be accomplished using three rotations and Euler angles by
$
  R(alpha,beta,gamma) = R_z (alpha) R_y (beta) R_z (gamma)
$
Applying this to spin $1/2$ we find
$
  D (alpha,beta,gamma) = D_z (alpha) D_y (beta) D_z (gamma)
$
which has the representation
$
  D (alpha,beta,gamma) =^dot mat(e^(-i(alpha+gamma)\/2) cos beta/2, -e^(-i(alpha-gamma)\/2) sin beta/2; e^(i(alpha-gamma)\/2) sin beta/2, e^(i(alpha+gamma)\/2) cos beta/2)
$
This has the form of $U(a,b)$ and is called the $j = 1/2$ irreducible representation of $D (alpha,beta,gamma)$. This guy has matrix elements
$
  D_(m'm)^(1\/2) (alpha,beta,gamma) = braket(j=1/2\,m', exp((-i J_z alpha)/hbar) exp((-i J_y beta)/hbar) exp((-i J_z gamma)/hbar), j =1/2\,m)
$

== Eigenvalues and eigenstates
We define
$
  bold(J)^2 equiv J_x J_x + J_y J_y + J_z J_z
$
this commutes with all $J_k$
$
  [bold(J)^2, J_k] = 0
$
Then we can find simultaneous eigenstates for $bold(J)^2$ and $J_z$. We write
$
  bold(J)^2 ket(a\,b) = a ket(a\,b)",  " J_z ket(a\,b) = b ket(a\,b)
$
To determine $a, b$ we define ladder operators
$
  J_plus.minus equiv J_x plus.minus i J_y
$
these satisfy
$
  [J_+ ,J_-] = 2 hbar J_z",  " [J_z,J_plus.minus] = plus.minus hbar J_z
$
which can be easily proven. Also $[bold(J)^2, J_plus.minus] = 0$. Consider
$
  J_z (J_plus.minus ket(a\,b)) &= ([J_z, J_plus.minus]+J_plus.minus J_z) ket(a\,b) \
  &= (b plus.minus hbar) (J_plus.minus ket(a\,b))
$
so applying $J_plus.minus$ to a $J_z$ eigenket gives a new eigenket of $J_z$ with eigenvalue $b plus.minus hbar$ and
$
  bold(J)^2 (J_plus.minus ket(a\,b)) & = J_plus.minus bold(J)^2 ket(a\,b) \
                                     & = a (J_plus.minus ket(a\,b))
$
So $J_plus.minus ket(a\,b)$ are simultaneous eigenkets of $bold(J)^2$ and $J_z$ with eigenvalues $a$ and $b plus.minus hbar$. Then we can write
$
  J_plus.minus ket(a","b) = c_plus.minus ket(a"," b plus.minus hbar)
$

Consider applying $J_+$ $n$ times to a simultaneous eigenket. Then we obtain an eigenket with the $J_z$ eigenvalue increased by $n hbar$ and the $bold(J)^2$ eigenvalue unchanged. This process has an upper limit by
$
  a >= b^2
$

#proof[
  To see this note
  $
    bold(J)^2 - J_z^2 & = 1/2 (J_+ J_- + J_- J_+) = 1/2 (J_+ J_+^dagger + J_+^dagger J_+)
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

Then there is some $b_"max"$ where
$
  J_+ ket(a","b_"max") & = 0 \
                       & => J_- J_+ ket(a","b_"max") = 0
$
with
$
  J_- J_+ & = J_x^2 + J_y^2 - i (J_y J_x - J_x J_y) \
          & = bold(J)^2 - J_z^2 - hbar J_z
$
this implies
$
  (bold(J)^2 - J_z^2 - hbar J_z) ket(a","b_"max") = 0
$
We obtain
$
  a = b_"max" (b_"max" + hbar)
$
and similarly
$
  a = b_"min" (b_"min" - hbar)
$
by comparison $b_"max" = - b_"min" > 0$ so the allowed values of $b$ lie in
$
  -b_"max" <= b <= b_"max"
$
We should be able to reach $b_"max"$ from $b_"min"$ by applying $J_+$ so
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
We have found
$
  bold(J)^2 ket(j\,m) & = j(j+1)hbar^2 ket(j\,m) \
        J_z ket(j\,m) & = m hbar ket(j\,m)
$
with $j$ being either integer or half-integer. Which is quite nice!

We assume $ket(j\,m)$ is normalized. Then
$
  braket(j'\,m', bold(J)^2, j\,m) & = j(j+1)hbar^2 delta_(j'j) delta_(m'm) \
        braket(j'\,m', J_z, j\,m) & = m hbar delta_(j'j) delta_(m'm)
$
Consider
$
  braket(j\,m, J_+^dagger J_+, j\,m) & = hbar^2 {j(j+1)-m^2-m}
$
we have seen
$
  J_+ ket(j\,m) = c^+_(j m) ket(j\,m+1)
$
so by comparison
$
  abs(c^+_(j m))^2 & = hbar^2 [j(j+1)-m(m+1)] \
                   & = hbar^2 (j-m)(j+m+1)
$
and we obtain
$
  J_+ ket(j\,m) & = sqrt((j-m)(j+m+1)) hbar ket(j\,m+1) \
  J_- ket(j\,m) & = sqrt((j+m)(j-m+1)) hbar ket(j\,m-1)
$
so
$
  braket(j'\,m', J_plus.minus, j\,m) = sqrt((j minus.plus m)(j plus.minus m + 1)) hbar delta_(j'j) delta_(m', m plus.minus 1)
$
We can also determine the matrix elements of $D_(m'm)^(j) (R)$
$
  D_(m'm)^(j) (R) = braket(j\,m', exp((-i bold(J) dot hat(bold(n))phi.alt)/hbar), j\,m)
$
rotations do not change $j$ since $[D,bold(J)^2]=0$ meaning
$
  bold(J)^2 D ket(j\,m) & = D (R) bold(J)^2 ket(j\,m) \
                        & = j(j+1) hbar^2 (D (R) ket(j\,m))
$
All values of $j$ define their own group. The identity is $phi.alt = 0$ and the inverse is given by $phi.alt -> - phi.alt$.

Consider
$
  ket(j\,m) -> D (R) ket(j\,m)
$
we expand the rotated state
$
  D (R) ket(j\,m) = sum_m' ket(j\,m') braket(j\,m', D(R), j\,m) = sum_m' ket(j\,m') D_(m'm)^(j) (R)
$
then we could write $D_(m' m)^(j)$ in terms of Euler angles but this is ass. We would find it fully determined by
$
  d_(m'm)^(j) (beta) equiv braket(j","m', exp((-i J_y beta)/hbar), j","m)
$

== Orbital Angular Momentum
We define
$
  bold(L) = bold(x) times bold(p)
$
which generates rotations if $bold(p)$ generates translations.

$
  dots #h(2em) ("see Sakurai")
$

The eigenfunctions of a spinless particle in a spherically symmetric potential can be written as
$
  braket(bold(x)', n\,l\,m) = R_(n l) (r) underbracket(Y_l^m (theta,phi.alt), "spherical harmonics")
$
we consider
$
  braket(hat(bold(n)), l\,m) = Y_l^m (theta, phi.alt) = Y_l^m (hat(bold(n)))
$
Assuming the Hamiltonian is spherically symmetric then it commutes with $L_z$ and $bold(L)^2$. Then
$
  L_z ket(l","m) &= m hbar ket(l","m) \
  - i hbar pdv(, phi.alt) braket(hat(bold(n)), l","m) &= m hbar braket(hat(bold(n)), l","m) \
  - i hbar pdv(, phi.alt) Y_l^m (theta,phi.alt) &= m hbar Y_l^m (theta, phi.alt)
$
so the $phi.alt$ dependence is of the form $tilde e^(i m phi.alt)$. Likewise
$
  bold(L)^2 ket(l","m) &= l(l+1)hbar^2 ket(l","m) \
  &=> 0 = [1/(sin theta) pdv(, theta) (sin theta pdv(, theta)) + 1/(sin^2 theta) pdv(, phi.alt, 2) + l(l+1)] Y_l^m
$
this is simply the equation defining $Y_l^m$.

We can find $Y_l^m$ by intelligent application of $L_-$ giving (see Sakurai)
$
  Y_l^m (theta,phi.alt) = (-1)^l/(2^l l!) sqrt((2l+1)/(4 pi) ((1+m)!)/((1-m)!)) e^(i m phi.alt) 1/(sin^m theta) d^(l-m)/(dd((cos theta)^(l-m))) (sin theta)^(2 l)
$
for $m >= 0$ and
$
  Y_l^(-m) (theta, phi.alt) = (-1)^m [Y_l^m (theta, phi.alt)]^*
$
We do not allow half-integer $l$ and therefore half-integer $m$. We would get
$
  e^(i m (2 pi)) = -1
$
after a $2 pi$ rotation if we did making the wavefunction doubly-valued which we don't allow.

== The radial equation
We assume a Hamiltonian of the form
$
  H = bold(p)^2/(2 m) + V (r)";  " r^2 = bold(x)^2
$
Then
$
  [bold(L),H]=[bold(L)^2,H] = 0
$
due to symmetry.

Then for the eigenstates $ket(E l m)$ we have
$
          H ket(E l m) & = E ket(E l m) \
  bold(L)^2 ket(E l m) & = l(l+1) hbar^2 ket(E l m) \
        L_z ket(E l m) & = m hbar ket(E l m)
$
Acting with $bra(bold(x)') dot (dots)$ on the first and using $braket(bold(x)', E l m) = R_(E l) Y_l^m$ we find the radial equation
$
  [- hbar^2/(2 m r^2) dv(, r) (r^2 dv(, r)) + (l(l+1)hbar^2)/(2m r^2) + V(r) ] R_(E l) (r) = E R_(E l) (r)
$
this is nice! The full eigenfunction is still $R_(E l) Y_l^m$.

We introduce $u_(E l) (r) equiv r R_(E l) (r)$ giving
$
  - hbar^2/(2 m) dv(u_(E l), r, 2) + [(l(l+1)hbar^2)/(2 m r^2) + V (r)] u_(E l) (r) = E u_(E l) (r)
$
with
$
  1 = integral dd(r) r^2 R_(E l)^* (r) R_(E l) (r) = integral dd(r) u_(E l)^* (r) u_(E l) (r)
$
This looks exactly like the time-independent Schrödinger wave equation but with the effective potential
$
  V_"eff" (r) = V(r) + (l(l+1)hbar^2)/(2 m r^2)
$

== Addition of angular momentum
The baseket for a spin $1/2$ particle is in the direct-product space spanned by ${ket(bold(x)')}$ and ${ket(+),ket(-)}$ meaning
$
  ket(bold(x)'\,plus.minus) = ket(bold(x)') times.circle ket(plus.minus) in cal(H)_"orb" times.circle cal(H)_"spin"
$
with anything in ${ket(bold(x)')}$ commuting with ${ket(+),ket(-)}$.

The rotation operator is the same but $bold(J) = bold(L) + bold(S) = bold(L) times.circle 1 + 1 times.circle bold(S)$ so
$
  D (R) = D^(("orb")) (R) times.circle D^(("spin")) (R) = exp((-i bold(L) dot hat(bold(n))phi.alt)/hbar) times.circle exp((-i bold(S) dot hat(bold(n)) phi.alt)/hbar)
$
since they commute. By $(A times.circle B) (x times.circle y) = (A x) times.circle (B y)$ it acts like
$
  D (R) (ket(psi)_"orb" times.circle ket(chi)_"spin") = (D^(("orb")) (R) ket(psi)_"orb") times.circle (D^(("spin")) (R) ket(chi)_"spin")
$
The wave function is just
$
  braket(bold(x)'","plus.minus, alpha) = psi_plus.minus (bold(x)')
$

Consider a system of two spin $1/2$ particles. Then the total spin operator is $ bold(S) = bold(S)_1 + bold(S)_2 = bold(S)_1 times.circle 1 + 1 times.circle bold(S)_2 $
again the different spin operators commute but within each space we have
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
We can expand a spin-state in terms of either the eigenkets of $bold(S)^2$ and $S_z$ or of $S_(1 z)$ and $S_(2 z)$.

1. For the ${m_1,m_2}$ representation we have eigenkets ${ket(+ +), ket(+ -), ket(- +), ket(- -)}$ with $m_i = plus.minus 1/2$.

2. For the ${s,m}$ representation, we have eigenkets,
$
  {ket(s=1","m=plus.minus 1","0), ket(s=0","m=0)}
$
with $s=1$ being the triplet and $s=0$ being the singlet.

The basekets are related by
$
   ket(s=1","m=1) & = ket(+ +) \
   ket(s=1","m=0) & = 1/sqrt(2) (ket(+ -) + ket(- +)) \
  ket(s=1","m=-1) & = ket(- -) \
   ket(s=0","m=0) & = 1/sqrt(2) (ket(+-) - ket(- +))
$
These are found by applying $S_-$ intelligently
$
  S_- & equiv S_(1 -) + S_(2 -) \
      & = (S_(1 x) - i S_(1 y)) + (S_(2 x) - i S_(2 y))
$
as an example
$
  S_- ket(s=1"," m=1) &= (S_(1 -) +S_(2 -)) ket(+ +) \
  sqrt((1+1)(1-1+1)) ket(s=1","m=0) &= sqrt((1/2+1/2)(1/2-1/2+1)) (ket(- +) + ket(+ -)) \
$
immediately simplifying to the second.

The coefficients of the RHS are called Clebsch-Gordan coefficients and are the elements of the transformation matrix for ${m_1,m_2} -> {s,m}$.

== Clebsch-Gordan coefficients
Consider two angular-momentum operators $bold(J)_1$ and $bold(J)_2$ satisfying the usual commutation relations and $[J_(1 k), J_(2 l)] = 0$. Then
$
  (1 - (i bold(J)_1 dot hat(bold(n)) dd(phi.alt, d: delta))/hbar) times.circle ((1 - (i bold(J)_2 dot hat(bold(n)) dd(phi.alt, d: delta))/hbar)) = 1 - (i (bold(J)_1 times.circle 1 + 1 times.circle bold(J)_2) dot hat(bold(n)) dd(phi.alt, d: delta))/hbar
$
we define the total angular momentum
$
  bold(J) equiv bold(J)_1 times.circle 1 + 1 times.circle bold(J)_2 = bold(J)_1 + bold(J)_2
$
So
$
  D_1 (R) times.circle D_2 (R) = exp((-i bold(J)_1 dot hat(bold(n))phi.alt)/hbar) times.circle exp((-i bold(J)_2 dot hat(bold(n))phi.alt)/hbar)
$
The total angular-momentum satisfies $[J_i,J_j] = i hbar epsilon_(i j k) J_k$ so $bold(J)$ is also an angular-momentum operator and can be treated as the generator of the entire system. This gives us two options for the choice of basekets

1. Simultaneous eigenkets of $bold(J)_1^2, bold(J)_2^2, J_(1 z)$, and $J_(2 z)$ denoted by $ket(j_1 j_2\; m_1 m_2)$ with defining equations
$
  bold(J)_1^2 ket(dots) & = j_1(j_1+1) hbar^2 ket(dots) \
      J_(1 z) ket(dots) & = m_1 hbar ket(dots) \
  bold(J)_2^2 ket(dots) & = j_2(j_2+1) hbar^2 ket(dots) \
      J_(2 z) ket(dots) & = m_2 hbar ket(dots)
$

2. Simultaneous eigenkets of $bold(J)^2, bold(J)_1^2, bold(J)_2^2$, and $J_z$. These mutually commute since
$
  bold(J)^2 = bold(J)_1^2 + bold(J)_2^2 + 2 J_(1 z) J_(2 z) + J_(1 +) J_(2 -) + J_(1-) J_(2 +)
$
in this case we denote the basekets by $ket(j_1 j_2 \; j m)$ with
$
  bold(J)_1^2 ket(dots) & = j_1 (j_1+1) hbar^2 ket(dots) \
  bold(J)_2^2 ket(dots) & = j_2 (j_2+1) hbar^2 ket(dots) \
    bold(J)^2 ket(dots) & = j(j+1) hbar^2 ket(dots) \
          J_z ket(dots) & = m hbar ket(dots)
$
typically this is just written as $ket(j m)$.

We consider a unitary transformation connecting these by
$
  ket(j_1 j_2 \; j m) = sum_(m_1, m_2) ket(j_1 j_2 \; m_1 m_2) underbracket(braket(j_1 j_2 \; m_1 m_2, j_1 j_2 \; j m), "Clebsch-Gordan")
$
where we just insert an identity. The matrix elements of the RHS are the Clebsch-Gordan coefficients. These are very important.

For $m eq.not m_1 + m_2$ they vanish.
#proof[
  Consider
  $
    (J_z - J_(1 z) - J_(2 z)) ket(j_1 j_2";" j m) = 0
  $
  then acting with $bra(j_1 j_2 ";" m_1 m_2) dot (dots)$ we find
  $
    (m - m_1 - m_2) braket(j_1 j_2 ";" m_1 m_2, j_1 j_2 ";" j m) = 0
  $
  so unless the Clebsch-Gordan coefficients vanish we have $m = m_1 + m_2$.
]

They do not vanish unless $abs(j_1 - j_2) <= j <= j_1 + j_2$.

We take the Clebsch-Gordan coefficients to be real. So they are their own inverses since they form a unitary matrix. Similarly we obtain
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
for the recursion relations see Sakurai.

Let $D^((j_1)) (R)$ be the rotation operator in the space spanned by the eigenket with eigenvalue $j_1$. Then we can write the Clebsch-Gordan series
$
  D^((j_1))_(m_1 m'_1) (R) D^((j_2))_(m_2 m'_2) (R) &= sum_(j,m,m') braket(j_1 j_2 \; m_1 m_2, j_1 j_2\; j m) \
  &times braket(j_1 j_2\;m'_1 m'_2, j_1 j_2 \; j m') D_(m m')^((j)) (R)
$
with $j: abs(j_1 -j_2) -> j_1 + j_2$. This comes from $D^((j_1)) times.circle D^((j_2))$ which is reducible.

#proof[
  The LHS is equivalent to
  $
    braket(j_1 j_2 \; m_1 m_2, D (R), j_1 j_2 \; m'_1 m'_2) &= braket(j_1 m_1, D (R), j_1 m'_1) braket(j_2 m_2, cal(D) (R), j_2 m'_2) \
    &= D^((j_1))_(m_1 m'_1) (R) D^((j_2))_(m_2 m'_2) (R)
  $
  we can also write
  $
    braket(j_1 j_2 \; m_1 m_2, D (R), j_1 j_2 \; m'_1 m'_2) &= sum_(j m,j'm') braket(j_1 j_2 \; m_1 m_2, j_1 j_2 \; j m) braket(j_1 j_2 \; j m, D (R), j_1 j_2 \; j' m') \
    & times braket(j_1 j_2\;j' m', j_1 j_2 \; m'_1 m'_2) \
    &= sum_(j m, j' m') braket(j_1 j_2 \; m_1 m_2, j_1 j_2 \; j m) D^((j))_(m m') (R) delta_(j j') \ & times braket(j_1 j_2 \; m'_1 m'_2, j_1 j_2 \; j' m')
  $
  which is the RHS.
]

== Schwinger's model
We consider two simple harmonic oscillators and denote these by the plus type and the minus type. We have accompanying annihilation and creation operators $a_+$ and $a_+^dagger$ etc. Similarly we define
$
  N_+ equiv a_+^dagger a_+";  " N_- equiv a_-^dagger a_-
$
we assume the usual commutation relations hold for each
$
   [a_+,a_+^dagger] & = 1 \
          [N_+,a_+] & = - a_+ \
  [N_+, a_+^dagger] & = a_+^dagger
$
and similarly for $a_-$. We assume they commute with eachother $[a_+,a_-^dagger]=[a_-,a_+^dagger]=0$. We call the oscillators decoupled or independent. $N_-$ and $N_+$ commute so there exists a simultaneous eigenket with eigenvalues $n_+$ and $n_-$. We have
$
  N_+ ket(n_+\,n_-) = n_+ ket(n_+\,n_-)";  " N_- ket(n_+\,n_-) = n_- ket(n_+\,n_-)
$
likewise we have
$
  a_+^dagger ket(n_+\,n_-) = sqrt(n_++1) ket(n_+ +1\,n_-)";  " a_+ ket(n_+\,n_-) = sqrt(n_+) ket(n_+-1\,n_-)
$
and similarly for $a_-$ and $a_-^dagger$. From
$
  a_+ ket(0\,0) = 0";  " a_- ket(0\,0) = 0
$
we obtain by succesive application
$
  ket(n_+\,n_-) = ((a_+^dagger)^(n_+) (a_-^dagger)^(n_-))/(sqrt(n_+ !) sqrt(n_- !)) ket(0\,0)
$
We define
$
  J_+ equiv hbar a_+^dagger a_-";  " J_- equiv hbar a_-^dagger a_+
$
and
$
  J_z equiv hbar/2 (a_+^dagger a_+ - a_-^dagger a_-) = hbar/2 (N_+ - N_-)
$
These satisfy
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
We define $N equiv N_+ + N_-$ then
$
  bold(J)^2 equiv J_z^2 + 1/2 (J_+ J_- + J_- J_+) = (hbar^2/2) N (N/2 + 1)
$

We associate spin up with one quantum unit of the plus-type and spin down with one quantum unit of the minus-type. $n_+$ and $n_-$ are then the amount of spin ups and downs respectively. Then $J_+$ destroys one unit of spin down with $z$-component $-hbar/2$ and creates one unit of spin up with $z$-component $hbar/2$. So it increases by $hbar$ and a similar interpretation holds for $J_-$. $J_z$ just counts $hbar/2$ times the difference of $n_+$ and $n_-$.

We find
$
  J_+ ket(n_+\,n_-) & = sqrt(n_- (n_++1)) hbar ket(n_+ + 1\, n_- -1) \
  J_- ket(n_+\,n_-) & = sqrt(n_+ (n_- + 1)) hbar ket(n_+ -1\, n_- + 1) \
  J_z ket(n_+\,n_-) & = 1/2 (n_+-n_-) hbar ket(n_+\,n_-)
$
letting $n_+ -> j + m$ and $n_- -> j - m$ these reduce to the relations we are used to. And the eigenvalue of $bold(J)^2$ becomes $hbar^2 j(j+1)$. This is expected since we constructed $J_plus.minus$ and $J_z$ such that they satisfy the same commutation relations. We can then write
$
  j equiv (n_+ + n_-)/2";  " m equiv (n_+ - n_-)/2
$
and we can write the general state as
$
  ket(j\,m) = ((a_+^dagger)^(j+m) (a_-^dagger)^(j-m))/sqrt((j+m)!(j-m)!) ket(0)
$
Letting $j = m$ we find the special case
$
  ket(j\,j) = ((a_+^dagger)^(2 j))/sqrt((2j)!) ket(0)
$
which corresponds to the physical case where the eigenvalue of $J_z$ is as large as possible for a given $j$. We can treat this state like $2 j$ spin $1/2$ particles with all spins pointing up. Generally some object with high $j$ can be treated as being made up from $j+m$ spin up particles and $j-m$ spin down particles.

== Landau levels
The Hamiltonian in the presence of a magnetic field is given by
$
  H = 1/(2 m) (bold(p)- e/c bold(A))^2 "with" bold(B) = nabla times bold(A)
$
and the equation of motion is
$
  m dot.double(bold(x)) = e/c dot(bold(x)) times bold(B)
$
We let $bold(B) = B_0 hat(bold(z))$ and use the symmetric gauge
$
  bold(A) = - 1/2 bold(x) times bold(B) = B_0/2 vec(-y, x, 0)
$
Then
$
  H &= 1/(2 m) vec(p_x + (e B_0)/(2 c) y, p_y - (e B_0)/(2 c) x)^2 \
  &= underbrace(1/(2 m) (p_x^2 + p_y^2) + 1/2 m (omega_L/2)^2 (x^2+y^2), "2D SHO with" omega= omega_L\/2) - omega_L/2 underbrace((x p_y - y p_x), L_z)
$
We define
$
  x_i = sqrt(hbar/(m omega_L)) (a_i + a_i^dagger)
$
and similarly for $p_i$. Then
$
  a_plus.minus = 1/sqrt(2) (a_1 minus.plus i a_2)";  " [a_+,a_-] = 0 dots
$
and $[a_+,a_+^dagger] = 1 dots$ so these are uncoupled and
$
  L_z = hbar(a_+^dagger a_+ - a_-^dagger a_-)
$
Then
$
  H &= hbar (omega_L)/2 (N_+ + N_- + 1) - hbar omega_L/2 (a_+^dagger a_+ - a_-^dagger a_-) \
  &= hbar omega_L (a_-^dagger a_- + 1/2) -> E_n = hbar omega_L (n_- + 1/2)
$
and $[H,L_z] = 0$ with eigenket $ket(n_+\,n_-)$ so
$
  L_z = hbar (n_+-n_-)
$
The found energies are the Landau levels! These have infinite degeneracy since to each $n_-$ we have infinite $n_+$ leading to different $L_z$.

== Bell's inequality
We consider a two-electron system in a spin-singlet state so with total spin zero. We write this state as
$
  ket("singlet") = 1/sqrt(2) (ket(hat(bold(z))+\;hat(bold(z))-)- ket(hat(bold(z))-\;hat(bold(z))+))
$
with the obvious interpretation. The two electrons are correlated in the sense that if we measure the spin of one electron and collapse to $ket(hat(bold(z))+\;hat(bold(z))-)$ then the second electron is forced to take some value. This correlation happens even at distance. A relevant example is the decay $eta -> mu^+ + mu^-$ where $eta$ has spin zero and each $mu^(plus.minus)$ has spin half. We imagine measuring the spin of $mu^+$ at some detector A and measureing the spin of $mu^-$ at some detector B. These measurements would be correlated in the obvious manner even at macroscopic distance.

Recall for this type of system
$
  ket(hat(bold(x)) plus.minus) = 1/sqrt(2) (ket(hat(bold(z))+) plus.minus ket(hat(bold(z))-))
$
and similarly with $hat(bold(z))<--> hat(bold(x))$. Then one can show by substitution that
$
  ket("singlet") = 1/sqrt(2) (ket(hat(bold(x))-";"hat(bold(x))+)-ket(hat(bold(x))+";"hat(bold(x))-))
$
which is the same apart form an overall phase. We consider measuring $S_x$ and/or $S_z$ at A and B.

If we measure $S_z$ at A and $S_x$ at B then the measurements are randomly correlated. If we measure the same at both A and B then we find total anti-correlation. Lastly if we do no measurement at A then the measurement at B is random. It seems like the second particle knows which spin component of the first particle is being measured instantaneously. This seems to contradict Einstein's locality principle and the problem itself is known as the Einstein-Podolsky-Rosen paradox. Some ascribe the aforementioned measurement behaviour to be a result of hidden variables. We would like to know if these hidden variable theories make predictions that can distinguish them from quantum mechanics. This is what Bell's inequality does.

We derive the Bell's inequality using a framework by Wigner. We consider an ensemble of spin $1/2$ particles with zero total angular momentum. We define different species of particles. Some have $(hat(bold(z))+,hat(bold(x))-)$ while other have $(hat(bold(z))-,hat(bold(x))-)$ etc. To ensure the angular momentum vanishes any particular pair must match so if particle 1 has $(hat(bold(z))+,hat(bold(x))-)$ then particle 2 must have $(hat(bold(z))-,hat(bold(x))+)$. For two particles we have four possible pairs with equal probability $1/4$. This leads to the same measurement correlations as the previous.

Consider the arbitrary unit vectors ${hat(bold(a)), hat(bold(b)), hat(bold(c))}$ then we have e.g. the type $[hat(bold(a))-,hat(bold(b))+,hat(bold(c))+]$ etc. and again we must necessarily have matching pairs. So we  have 8 types with each type having a population $N_i$. Suppose A measures $bold(S)_1 dot hat(bold(a))$ to be $+$ and B measures $bold(S)_2 dot hat(bold(b))$ to be $+$ then this pair belongs to either $N_3$ or $N_4$. We have
$
  N_3 + N_4 <= (N_2 +N_4) + (N_3 + N_7)
$
and
$
  P(hat(bold(a))+";"hat(bold(b))+) = (N_3 +N_4)/(sum N_i)
$
with the obvious meaning. Similarly
$
  P(hat(bold(a))+";"hat(bold(c))+) & = (N_2 + N_4)/(sum N_i) \
  P(hat(bold(c))+";"hat(bold(b))+) & = (N_3+N_7)/(sum N_i)
$
immediately giving Bell's inequality
$
  P(hat(bold(a))+";"hat(bold(b))+) <= P(hat(bold(a))+";"hat(bold(c))+) + P(hat(bold(c))+";"hat(bold(b))+)
$
for hidden variable theories this inequality is satisfied.

We now show that this inequality is violated in quantum mechanics. We characterize all singlet states by $ket("singlet")$ defined before. Then we can calculate all terms in Bell's inequality. We define $bold(hat(a)) dot bold(hat(b)) = cos theta_"ab"$ and take both to be in the $x y$-plane with $bold(hat(a))$ being $bold(hat(x))$ without loss of generality. Then
$
  ket(S_a";"plus.minus) equiv ket(S_x";"plus.minus) = 1/sqrt(2) (ketbra(+)-ketbra(-))
$
we rotate this by $theta_"ab"$ about the $z$-axis
$
  ket(S_b";"+) &= exp((-i S_z theta_"ab")/hbar) ket(S_x";"+) \
  &= exp((-i theta_"ab")/2) 1/sqrt(2) ket(+) + exp((i theta_"ab")/2) 1/sqrt(2) ket(-)
$
Then
$
  braket(S_a";"+, S_b";"+) & = cos theta_"ab"/2
$
so
$
  P(S_a +";" S_b +) = 1/2 cos^2 theta_"ab"/2
$
with the factor $1/2$ coming from $S_a$ being $+$.

We find
$
  sin^2 theta_"ab"/2 <= sin^2 theta_"ac"/2 + sin^2 theta_"cb"/2
$
this should hold for all angles. We try $theta_"ab" = 2 theta$ and $theta_"ac" = theta_"bc" = theta$ then it is violated for $0 < theta < pi\/2$!
