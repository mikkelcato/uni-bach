//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *solved problems qp*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    problems from _Modern Quantum Mechanics_ by Sakurai and various problem sets.
  ],
)

#set heading(numbering: "1.1.1")
// for the desired level
#show selector(heading.where(level: 2)): set heading(numbering: none)
// update the counter because it reset to 0
#counter(heading).step(level: 1)

= Fundamentals
== 1.2
Prove
$
  [A B, C D] = - A C {D, B} + A {C,B} D - C {D,A} B + {C,A} D B
$
we just expand the left-hand side
$
  [A B, C D] & = A B C D - C D A B \
             & = A B C D - C D A B + A C B D - A C B D \
             & = A {B, C} D - C D A B - A C B D \
             & = A {B, C} D - C D A B - A C B D - C A D B + C A D B \
             & = A {B, C} D - C {A, D} B - A C B D + C A D B \
             & = A {B, C} D - C {A, D} B - A C B D + C A D B - A C D B + A C D B \
             & = A {B, C} D - C {A, D} B - A C {B, D} + {A, C} D B
$


== 1.6
Using the rules of bra-ket algebra, prove or evaluate the following:

$ tr(X Y) = tr(Y X) $

We have by definition
$
  tr (X Y) & = sum_a' braket(a', X Y, a') \
           & = sum_a' sum_b' braket(a', X, b') braket(b', Y, a') \
           & = sum_a' sum_b' braket(b', Y, a') braket(a', X, b') \
           & = sum_b' braket(b', Y X, b') equiv tr (Y X)
$

$ (X Y)^dagger = Y^dagger X^dagger $

We define
$
  X Y ket(alpha) = X ket(beta)
$
this has a dual $bra(beta) X^dagger$, likewise $Y ket(alpha)$ has a dual $bra(alpha) Y^dagger$ so subbing in we get
$
  X Y ket(alpha) <-> bra(alpha) Y^dagger X^dagger
$
meaning that $(X Y)^dagger = Y^dagger X^dagger$.

$exp(i f(A)) = ?$ in ket-bra form, where $A$ is a Hermitian operator whose eigenvalues are known.

We have
$
  e^(i f(A)) & = e^(i f(A)) sum_a' ketbra(a') \
             & = sum_a' e^(i f(A)) ketbra(a') \
             & =_"herm." sum_a' e^(i f(a')) ketbra(a')
$

$sum_a' psi_a'^* (bold(x)') psi_a' (bold(x)'')$ with $psi_a' (bold(x)') = braket(bold(x)', a')$.

We have
$
  sum_a' psi_a'^* (bold(x)') psi_a' (bold(x)'') &= sum_a' braket(a', bold(x)') braket(bold(x)'', a') \
  &= sum_a' braket(bold(x)'', a') braket(a', bold(x)') \
  &= braket(bold(x)'', bold(x)') = delta^3 (bold(x)''-bold(x)')
$

== 1.7
Consider two kets $ket(alpha)$ and $ket(beta)$. Suppose $braket(a', alpha), braket(a'', alpha),dots$ and $braket(a', beta), braket(a'', beta),dots$ are all known, where $ket(a'),dots$ form a complete set of base kets. Find the matrix representation of the operator $ketbra(alpha, beta)$ in this basis.

Let $X = ketbra(alpha, beta)$, we want the matrix elements of this operator,
$
  braket(a^(i), X, a^(j)) = braket(a^(i), alpha) braket(beta, a^(j))
$
these are both known and we are done.

We now consider a spin-$1\/2$ system and let $ket(alpha)$ and $ket(beta)$ be $ket(+)$ and $ket(S_x";"+)$ respectively. Write down the square matrix corresponding to $ketbra(alpha, beta)$ in the usual basis.

We just find the matrix elements, with $ ket(S_x";"plus.minus) = 1/sqrt(2) ket(+) plus.minus 1/sqrt(2) ket(-) $ e.g.
$
  X_(0,0) = braket(+, +) braket(S_x";"+, +) = 1/sqrt(2)
$
the rest are $X_(0,1) = 1/sqrt(2)$, $X_(1,0) = X_(1,1) = 0$. So we have
$
  ketbra(+, S_x";"+) = 1/sqrt(2) mat(1, 1; 0, 0)
$

== 1.8
Suppose $ket(i)$ and $ket(j)$ are eigenkets of some Hermitian operator $A$. Under what condition can we conclude that $ket(i)+ket(j)$ is also an eigenket of $A$?

By definition $A ket(i) = i ket(i)$ and $A ket(j) = j ket(j)$ so
$
  A (ket(i)+ket(j)) = i ket(i) + j ket(j) =^! a (ket(i)+ket(j))
$
for this to be true we require $a = i = j$.

== 1.10
Trivial

== 1.11
Construct $ket(bold(S) dot hat(bold(n))";"+)$ so that
$
  bold(S) dot hat(bold(n)) ket(bold(S)dot hat(bold(n))";"+) = hbar/2 ket(bold(S) dot hat(bold(n))";"+)
$
from the figure
$
  hat(bold(n)) = cos alpha sin beta bold(e)_x + sin alpha sin beta bold(e)_y + cos beta bold(e)_z
$

The operator then becomes
$
  bold(S) dot hat(bold(n)) &= cos alpha sin beta S_x + sin alpha sin beta S_y + cos beta S_z \
  &= hbar/2 mat(cos beta, e^(-i alpha) sin beta; e^(i alpha) sin beta, - cos beta)
$
the eigenvector is given by
$
  a cos beta + b e^(- i alpha) sin beta & = a \
    a e^(i alpha) sin beta - b cos beta & = b
$
which can be written as
$
  a (1-cos beta) & = b e^(- i alpha) sin beta \
      a sin beta & = b e^(-i alpha) (1 + cos beta)
$
the first one
$
  a & = -e^(-i alpha) (sin beta)/(cos beta - 1) b \
    & = e^(-i alpha) cot(beta/2) b
$
and
$
             a (1-cos beta) & = b e^(- i alpha) sin beta \
  abs(a)^2 (1 - cos beta)^2 & = (1 - abs(a)^2) sin^2 beta \
    4 abs(a)^2 sin^4 beta/2 & = 4 (1 - abs(a)^2) sin^2 beta/2 cos^2 beta/2 \
      abs(a)^2 sin^2 beta/2 & = (1-abs(a)^2)cos^2 beta/2 \
                   abs(a)^2 & = cos^2 beta/2 => a = cos beta/2
$
then
$
  b =e^(i alpha) sin beta/2
$
so we end up with
$
  ket(bold(S) dot hat(bold(n))";"+) = cos beta/2 ket(+) + e^(i alpha) sin beta/2 ket(-)
$
by $beta -> beta + pi$ we get the spin-down state
$
  ket(bold(S) dot hat(bold(n))";"-) = - sin beta/2 ket(+) + e^(i alpha) cos beta/2 ket(-)
$

== 1.12
The Hamiltonian for a two-state system is given by
$
  H = a (ketbra(1, 1)-ketbra(2, 2)+ketbra(1, 2)+ketbra(2, 1))
$
find the eigenvalues and correspond eigenkets.

We identify $ket(1)=ket(+)$ and $ket(2)=ket(-)$. Then
$
  H = (2 a)/hbar (S_x + 0 + S_z) = (2 sqrt(2) a)/hbar bold(S) dot hat(bold(n))
$
with
$
  hat(bold(n)) = 1/sqrt(2) vec(1, 0, 1) => alpha = 0", " beta = pi/4
$
so we immediately get
$
  E_(plus.minus) = plus.minus (2 sqrt(2) a)/hbar hbar/2 = plus.minus sqrt(2) a
$
with eigenkets
$
  ket(E_+) & = cos pi/8 ket(+) + sin pi/8 ket(-) \
  ket(E_-) & = - sin pi/8 ket(+) + cos pi/8 ket(-)
$

== 1.13
A two-state system is characterized by
$
  H = H_(11) ketbra(1, 1) + H_(22) ketbra(2, 2) + H_(12) (ketbra(1, 2) + ketbra(2, 1))
$

We can rewrite
$
  H &= (H_(11)+H_(22))/2 {ketbra(1, 1)+ketbra(2, 2)} + (H_(11)-H_(22))/2 {ketbra(1, 1)-ketbra(2, 2)} + H_(12) {ketbra(1, 2)+ketbra(2, 1)} \
  &= A II + (2B)/hbar S_z + (2C)/hbar S_x = A II + 2/hbar (C S_x + 0 + B S_z) \
  &= A II + 2/hbar sqrt(B^2+C^2) bold(S) dot hat(bold(n))
$
where
$
  hat(bold(n)) = 1/sqrt(B^2+C^2) vec(C, 0, B)
$
so $alpha = 0$ and $tan beta = C\/B = 2 H_(12) \/(H_(11)-H_(22))$. The eigenvalue of the first term is just $A$ and the eigenvalue of the second term is $plus.minus sqrt(B^2+C^2)$ giving
$
  E_(plus.minus) = A plus.minus sqrt(B^2+C^2)
$

== 1.17
Let $A$ and $B$ be observables. Suppose the simultaneous eigenkets of $A$ and $B$ ${ket(a","b)}$ form a complete orthonormal set of base kets. Can we always conclude that $ [A,B] = 0 $

By definition $A ket(a","b) = a ket(a","b)$ and $B ket(a","b) = b ket(a","b)$ then
$
  A B & = A B sum_(a b) ketbra(a","b) \
      & = sum_(a b) A B ketbra(a","b) \
      & = sum_(a b) a b ketbra(a","b) \
      & = sum_(a b) b a ketbra(a","b) \
      & = sum_(a b) b A ketbra(a","b) \
      & = sum_(a b) B A ketbra(a","b) \
      & = B A => [A, B] = 0
$

== 1.18
Two Hermitian operators anticommute
$
  {A, B} = A B + B A = 0
$
is it possible to have a simultaneous eigenket of $A$ and $B$.

By assumption $A B = - B A$, if we had a simultaneous eigenket $ket(a","b)$ then
$
  A B ket(a","b) & = A b ket(a","b) \
                 & = a b ket(a","b) \
                 & =^"and" - B A ket(a","b) \
                 & = - a b ket(a","b)
$
so $a b = - a b$, meaning either $a = 0, b = 0$ or $a = b = 0$.

== 1.23
Particle in a box---evaluate $x"-"p$ uncertainty product.

Eigenfunctions are
$
  braket(x, n) = u_n (x) = sqrt(2/a) sin (n pi x)/a
$
so
$
  expval(x) = braket(n, x, n) &= integral_0^a braket(n, x, x')braket(x', n) dd(x)' \
  &= integral_0^a braket(n, x')x' braket(x', n) dd(x)' \
  &= 2/a integral_0^a x' sin^2 (n pi x')/a dd(x)' = a/2
$
and
$
  expval(x^2) &= 2/a integral_0^a x'^2 sin^2 (n pi x')/a dd(x)' = a^2/6 (2 - 3/(n^2 pi^2))
$
so
$
  expval((Delta x)^2) = a^2/12 (1 - 6/(n^2 pi^2))
$
and
$
  expval(p) = braket(n, p, n) &= integral_0^a braket(n, x') braket(x', p, n) dd(x)' \
  &= integral_0^a braket(n, x') (- i hbar) pdv(, x') braket(x', n) dd(x)' \
  &= 0
$
and
$
  expval(p^2) = braket(n, p^2, n) &= integral_0^a braket(n, x') braket(x', p^2, n) dd(x)' \
  &= (-i hbar)^2 integral_0^a braket(n, x') pdv(, x', 2) braket(x', n) dd(x)' \
  &= - (i hbar)^2 ((n pi)/a)^2 2/a integral_0^a sin^2 (n pi x')/a dd(x)' \
  &= (hbar^2 n^2 pi^2)/(a^2) = expval((Delta p)^2)
$
so
$
  expval((Delta x)^2) expval((Delta p)^2) &= (hbar^2 n^2 pi^2)/a^2 a^2/12 (1 - 6/(n^2 pi^2)) \
  &= hbar^2/4 ((n^2 pi^2)/3 - 2)
$

== 1.28
Construct the transformation matrix that connect the $S_z$ diagonal basis  to the $S_x$ diagonal basis. Show that your result is consisten with the general relation:
$
  U = sum_r ketbra(b^((r)), a^((r)))
$

We have
$
  ket(+) = vec(1, 0)",  " ket(-) = vec(0, 1)
$
$
  ket(S_x";"+) = 1/sqrt(2) vec(1, 1)",  " ket(S_x";"-) = 1/sqrt(2) vec(1, -1)
$
to find $U$ we find elements $braket(a^((i)), b^((j)))$:
$
  U &= mat(braket(+, S_x";"+), braket(+, S_x";"-); braket(-, S_x";"+), braket(-, S_x";"-)) \
  &= 1/sqrt(2) mat(1, 1; 1, -1)
$
and we need to check
$
  ketbra(+, S_x";"+) + ketbra(-, S_x";"-) &= 1/sqrt(2) vec(1, 0) vecrow(1, 1)+1/sqrt(2) vec(0, 1) vecrow(1, -1) \
  &= 1/sqrt(2) mat(1, 1; 1, -1)
$

== 1.29
Suppose $f(A)$ is a function of a Hermitian operator $A$ with the property $A ket(a') = a' ket(a')$. Evaluate $braket(b'', f(A), b')$ when the transformation matrix from the $a'$ basis to the $b'$ basis is known.

By definition
$
  U = sum_l ketbra(b^((l)), a^((l)))",  " U ket(a^((l))) = ket(b^((l)))
$
with matrix elements $braket(a^((l)), U, a^((k))) = braket(a^((l)), b^((k)))$, which we know. Then calculate
$
  braket(b'', f(A), b') & = sum_a' braket(b'', f(A), a')braket(a', b') \
                        & = sum_a' f(a') braket(b'', a')braket(a', b')
$
both of the brakets are just matrix elements so they are known, and we are done.

Using the continuum analogue of this result, evaluate
$
  braket(bold(p)'', F(r), bold(p)')
$
where $r^2 = x^2 + y^2 + z^2$.

We have
$
  braket(bold(p)'', F(r), bold(p)') &= integral dd(x', 3) braket(bold(p)'', F(r), bold(x)') braket(bold(x)', bold(p)') \
  &= integral dd(x', 3) F(bold(x)') braket(bold(p)'', bold(x)')braket(bold(x)', bold(p)') \
  &= 1/(2 pi hbar)^3 integral dd(x', 3) F(r') e^((i(bold(p)'-bold(p)'') dot bold(x)')\/hbar)
$
defining $bold(q) = bold(p)'-bold(p)''$ along $hat(bold(z))$ we get
$
  bold(q) dot bold(x)' = q r' cos theta
$
then the angular part of the integral can be done
$
  braket(bold(p)'', F(r), bold(p)') &= 1/(2 pi hbar)^3 integral_0^(2 pi) dd(phi) integral_0^pi dd(theta) integral_0^oo dd(r') r'^2 sin theta F(r') e^((i q r' cos theta)\/hbar) \
  &= 1/(2 pi^2 hbar^2 q) integral_0^oo dd(r') F(r') r' sin((q r')/hbar)
$
where the integral can be solved by $u = cos theta$.

== 1.31
Verify
$
  [x_i, G(bold(p))] = i hbar pdv(G, p_i)",  " [p_i,F(bold(x))] = - i hbar pdv(F, x_i)
$
if functions $F$ and $G$ can be expressed as power series in their arguments---and knowing the canonical commutation relations.

We start with the first. By assumption we can write
$
  G(bold(p)) = sum_(m,n,l) c_(m,n,l) p_x^m p_y^n p_z^l
$
we know
$
  [x_i,p_j] =^(i eq.not j) 0",  " [x_i,p_i] = x_i p_i - p_i x_i = i hbar
$
so we just need to calculate
$
  [x_i,p_i^n]
$
without loss of generality we pick $x_i = x$ and $p_i = p_x = p$, then
$
  [x, p^n] & = x p^n - p^n x \
           & = (x p) p^(n-1) - p^n x \
           & = (i hbar + p x) p^(n-1) - p^n x \
           & = i hbar p^(n-1) + p (x p) p^(n-2) - p^n x \
           & = i hbar p^(n-1) + p (i hbar + p x) p^(n-2) - p^n x \
           & = 2 i hbar p^(n-1) + p^2 x p^(n-2) - p^n x \
           & = n i hbar p^(n-1) + p^n x - p^n x \
           & = n i hbar p^(n-1) = i hbar pdv(p^(n), p)
$
so
$
  [x_i,p_i^n] = i hbar pdv(p_i^n, p_i)
$
since $[x_i, p_j] = 0 = [p_i,p_j]$ the result immediately follows. Similarly for the second term here we need to calculate
$
  [p,x^n] & = p x^n - x^n p \
          & = (p x) x^(n-1) - x^n p \
          & = (x p - i hbar) x^(n-1) - x^n p \
          & = - i hbar x^(n-1) + x (p x) x^(n-2) -x^n p \
          & = - n i hbar x^(n-1) + x^n p - x^n p \
          & = - n i hbar x^(n-1) = - i hbar pdv(x^n, x)
$
so
$
  [p_i,x_i^n] = - i hbar pdv(x_i^n, x)
$
and the result follows since $[x_i,x_j] = 0$.

Evaluate $[x^2, p^2]$.

$
  [x^2,p^2] & = x^2 p^2 - p^2 x^2 \
            & = x x p^2 - p^2 x x - x p^2 x + x p^2 x \
            & = x [x, p^2] + [x,p^2] x \
            & = x (2 i hbar p) + (2 i hbar p) x \
            & = 2 i hbar {x,p}
$

#pagebreak()
= Quantum Dynamics
== eq. 2.15
Show
$
  cal(U) (t_0 + dd(t), t_0) = 1 - i Omega dd(t)
$
is unitary and satisfies the composition property.

Take
$
  cal(U) cal(U)^dagger & = (1 - i Omega dd(t))(1 + i Omega^dagger dd(t)) \
                       & = (1 - i Omega dd(t)) (1 + i Omega dd(t)) \
                       & = 1 + Omega^2 dd(t)^2 tilde.eq 1
$
and
$
  cal(U) (t_0 + dd(t)_1,t_0) cal(U) (t_0 + dd(t)_2, t_0) &= (1 - i Omega dd(t)_1) (1 - i Omega dd(t)_2) \
  &tilde.eq 1 -i Omega dd(t)_1 - i Omega dd(t)_2 \
  &= 1 - i Omega (dd(t)_1 + dd(t)_2) \
  &= cal(U) (t_0+(dd(t)_1+dd(t)_2),t_0)
$

== eq. 2.82
Show
$
  expval(bold(x)) -> expval(bold(x)) + expval(dd(bold(x))')
$
in Heisenberg and Schrödinger picture.

In Heisenberg $bold(x)$ becomes
$
  bold(x) &arrow (1 + (i bold(p) dot dd(bold(x))')/hbar) bold(x) (1 - (i bold(p) dot dd(bold(x)'))/hbar) \
  &= bold(x) + (i/hbar) [bold(p) dot dd(bold(x)'),bold(x)] \
  &= bold(x) + dd(bold(x)')
$
so
$
  expval(bold(x)) = expval(bold(x), a) -> expval(cal(J)^dagger bold(x) cal(J), a) = expval(bold(x)+dd(bold(x)'), a) = expval(bold(x)) + expval(dd(bold(x)'))
$

In Schrödinger picture
$
  cal(J) (dd(bold(x)')) ket(a) = (1 - (i bold(p) dot dd(bold(x)'))/hbar) ket(a)
$
so
$
  braket(alpha (1 + (i bold(p) dot dd(bold(x)'))/hbar), bold(x), (1 - (i bold(p) dot dd(bold(x)'))/hbar) alpha) = expval(bold(x)) + expval(dd(bold(x)'))
$

== eq. 2.369

== 2.1
Spin precession can also be solved in the Heisenberg picture. Consider
$
  H = - ((e B)/(m c)) S_z = omega S_z
$
write the Heisenberg equations of motion for the time-dependent operators $S_z (t)$, $S_y (t)$ and $S_z (t)$, and solve them.

The time-evolution operator becomes
$
  cal(U) = exp((- i S_z omega t)/hbar)
$
note that
$
  H = cal(U)^dagger cal(U) H = cal(U)^dagger H cal(U) = omega cal(U)^dagger S_z cal(U) = omega S_z^((H))
$
so
$
  dv(S_x^((H)), t) &= 1/(i hbar) [S_x^((H)),H] = omega/(i hbar) [S_x^((H)),S_z^((H))] = -omega S_y^((H)) \
  dv(S_y^((H)), t) &= omega S_x^((H)) \
  dv(S_z^((H)), t) &= 0
$
the first and second gives
$
  dv(S_x^((H)), t, 2) = - omega^2 S_x^((H))" and " dv(S_y^((H)), t, 2) = - omega^2 S_y^((H))
$
this is basic harmonic motion
$
  S_x^((H)) (t) & = - S_y^((H)) (0) sin omega t + S_x^((H)) (0) cos omega t \
  S_y^((H)) (t) & = S_x^((H)) (0) sin omega t + S_y^((H)) (0) cos omega t
$

== 2.2
Suppose that
$
  H = H_(11) ketbra(1) + H_(22) ketbra(2) + H_(12) ketbra(1, 2)
$
what principle is violated?

We have $H^dagger eq.not H$, so it is not Hermitian due to the cross-term.

We assume $H_(11) = H_(22) = 0$, so
$
  H = H_(12) ketbra(1, 2) => H^2 = H_(12) ketbra(1, 2) ketbra(1, 2) = 0 => H^n = 0", " n >= 2
$
this makes the time-evolution operator
$
  cal(U) (t) = 1 - i/hbar H t = 1 - i/hbar t H_(12) ketbra(1, 2)
$
consider $ket(alpha) = ket(2)$ at $t = 0$, then for $t > 0$,
$
  ket(alpha";"t) = cal(U) (t) ket(2) = ket(2) - i/hbar t H_(12) ket(1)
$
so
$
  braket(alpha";"t) = 1 + H_(12)^2/hbar^2 t^2
$
which is nonsense.

== 2.3
We let $alpha = 0$ then
$
  hat(bold(n)) = sin beta bold(e)_x + cos beta bold(e)_z
$
and
$
  ket(alpha) = ket(bold(S) dot hat(bold(n))","+) = cos beta/2 ket(+) + sin beta/2 ket(-)
$
and
$
  H = omega S_z
$
at time $t$
$
  ket(alpha";"t) & = cal(U) (t) ket(alpha) \
  &= exp[(- i omega t)/2] cos beta/2 ket(+) + exp[(i omega t)/2] sin beta/2 ket(-)
$
the probability for
$
  ket(S_x ";" t) = 1/sqrt(2) ket(+) + 1/sqrt(2) ket(-)
$
is
$
  P & = abs(braket(S_x";"t, alpha";"t))^2 \
    & = 1/2 (1 + sin beta cos omega t)
$
the expectation value
$
  expval(S_x, alpha";"t) &= braket(hbar/2 ketbra(+, -) + hbar/2 ketbra(-, +), alpha";"t) \
  &= hbar/2 sin beta cos omega t
$
for $beta -> 0 => expval(S_x) = 0$ and $beta -> pi\/2 => expval(S_x) = (hbar\/2) cos omega t => "precession in "x y"-plane"$.

== 2.10
Let $ket(a')$ and $ket(a'')$ be eigenstates of a Hermitian operator $A$ with eigenvalues $a'$ and $a''$ ($a' eq.not a''$). The Hamiltonian is given by
$
  H = ket(a') delta bra(a'') + ket(a'') delta bra(a')
$

a. Write down the eigenstates and eigenvalues of $H$.

We have
$
  H = mat(0, delta; delta, 0)
$
so
$
  det(H) = lambda^2 - delta^2 = 0 => lambda = plus.minus delta
$
for the eigenstates
$
  mat(0, delta; delta, 0) vec(a, b) = delta vec(a, b)
$
giving
$
  delta b = delta a => a = b
$
so
$
  ket(delta) = 1/sqrt(2) vec(1, 1) = 1/sqrt(2) ket(a') + 1/sqrt(2) ket(a'')
$
and
$
  mat(0, delta; delta, 0) vec(a, b) = -delta vec(a, b)
$
giving
$
  delta b =- delta a => a = - b
$
so
$
  ket(-delta) = 1/sqrt(2) ket(a') - 1/sqrt(2) ket(a'')
$

b. Suppose the system is in $ket(a')$ at $t = 0$. Write the state vector in the Schrödinger picture.

$H$ is time-independent so
$
  cal(U) (t,0) = exp((- i H t)/hbar)
$
using the previous we can write
$
   ket(a') & = 1/sqrt(2) ket(delta) + 1/sqrt(2) ket(-delta) \
  ket(a'') & = 1/sqrt(2) ket(delta) - 1/sqrt(2) ket(-delta)
$
so
$
  ket(alpha";"t) &= exp((- i H t)/hbar) ket(a') \
  &= 1/sqrt(2) exp((-i H t)/hbar) (ket(delta) + ket(-delta)) \
  &= 1/sqrt(2) exp((- i delta t)/hbar) ket(delta) + 1/sqrt(2) exp((i delta t)/hbar) ket(-delta)
$

c. Probability to find the system in $ket(a'')$ at some $t$.

$
  abs(braket(a'', alpha";"t))^2 &= 1/4 abs((bra(delta) - bra(-delta))(exp((-i delta t)/hbar) ket(delta) + exp((i delta t)/hbar) ket(-delta)))^2 \
  &= 1/4 abs(exp((-i delta t)/hbar) - exp((i delta t)/hbar))^2 \
  &= 1/4 abs(- 2 i sin (delta t)/hbar)^2 = sin^2 (delta t)/hbar
$
== 2.12
A one-dimensional simple harmonic oscillator with frequency $omega$ is in an initial state
$
  ket(alpha) = 1/sqrt(2) ket(0) + (e^(i delta))/sqrt(2) ket(1)
$
with $delta$ real.

Note that $H ket(0) = hbar omega\/2$ and $H ket(1) = 3 hbar omega\/2$.

a. Find $braket(x', alpha";"t)$ and the expectation values $expval(x)$ and $expval(p)$ in the state $ket(alpha";"t)$ (Schrödinger picture).

We have
$
  ket(alpha";"t) &= 1/sqrt(2) exp((- i H t)/hbar) ket(0) + e^(i delta)/sqrt(2) exp((-i H t)/hbar) ket(1) \
  &= 1/sqrt(2) exp((- i omega t)/2) ket(0) + e^(i delta)/sqrt(2) exp((-3i omega t)/2) ket(1) \
  &= exp((-i omega t)\/2)/sqrt(2) {ket(0) + exp(- i (omega t -delta)) ket(1)} \
  braket(x', alpha";"t)&= exp[(-i omega t)/2]/sqrt(2) {braket(x', 0) + exp[- i (omega t -delta)] braket(x', 1)}
$
from the book
$
  braket(x', 0) = 1/(pi^(1\/4) sqrt(x_0)) exp[-1/2 (x'/x_0)^2]
$
and
$
  braket(x', 1) & = (1/(sqrt(2)x_0)) (x'-x_0^2 dv(, x')) braket(x', 0) \
                & = (sqrt(2) x')/x_0 braket(x', 0)
$
so
$
  braket(x', alpha";"t) = exp((-i omega t)/2)/sqrt(2) [1 + (sqrt(2) x')/x_0 exp[-i (omega t-delta)]] 1/(pi^(1\/4) sqrt(x_0)) exp[-1/2 (x'/x_0)^2]
$
giving
$
  expval(x) & = braket(alpha";"t, x, alpha";"t) \
            & = integral_(-oo)^oo dd(x)' x' abs(braket(x', alpha";"t))^2 \
            & dots = x_0/sqrt(2) cos(omega t- delta)
$
and
$
  expval(p) & = dots = - hbar/(sqrt(2) x_0) sin(omega t - delta)
$


b. Do the same in the Heisenberg picture.

We have
$
  x = x_0/sqrt(2) (a + a^dagger)" and "p = i hbar/(sqrt(2) x_0) (-a+a^dagger)
$
so
$
  x(t) &= x_0/sqrt(2) [a(t) + a^dagger (t)] = x_0/sqrt(2) [a(0) exp(-i omega t) + a^dagger (0) exp(i omega t)] \
  p(t) &= i hbar/(sqrt(2) x_0) [- a(t) + a^dagger (t)] = i hbar/(sqrt(2) x_0) [-a(0) exp(-i omega t) + a^dagger (0) exp(i omega t)]
$
and we obtain
$
  expval(x) & = braket(alpha, x(t), alpha) \
            & = x_0/sqrt(2) cos(omega t - delta)
$
$
  expval(p) & = braket(alpha, p(t), alpha) \
            & = - hbar/(sqrt(2) x_0) sin(omega t - delta)
$
the replacement of $a arrow a(t)$ whne $x arrow x(t)$, works because any such algebraic relation between operators that hold at $t = 0$ will be satisfied for all $t$ in the Heisenberg picture. Since
$
  x (t)^((H)) = cal(U)^dagger x cal(U) = c ( cal(U)^dagger a cal(U) + cal(U)^dagger a^dagger cal(U)) = c ( a(t)^((H))+a^(dagger (H)) (t))
$
these are just way simpler to calculate since we don't need to do some weird integral.

== 2.14
Consider a particle subject to a one-dimensional simple harmonic oscillator potential. Suppose at $t = 0$ the state vector is given by
$
  exp((- i p a)/hbar) ket(0)
$
where $ket(0)$ has $expval(x) = 0 = expval(p)$. Evaluate $expval(x)$ for $t >= 0$ using the Heisenberg picture.

Recall
$
  x(t) = x cos omega t + p/(m omega) sin omega t
$
with $x = x(0)$ and $p = p(0)$.

$
  expval(x)_t &= braket(alpha, x(t), alpha) \
  &= braket(0, exp((i p a)/hbar) (x cos omega t + p/(m omega) sin omega t) exp((-i p a)/hbar), 0) \
  &=^(p "commute" p) braket(0, exp((i p a)/hbar) x exp((-i p a)/hbar), 0) cos omega t + braket(0, p, 0) 1/(m omega) sin omega t \
  &= braket(0, exp((i p a)/hbar)([x,exp((- i p a)/hbar)] + exp((- i p a)/hbar)x), 0) cos omega t \
  &= braket(0, exp((i p a)/hbar) (i hbar) (- i a)/hbar exp((- i p a)/hbar), 0) cos omega t \
  &= a cos omega t
$

== 2.18
Consider the correlation function
$
  C (t) = expval(x(t) x(0))
$
where $x(t)$ is the position operator in the Heisenberg picture. Find it for the ground state of a one-dimensional simple harmonic oscillator.

We have
$
  x(t) = x(0) cos omega t + (p(0))/(m omega) sin omega t
$
in the ground state
$
  C (t) &= braket(0, x(t) x, 0) \
  &= braket(0, (x cos omega t + p/(m omega) sin omega t) x, 0) \
  &= braket(0, x^2, 0) cos omega t + braket(0, p x, 0) (sin omega t)/(m omega)
$
using
$
  x & = sqrt(hbar/(2 m omega)) (a + a^dagger) \
  p & = i sqrt((m hbar omega)/2) (-a + a^dagger)
$
we obtain
$
  C(t) &= hbar/(2 m omega) cos omega t braket(0, (a+a^dagger)(a+a^dagger), 0) \
  &+ i hbar/(2 m omega) sin omega t braket(0, (-a + a^dagger)(a+a^dagger), 0) \
  &= hbar/(2 m omega) { cos omega t braket(0, a a^dagger, 0) - i sin omega t braket(0, a a^dagger, 0)} \
  &= hbar/(2 m omega) { cos omega t - i sin omega t } \
  &= hbar/(2 m omega) e^(- i omega t)
$


== 2.36
Use the WKB method to approximate energy eigenvalues for the one-dimejnsional simple harmonic oscillator potential $V(x) = m omega^2 x^2 \/2$.

We have
$
  E = (m omega^2 x^2)/2 => x_"turn" = plus.minus sqrt((2 E)/(m omega^2))
$
quantization condition
$
  integral_(x_1)^(x_2) sqrt(2m (E - m omega^2 x^2\/2)) dd(x) &= (n+1/2) pi hbar \
  integral_(x_1)^(x_2) sqrt(2 m (m omega^2)\/2 (2 E\/m omega^2 - x^2)) dd(x) &= (n + 1/2) pi hbar \
  m omega integral_(- sqrt(2E\/m omega^2))^(sqrt(2E\/m omega^2)) sqrt(2E\/m omega^2 - x^2) dd(x) &= (n+1/2) pi hbar
$
let
$
  x = sqrt((2 E)/(m omega^2)) sin theta => dd(x) = sqrt((2 E)/(m omega^2)) cos theta dd(theta)
$
so we obtain
$
  (2 E)/(omega) integral_(- pi\/2)^(pi\/2) sqrt(1 - sin^2 theta) cos theta dd(theta) &= (n + 1/2) pi hbar \
  (2 E)/omega integral_(-pi\/2)^(pi\/2) cos^2 theta dd(theta) &= (n+1/2) pi hbar \
  (2 E)/omega pi/2 &= (n+1/2) pi hbar \
  &=> E = (n +1/2) hbar omega
$


== 2.42
We have
$
  Z = integral dd(x', 3) K (bold(x)',t;bold(x)',0) = sum_a' exp[- beta E_a']
$
now
$
  pdv(Z, beta) = - sum_a' E_a' exp[- beta E_a']
$
so
$
  - 1/Z pdv(Z, beta) & = (sum_a' E_a' exp[-beta E_a'])/(sum_a' exp[-beta E_a']) \
  &= (sum_a' E_a' exp[-beta (E_a' - E_0)])/(sum_a' exp[- beta (E_a'-E_0)]) \
  &= (E_0 + E_1 exp[- beta(E_1-E_0)] + dots)/(1 + exp[-beta(E_1-E_0) + dots])
$
now consider $E_i - E_0 > 0$ since $E_0$ is by definition the smallest energy. Thus in the limit $beta -> oo$ every exponential vanishes and we get
$
  lim_(beta -> oo) (- 1/Z pdv(Z, beta)) & = E_0
$

== 2.43
For a free particle
$
  cal(U)^dagger = exp[i/hbar bold(p)^2/(2m) t]
$
so the propagator is
$
  braket(bold(p)''","t, bold(p)'","t_0)_H &= braket(bold(p)'', cal(U) cal(U)^dagger, bold(p)') \
  &= braket(bold(p)'', exp[- i/hbar bold(p)^2/(2m) t] exp[i/hbar bold(p)^2/(2m) t_0], bold(p)') \
  &= braket(bold(p)'', exp[- i/hbar bold(p^2)/(2m) (t-t_0)], bold(p)') \
  &= exp[- i/hbar bold(p)'^2/(2 m) (t-t_0)] braket(bold(p)'', bold(p)') \
  &= exp[- i/hbar bold(p)'^2/(2 m) (t-t_0)] delta^3 (bold(p)''-bold(p)')
$
