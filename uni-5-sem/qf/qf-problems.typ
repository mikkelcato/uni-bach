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

== 1.11

== 1.12

== 1.13

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

