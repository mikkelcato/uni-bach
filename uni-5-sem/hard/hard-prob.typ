//**** init-ting
#import "@preview/physica:0.9.5": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *problems in solid state physics*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Worked problem sets from class.
  ],
)

= Lecture 1
== Tetrahedral angles
We have
$
  bold(a)_1 & = a/2 vecrow(1, 1, -1) \
  bold(a)_2 & = a/2 vecrow(-1, 1, 1)
$
then
$
  bold(a)_1 dot bold(a)_2 & = - a^2/4 \
                          & = (sqrt(3)/2 a)^2 cos theta \
                          & = 3/4 a^2 cos theta \
                          & => cos theta = - 1/3 \
                    theta & tilde.eq 109.47 degree
$

== Indices of planes
We define
$
  bold(a)_1 & = vecrow(a, 0, 0) \
  bold(a)_2 & = vecrow(0, a, 0) \
  bold(a)_3 & = vecrow(0, 0, a)
$
in the other basis
$
  bold(a)'_1 & = a vecrow(1/2, 1/2, 0) \
  bold(a)'_2 & = a vecrow(0, 1/2, 1/2) \
  bold(a)'_3 & = a vecrow(1/2, 0, 1/2)
$
the Miller indices $(100)$ correspond to the $x = 1$ plane in the $bold(a)_i$ basis---so we require the $x$-component to be $1$. Of those listed $2 bold(a)'_1$ and $2 bold(a)'_3$ intersect this plane so we get $(1/2 0 1/2)$.

Similarly $(001)$ is the $z = 1$ plane, and this corresponds with $(0 1/2 1/2)$.

== hcp structure
We look at the tetrahedron with a basis made from $3$ atoms at the lower layer and tip at the atom in the middle at a height of $h = c/2$---and the distance between any two atoms in the base is of course $a$. Let $d$ be the distance from any atom to the middle of the base (point directly below tip), then by Pythagoras
$
  d^2 + (c/2)^2 = a^2
$
similarly taking $d$ to be the hypotenuse of a triangle with one side $a\/2$, since the base is nice the every angle is $2 theta = 60 degree$, so using this triangle we get
$
  cos 30 degree = a/2 1/d = a/(2 d)
$
plugging this in for $d$ we find
$
      a^2 & = a^2/(4 cos^2 30 degree) + c^2/4 \
      c^2 & = a^2 (4 - 1/(cos^2 30 degree)) \
  (c/a)^2 & = 4 - 4/3 \
      c/a & = sqrt(8/3)
$
#pagebreak()
= Lecture 2
== Interplanar separation
Consider a plane $h k l$ in a crystal.

a. Prove that $bold(G) = h bold(b)_1 + k bold(b)_2 + l bold(b)_3$ is perpendicular to the plane.

The plane contains $bold(r)_1$ and $bold(r)_2$ then it contains $Delta bold(r) = bold(r)_2-bold(r)_1$. Then
$
  bold(G) dot Delta bold(r) = 2 pi (h Delta u + k Delta v + l Delta w ) = 2 pi (0) = 0
$


b. Prove the distance between adjacent parallel planes is $d(h k l) = 2pi \/ abs(bold(G))$.


By definition $bold(G) dot bold(r) = 2 pi m$ adjacent plates have
$
  bold(G) dot bold(d) = 2pi => d = (2 pi)/abs(bold(G))
$


c. Show for a sc lattice that $d^2 = a^2\/(h^2+k^2+l^2)$

We have
$
  bold(b)_i = (2 pi)/a hat(bold(x_i))
$
so
$
  abs(bold(G)) = (2 pi)/a sqrt(h^2+k^2+l^2) => d = a/(sqrt(h^2+k^2+l^2))
$

== Hexagonal space lattice
We have
$
  bold(a)_1 & = (sqrt(3)a)/2 hat(bold(x)) + a/2 hat(bold(y)) \
  bold(a)_2 & = - (sqrt(3)a)/2 hat(bold(x)) + a/2 hat(bold(y)) \
  bold(a)_3 & = c hat(bold(z))
$

a. Show the volume is $sqrt(3)a^2 c \/2$.

We have $V_c = abs(bold(a)_1 dot bold(a)_2 times bold(a)_3)$.

$
  bold(a)_1 dot bold(a)_2 times bold(a)_3 &= bold(a)_1 dot ( (sqrt(3) a c)/2 hat(bold(y)) + (a c)/2 hat(bold(x))) \
  &= (sqrt(3) a^2 c)/4 + (sqrt(3)a^2 c)/4 = (sqrt(3) a^2 c)/2
$

b. Find the primitive reciprocal lattice vectors.

$
  bold(b)_1 &= 2 pi/V_c (bold(a)_2 times bold(a)_3) \
  &= 2 pi 2/(sqrt(3) a^2 c) ((sqrt(3) a c)/2 hat(bold(y)) + (a c)/2 hat(bold(x))) \
  &= (2 pi)/a (hat(bold(x))/sqrt(3) + hat(bold(y)))
$
$
  bold(b)_2 & = (4 pi)/(sqrt(3) a^2 c) (bold(a)_3 times bold(a)_1) \
            & = (2 pi)/(a) (- hat(bold(x))/sqrt(3)+hat(bold(y)))
$
$
  bold(b)_3 & = (4 pi)/(sqrt(3) a^2 c) (bold(a)_1 times bold(a)_2) \
            & = (4 pi)/(sqrt(3) a^2 c) (sqrt(3) a^2)/2 hat(bold(z)) \
            & = (2pi)/(c) hat(bold(z))
$

c. Brillouin zone.

== Width of diffraction maximum
a. Show that
$
  abs(F)^2 = (sin^2 [M/2 (bold(a) dot Delta bold(k))])/(sin^2 [1/2 (bold(a) dot Delta bold(k))])
$

we have
$
  F = (1 - exp[-i M (bold(a) dot Delta bold(k))])/(1 - exp[-i (bold(a) dot Delta bold(k))])
$
so
$
  F^* F &= [(1 - exp[i M (bold(a) dot Delta bold(k))])/(1 - exp[i (bold(a) dot Delta bold(k))])][(1 - exp[-i M (bold(a) dot Delta bold(k))])/(1 - exp[-i (bold(a) dot Delta bold(k))])] \
  &= (2 - (exp[i M (bold(a) dot Delta bold(k))] + exp[-i M (bold(a) dot Delta bold(k))]))/(2 -( exp[i(bold(a) dot Delta bold(k))]+exp[-i(bold(a) dot Delta bold(k))])) \
  &= (1 - cos ( M(bold(a) dot Delta bold(k))))/(1 - cos (bold(a) dot Delta bold(k))) \
  &= (sin^2 [M/2 (bold(a) dot Delta bold(k))])/(sin^2 [1/2 (bold(a) dot Delta bold(k))])
$

b. Max at $bold(a) dot Delta bold(k) = 2 pi h$ with $h in ZZ$. $-> Delta bold(k)$ so $bold(a) dot Delta bold(k) = 2pi h + epsilon$ such that $epsilon$ gives first zero in $sin[M/2 (bold(a) dot Delta bold(k))]$. Show that $epsilon = 2 pi\/M$.

We have
$
  0 = sin[M/2 (2 pi h + epsilon)] & = sin(pi h M + (epsilon M)/2) \
                                  & = (-1)^(h M) sin (epsilon M)/2 \
                                  & = sin (epsilon M)/2 => epsilon = (2 pi)/M
$
