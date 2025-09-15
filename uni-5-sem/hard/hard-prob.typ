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
