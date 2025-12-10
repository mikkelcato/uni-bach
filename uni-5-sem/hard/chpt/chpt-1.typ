//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Introduction
This course is split up in three-four parts with the first covering the mathematical description of crystal, in terms of the lattice and reciprocal lattice. The second part covers elastic waves and phonons in crystals. The third part covers electrons and energy bands. And the last part covers semiconductors and the Fermi surface.

#pagebreak()

= Crystal structure
Everything in solid state physics is basically a crystal. This part of the course covers the mathematical description of these.

== The lattice
An ideal crystal consists of the infinite repetition of identical groups of atoms (the basis). The set to which the basis is attached is called the lattice. This lattice can be defined by three translation vectors $bold(a)_1, bold(a)_2$ and $bold(a)_3$ such that the crystal looks the same at $bold(r)$ and $bold(r)'$ related by
$
  bold(r)' = bold(r) + bold(T) " with " bold(T) = u_1 bold(a)_1 + u_2 bold(a)_2 + u_3 bold(a)_3", " u_1, u_2, u_3 in ZZ
$<lattice-translation>
all $bold(r)'$ give the lattice.

// 8 subjects, 5 min

#let latt1 = cetz.canvas({
  import cetz.draw: *
  // lattice
  let n = 0
  while n < 4 {
    circle((n, 0), radius: 1pt, fill: rgb(black))
    circle((n + .2, 1), radius: 1pt, fill: rgb(black))
    circle((n + .4, 2), radius: 1pt, fill: rgb(black))
    circle((n + .6, 3), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  content((1 + .2, 1), anchor: "north-east", padding: .2em, [$bold(r)$])
  content((3, 2), anchor: "south-west", padding: .3em, [$bold(r)'$])

  line(
    (1 + .2, 1),
    (1 + .4, 2),
    mark: (end: ">"),
    name: "a2",
    stroke: 1.3pt,
    fill: rgb(black),
  )
  line(
    (1 + .2, 1),
    (2 + .2, 1),
    mark: (end: ">"),
    name: "a1",
    stroke: 1.3pt,
    fill: rgb(black),
  )

  content("a2.60%", anchor: "east", padding: .3em, [$bold(a)_2$])
  content("a1.80%", anchor: "north", padding: .3em, [$bold(a)_1$])
})

#figure(
  scale(100%, latt1),
  caption: [example of (primitive) lattice. (here $bold(r)' = bold(r) + 2 bold(a)_1 + bold(a)_2$)],
) <oblique-latt>

To add the basis we can write the position of the $j$th atom in the basis relative to some associated lattice point as
$
  bold(r)_j = x_j bold(a)_1 + y_j bold(a)_2 + z_j bold(a)_3
$ <basis-coord>
it should be clear that we can pick this point such that $0 <= x_j,y_j,z_j <= 1$.

#let basis1 = cetz.canvas({
  import cetz.draw: *
  let a = 3.8
  circle((0, 0), radius: 1pt, fill: rgb(black))
  line((0, -.2), (a, -.2), mark: (end: ">"), name: "a1", fill: rgb(black))
  line(
    (-.2, .2),
    (a / 2 - .2, a * 0.866 + .2),
    mark: (end: ">"),
    name: "a2",
    fill: rgb(black),
  )
  line((0, 0), (a, 0), stroke: (dash: "dashed", thickness: 0.5pt), name: "l1")
  line(
    (0, 0),
    (a / 2, a * 0.866),
    stroke: (dash: "dashed", thickness: 0.5pt),
    name: "l2",
  )
  line(
    (a, 0),
    (a + a / 2, a * 0.866),
    stroke: (
      dash: "dashed",
      thickness: 0.5pt,
    ),
    name: "l3",
  )
  line(
    (a / 2, a * 0.866),
    (a / 2 + a, a * 0.866),
    stroke: (
      dash: "dashed",
      thickness: 0.5pt,
    ),
    name: "l4",
  )

  content("a2.50%", anchor: "east", padding: .4em, [$bold(a)_2$])
  content("a1.50%", anchor: "north", padding: .4em, [$bold(a)_1$])

  circle((a / 3 + a / 6, a * 0.866 / 3), radius: 6pt, fill: teal, name: "c2")
  circle(
    (2 * (a / 3 + a / 6), 2 * a * 0.866 / 3),
    radius: 4pt,
    fill: aqua,
    name: "c1",
  )
  line("c1", "c2")
  line("c1", "l4")
  line("c1", "l3")
  line("c2", "l1")
  line("c2", "l2")

  content(
    "c1",
    anchor: "south-east",
    padding: .3em,
    [$2/3(bold(a)_1+bold(a)_2)$],
  )
  content(
    "c2",
    anchor: "north-west",
    padding: .4em,
    [$1/3 (bold(a)_1+bold(a)_2)$],
  )
})

#figure(scale(100%, basis1), caption: [example of (primitive) basis])

We can define primitive translation vectors $bold(a)_i$ such that every point on the lattice $bold(r)'$ can be written as $bold(r) + bold(T)$. So all points $bold(r)_1, bold(r)_2$ where the crystal looks the same satisfy @lattice-translation. The $bold(a)_i$ also define the primitive cell or unit cell with minimum-volume $V_"cell" = abs(bold(a)_1 dot bold(a)_2 times bold(a)_3)$.
#let latt2 = cetz.canvas({
  import cetz.draw: *
  // lattice
  let n = 0
  while n < 4 {
    circle((0, n), radius: 1pt, fill: rgb(black))
    circle((1, n), radius: 1pt, fill: rgb(black))
    circle((2, n), radius: 1pt, fill: rgb(black))
    circle((3, n), radius: 1pt, fill: rgb(black))
    circle((4, n), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  content((1, 1), anchor: "north-east", padding: .2em, [$bold(r)_1$])
  content((2, 2), anchor: "south-west", padding: .2em, [$bold(r)_2$])

  line(
    (1, 1),
    (1, 2),
    mark: (end: ">"),
    name: "a2",
    stroke: 1.3pt,
    fill: rgb(black),
  )
  line(
    (1, 1),
    (3, 1),
    mark: (end: ">"),
    name: "a1",
    stroke: 1.3pt,
    fill: rgb(black),
  )
  line(
    (1, 1),
    (2, 2),
    mark: (end: ">"),
    name: "T",
    stroke: (dash: "dashed", thickness: 1.3pt),
    fill: rgb(black),
  )

  content("a2.90%", anchor: "east", padding: .3em, [$bold(a)'_2$])
  content("a1.90%", anchor: "north", padding: .3em, [$bold(a)'_1$])
  content("T.100%", anchor: "north-west", padding: .3em, [$bold(T)$])
})
#figure(
  scale(100%, latt2),
  caption: [example of non-primitive (conventional) $bold(a)'_i$. Here $bold(r)_2 eq.not bold(r)_1 + bold(T)$ since the required $bold(T)$ doesn't exist.],
)

The $bold(a)_i$ are not unique so we can have different types of primitive cells, but the number of atoms in a primitive cell, the primitive basis, is the same, and is minimal---there is one lattice point per primitive cell. Any of these will fill all space by repetition of suitable $bold(T)$. Another type of primitive cell is the Wigner-Seitz cell---essentially a Voronoi cell.

=== Symmetries
As described above lattices can be mapped onto themselves using translations $bold(T)$ (by definition). We can also imagine other symmetry operations leaving the lattice invariant. The lattice point group is the collection of these symmetries. This group includes rotations about a point and reflections about a plane (and inversions). For rotations we can find lattices with one-, two-, three-, four-, and sixfold rotational symmetry. Others are not possible since the lattice must have translational symmetry. A lattice type with some specific symmetry(ies) is called a Bravais lattice.

In two dimensions we have five different Bravais lattices. A general lattice is oblique and has one- and twofold rotational symmetry as in @oblique-latt. To get lattices with nice symmetries we impose restrictions on our lattice vectors $bold(a)_1$ and $bold(a)_2$ (see Kittel).

In three dimensions the point symmetry group requires 14 different lattice types. These are all grouped into systems according to the type of cells they have. The one we care about the most is the cubic system with three different lattices, with restrictions: $a_1 = a_2 = a_3$ and $alpha = beta = gamma = 90 degree$. These are the simple-, body-centered-, and the face-centered-cubic (see Kittel).

=== Miller indices
As an aside we are typically interested in certain planes of a given crystal. The orientation of these planes is given by the Miller indices $(h k l)$. These are based on how the lattice vectors $bold(a)_i$ intercept with the plane. To determine the indices we find the intercepts as follows: if the plane intercepts the $bold(a)_1$ axis at the midpoint we take the intercept to be $1/2 bold(a)_1$. Take the same plane to never intercept $bold(a)_2$ and $bold(a)_3$ (intercept at infinity) then the Miller indices would simply be $\(h = (1/2)^(-1), k=0, l=0) = (2 0 0)$. For a negative index we denote it by $\(h overline(k) l)$. The indices $[u v w]$ indicate a direction and can be thought of as a vector. For the cubic system $(h k l) perp [h k l]$ so the direction is normal to the corresponding plane.

== The reciprocal lattice
Crystal structure can be studied using diffraction. For example by using the Bragg condition. Here we suppose incident waves are reflected from parallel planes of atoms with spacing $d$. These will interfere destructively and constructively, with constructive interference when the Bragg condition
$
  2 d sin theta = n lambda " for " n in NN
$
is satisfied. This is a result of lattice periodicity.

As we will see crystal periodicity is naturally described using the reciprocal lattice.

=== Fourier expansion
We know the crystal structure is invariant under $bold(T)$. Then any local property must be invariant under $bold(T)$.

As an example for the electron density $n(bold(r))$ we require $n(bold(r)+bold(T)) = n(bold(r))$. This hints at a Fourier series of the form
$
  n(bold(r)) = sum_bold(G) n_bold(G) e^(i bold(G) dot bold(r))
$
where $bold(G)$ is a set of vectors we have yet to define. We can do the inverse Fourier transform to get the amplitudes
$
  n_bold(G) = V_c^(-1) integral_"cell" dd(V) n(bold(r)) e^(-i bold(G) dot bold(r))
$

To determine $bold(G)$ we proceed by construction and define the reciprocal lattice basis vectors by
$
  bold(b)_1 = 2pi (bold(a)_2 times bold(a)_3)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)";  " bold(b)_2 = 2pi (bold(a)_3 times bold(a)_1)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)";  " bold(b)_3 = 2pi (bold(a)_1 times bold(a)_2)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)
$
given the $bold(a)_i$ are primitive then so are the $bold(b)_i$ and $bold(b)_i dot bold(a)_j = 2 pi delta_(i j)$. Using these we can define the reciprocal lattice by
$
  bold(G) = h bold(b)_1 + k bold(b)_2 + l bold(b)_3
$
with $h, k, l in ZZ$. Now consider
$
  e^(i bold(G) dot bold(T)) = e^(i [h u_1 2 pi + k u_2 2 pi + l u_3 2 pi]) = e^(2pi i[h u_1+k u_2+ l u_3]) = e^(2 m pi i) = 1
$
meaning we have
$
  n(bold(r) + bold(T)) = sum_bold(G) n_bold(G) e^(i bold(G) dot bold(r)) e^(i bold(G) dot bold(T)) = sum_bold(G) n_bold(G) e^(i bold(G) dot bold(r)) = n (bold(r))
$
this is quite nice! And tells us that $bold(G)$ is exactly the set of vectors for which plane waves respect crystal periodicity.

=== Laue condition
We now show that the set $bold(G)$ determine possible reflections. We consider an incoming wavevector $bold(k)$ and outgoing wavevector $bold(k)'$. The scattering vector is then $dd(bold(k), d: Delta) = bold(k)' - bold(k)$. Then
$
  "amplitude at detector" tilde sum_"waves" e^(- i dd(bold(k), d: Delta) dot bold(r))
$
in the continuous limit $sum -> integral dd(r, 3) n(bold(r))$, this gives the definition of the scattering amplitude
$
  F(dd(bold(k), d: Delta)) equiv integral dd(r, 3) n(bold(r)) e^(-i dd(bold(k), d: Delta) bold(r))
$
using the expansion for $n(bold(r))$ we can write
$
  F(dd(bold(k), d: Delta)) &= sum_bold(G) n_bold(G) integral dd(r, 3) e^(i (bold(G)-Delta bold(k)) dot bold(r)) \
  &=^"ideal crystal" sum_bold(G) n_bold(G) (2pi)^3 delta^((3)) (dd(bold(k), d: Delta)-bold(G))
$
this is obviously only non-zero for $bold(G) = dd(bold(k), d: Delta)$. This is called the Laue condition, and shows why $bold(G)$ determine possible reflections. Equivalently we have the Laue equations
$
  bold(a)_i dot Delta bold(k) = 2 pi v_i " with " v_i in ZZ
$
To recover the Bragg condition we consider elastic scattering with $abs(bold(k)')=abs(bold(k))$ giving
$
  bold(k)^2 & = overbracket(abs(bold(k)+bold(G))^2, "by Laue") = bold(k)^2 + bold(G)^2 + 2 bold(k) dot bold(G) \
  0 & = bold(G)^2 + 2 bold(k) dot bold(G) = bold(G)^2 - 2 abs(bold(k)) abs(bold(G)) sin theta \
  abs(bold(G)) & = 2 abs(bold(k)) sin theta = (4 pi)/lambda sin theta \
  lambda & = 2 d_"hkl" sin theta
$
which becomes the general Bragg condition with integer multiples of $bold(G)$. To understand why $abs(bold(G)) = 2 pi\/d_"hkl"$ consider
$
  e^(i bold(G) dot bold(T)) = 1 => bold(G) dot bold(T) = 2 pi times "integer"
$
picking $bold(T) = d hat(n)$ then gives $bold(G) dot (d hat(n)) = 2 pi$ if we take $d$ to be the spacing between planes.

The Bragg condition is only satisfied by $bold(k)$ lying on the edge of any Brillouin zone. This follows from
$
  1/2 bold(G)^2 = - bold(k) dot bold(G)
$
so $bold(k)$ must lie in a plane that is the perpendicular bisector of a line joining the origin to a reciprocal lattice point $bold(G)$. This is exactly how the Brillouin zones are defined.

=== The form factor
When $Delta bold(k) = bold(G)$ we can write
$
  F_bold(G) = N integral_"cell" dd(r, 3) n(bold(r)) e^(-i bold(G) dot bold(r)) equiv N S_bold(G)
$
where $S_bold(G)$ is the structure factor. We can also write
$
  n(bold(r)) = sum_(j-1)^s n_j (bold(r)-bold(r)_j)
$
where $n_j (bold(r)-bold(r)_j)$ is the contribution to $n(bold(r))$ from the $j$th atom (of $s$) in the basis. We obtain
$
  S_bold(G) &= sum_j integral dd(V) n_j (bold(r)-bold(r)_j) e^(-i bold(G) dot bold(r)) \
  &= sum_j e^(-i bold(G) dot bold(r)_j) integral dd(V) n_j (bold(rho)) e^(-i bold(G) dot bold(rho)) \
  &= sum_j f_j e^(- i bold(G) dot bold(r)_j)
$
with $bold(rho) = bold(r)-bold(r)_j$ and $f_j$ being the atomic form factor.

These are useful because $f_j$ contains all information about the electron distribution within a cell, while the rest contains information about where the basis is located.
