//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Introduction
/*
#cetz.canvas({
  import cetz.draw: *

})
*/

#pagebreak()
= Crystal structure
Everything in solid state physics is basically a crystal, this part of the course describes the basic mathematical description of these.

== The Lattice
An ideal crystal consists of the infinite repetition of identical groups of atoms---this group is referred to as the basis. The set to which the basis is attached is called the lattice. This lattice can be defined by three translation vectors $bold(a)_1, bold(a)_2$ and $bold(a)_3$, such that our crystal looks the same when viewed from $bold(r)$ and $bold(r)'$ related by
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
  caption: [example of (primitive) lattice, here $bold(r)' = bold(r) + 2 bold(a)_1 + bold(a)_2$.],
) <oblique-latt>

To add the basis we can write the position of the $j^"th"$ atom in the basis relative to some associated lattice point as
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

We can define primitive translation vectors $bold(a)_i$ such that every point on the lattice $bold(r)'$ can be written as $bold(r) + bold(T)$---so all points $bold(r)_1, bold(r)_2$ where the crystal looks the same satisfy @lattice-translation. The $bold(a)_i$ also define the primitive cell, or unit cell with minimum-volume $V_"cell" = abs(bold(a)_1 dot bold(a)_2 times bold(a)_3)$.
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
  caption: [example of non-primitive (conventional) $bold(a)'_i$, here $bold(r)_2 eq.not bold(r)_1 + bold(T)$ since the required $bold(T)$ doesn't exist.],
)

The $bold(a)_i$ are not unique so we can have different types of primitive cells, but the number of atoms in a primitive cell, the primitive basis, is the same, and is minimal---there is one lattice point per primitive cell. Any of these will fill all space by repetition of suitable $bold(T)$. Another type of primitive cell is the Wigner-Seitz cell---essentially a Voronoi cell.
== Symmetries and types
=== Point groups
Lattices can be mapped onto themselves by translations $bold(T)$ (the space group) and other symmetry operations. The lattice (crystallographic) point group is the collection of symmetries which when applied to a lattice point returns the lattice. This group includes rotations about an axis and reflections about a plane---an inversion would be a combination of these. For rotations we can find lattices with one-, two-, three-, four-, and sixfold rotational symmetry---others are not possible since our lattice is infinite due to translational symmetry.



=== Two dimensions
A lattice type is called a Bravais lattice, in two dimensions we have five different Bravais lattices. A general lattice is oblique and has one- and twofold rotational symmetry---see @oblique-latt. To get lattices with nice symmetries we impose restrictions on $bold(a)_1$ and $bold(a)_2$. Here we have four special lattices---see @spec-latt.
#let latt-spec-1 = cetz.canvas({
  import cetz.draw: *
  // square lattice
  let n = 0
  while n < 3 {
    circle((2 * n, 0), radius: 1pt, fill: rgb(black))
    circle((2 * n, 2), radius: 1pt, fill: rgb(black))
    circle((2 * n, 4), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  line(
    (0, 4),
    (0, 2),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a1-1",
    stroke: 1.3pt,
  )
  line(
    (0, 4),
    (2, 4),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a2-1",
    stroke: 1.3pt,
  )

  arc(
    "a2-1.20%",
    start: 0deg,
    stop: -90deg,
    radius: 12pt,
    mark: (end: ">", stroke: 1pt),
    name: "arc-1",
  )

  content("a1-1.40%", anchor: "east", padding: .3em, [$bold(a)_1$])
  content("a2-1.40%", anchor: "south", padding: .3em, [$bold(a)_2$])
  content("a2-1.25%", anchor: "north", padding: .3em, [$phi$])
})

#let latt-spec-3 = cetz.canvas({
  import cetz.draw: *
  // square lattice
  let n = 0
  while n < 3 {
    circle((2 * n, 0), radius: 1pt, fill: rgb(black))
    circle((2 * n, 1), radius: 1pt, fill: rgb(black))
    circle((2 * n, 2), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  line(
    (0, 2),
    (0, 1),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a1-1",
    stroke: 1.3pt,
  )
  line(
    (0, 2),
    (2, 2),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a2-1",
    stroke: 1.3pt,
  )

  arc(
    "a2-1.20%",
    start: 0deg,
    stop: -90deg,
    radius: 12pt,
    mark: (end: ">", stroke: 1pt),
    name: "arc-1",
  )

  content("a1-1.40%", anchor: "east", padding: .3em, [$bold(a)_1$])
  content("a2-1.40%", anchor: "south", padding: .3em, [$bold(a)_2$])
  content("a2-1.25%", anchor: "north", padding: .4em, [$phi$])
})

#let latt-spec-2 = cetz.canvas({
  import cetz.draw: *
  // square lattice
  let n = 0
  while n < 4 {
    circle((2 * n, 0), radius: 1pt, fill: rgb(black))
    circle((2 * n + 1, 1), radius: 1pt, fill: rgb(black))
    circle((2 * n, 2), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  line(
    (0, 2),
    (1, 1),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a2-1",
    stroke: 1.3pt,
  )
  line(
    (0, 2),
    (2, 2),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a1-1",
    stroke: 1.3pt,
  )

  line(
    (4, 2),
    (4, 0),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a2-2",
    stroke: 1.3pt,
  )
  line(
    (4, 2),
    (6, 2),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a1-2",
    stroke: 1.3pt,
  )


  arc(
    "a1-2.22%",
    start: 0deg,
    stop: -90deg,
    radius: 12pt,
    mark: (end: ">", stroke: 1pt),
    name: "arc-1",
  )

  content("a2-1.50%", anchor: "east", padding: .3em, [$bold(a)_2$])
  content("a1-1.50%", anchor: "south", padding: .3em, [$bold(a)_1$])

  content("a2-2.50%", anchor: "east", padding: .3em, [$bold(a)_2$])
  content("a1-2.50%", anchor: "south", padding: .3em, [$bold(a)_1$])

  content("a1-2.30%", anchor: "north", padding: .4em, [$phi$])
})

#let latt-spec-4 = cetz.canvas({
  import cetz.draw: *
  // square lattice
  let n = 0
  while n < 3 {
    circle((2 * n + 1, 0), radius: 1pt, fill: rgb(black))
    circle((2 * n, 2), radius: 1pt, fill: rgb(black))
    circle((2 * n + 1, 4), radius: 1pt, fill: rgb(black))
    n = n + 1
  }
  line(
    (0 + 1, 4),
    (0, 2),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a2-1",
    stroke: 1.3pt,
  )
  line(
    (0 + 1, 4),
    (3, 4),
    mark: (end: ">"),
    fill: rgb(black),
    name: "a1-1",
    stroke: 1.3pt,
  )

  arc(
    "a1-1.12%",
    start: 0deg,
    stop: -90deg,
    radius: 12pt,
    mark: (end: ">", stroke: 1pt),
    name: "arc-1",
  )

  content("a2-1.30%", anchor: "east", padding: .3em, [$bold(a)_2$])
  content("a1-1.30%", anchor: "south", padding: .3em, [$bold(a)_1$])
  content("a1-1.19%", anchor: "north", padding: .4em, [$phi$])
})

#subpar.grid(
  figure(
    scale(100%, latt-spec-1),
    caption: [square with $abs(bold(a)_1) = abs(bold(a)_2)$; $phi = 90 degree$],
  ),
  <a>,

  figure(
    scale(100%, latt-spec-4),
    caption: [hexagonal with $abs(bold(a)_1) = abs(bold(a)_2)$; $phi = 120 degree$],
  ),
  <b>,

  figure(
    scale(100%, latt-spec-3),
    caption: [rectangular with $abs(bold(a)_1) eq.not abs(bold(a)_2)$; $phi = 90 degree$],
  ),
  <c>,

  figure(
    scale(100%, latt-spec-2),
    caption: [centered rectangular, with primitive $bold(a)_i$, and unit cell with $phi=90 degree$],
  ),
  <d>,

  columns: (1fr, 1fr),
  caption: [special lattices in two dimensions],
  label: <spec-latt>,
)

=== Three dimensions
In three dimensions the point symmetry group requires 14 different lattice types. The generel lattice is called triclinic, the other 13 are special. These are all grouped into systems according to seven types of cells---see @latt3d.
#figure(
  table(
    columns: (auto, auto, auto),
    table.hline(stroke: 2pt),
    table.header(
      [system], [$hash$ of lattices], [restrictions on conventional cell]
    ),
    stroke: 0pt,
    table.hline(stroke: 1.5pt + gray),
    [triclinic],
    $1$,
    [$a_1 eq.not a_2 eq.not a_3$; $alpha eq.not beta eq.not gamma$],

    [monoclinic],
    $2$,
    [$a_1 eq.not a_2 eq.not a_3$; $alpha = gamma = 90 degree eq.not beta$],

    [orthorhombic],
    $4$,
    [$a_1 eq.not a_2 eq.not a_3$; $alpha eq beta eq gamma eq 90 degree$],

    [tetragonal],
    $2$,
    [$a_1 eq a_2 eq.not a_3$; $alpha eq beta eq gamma eq 90 degree$],

    [cubic], $3$, [$a_1 eq a_2 eq a_3$; $alpha eq beta eq gamma eq 90 degree$],
    [trigonal],
    $1$,
    [$a_1 eq a_2 eq a_3$; $alpha eq beta eq gamma < 120 degree"," eq.not 90 degree$],

    [hexagonal],
    $1$,
    [$a_1 eq a_2 eq.not a_3$; $alpha eq beta eq 90 degree", " gamma = 120 degree$],
    table.hline(stroke: 2pt),
  ),
  caption: [all 14 lattice types],
) <latt3d>

The cubic system is very important, the three lattices are called the simple cubic, body-centered cubic, and the face-centered cubic---what they look like should be obvious. Note that the position of a point within a cell is denoted using @basis-coord.

/*
#let latt-cubic = cetz.canvas({
  import cetz.draw: *
  // cubic lattice
  let n = 0
  while n < 3 {
    circle((5 * n, 0), radius: 1pt, fill: rgb(black))
    circle((5 * n, 2), radius: 1pt, fill: rgb(black))
    line((5 * n, 0), (5 * n, 2))
    line((5 * n, 0), (5 * n + 2, -0.4))
    line((5 * n, 0), (5 * n + 1, 1))
    line((5 * n, 2), (5 * n + 1, 3))
    line((5 * n, 2), (5 * n + 2, 1.6))
    line((5 * n + 1, 3), (5 * n + 3, 2.6))
    line((5 * n + 3, 2.6), (5 * n + 2, 1.6))
    line((5 * n + 3, .6), (5 * n + 2, -0.4))
    line((5 * n + 3, .6), (5 * n + 1, 1))

    circle((5 * n + 2, 0 - 0.4), radius: 1pt, fill: rgb(black))
    circle((5 * n + 2, 2 - 0.4), radius: 1pt, fill: rgb(black))
    line((5 * n + 2, -0.4), (5 * n + 2, 2 - 0.4))

    circle((5 * n + 1, 0 + 1), radius: 1pt, fill: rgb(black))
    circle((5 * n + 1, 2 + 1), radius: 1pt, fill: rgb(black))
    line((5 * n + 1, 1), (5 * n + 1, 3))

    circle((5 * n + 2 + 1, 0 - 0.4 + 1), radius: 1pt, fill: rgb(black))
    circle((5 * n + 2 + 1, 2 - 0.4 + 1), radius: 1pt, fill: rgb(black))
    line((5 * n + 3, 0.6), (5 * n + 3, 2.6))


    n = n + 1
  }
  circle((6.5, 1.2), radius: 2pt, fill: rgb(black), name: "bcc")
  line((6, 1), "bcc")
  line((5, 0), "bcc")
  line((7, -0.4), "bcc")
  line((7, 1.6), "bcc")
  line((8, 0.6), "bcc")
  line((8, 2.6), "bcc")
  line((5, 2), "bcc")
  line((6, 3), "bcc")
  content((6.3, -0.6), anchor: "north", [bcc])
  content((1.3, -0.6), anchor: "north", [sc])
  content((11.3, -0.6), anchor: "north", [fcc])

  circle((11, 0.75), radius: 2pt, fill: rgb(black), name: "1")
  circle((10.6, 1.4), radius: 2pt, fill: rgb(black), name: "2")
  circle((12.6, 1.1), radius: 2pt, fill: rgb(black), name: "3")
  circle((11.4, 2.2), radius: 2pt, fill: rgb(black), name: "4")
  circle((11.4, 0.3), radius: 2pt, fill: rgb(black), name: "5")
  circle((11.9, 1.7), radius: 2pt, fill: rgb(black), name: "6")

  line("1", (10, 0))
  line("1", (10, 2))
  line("1", (12, -.4))
  line("1", (12, 1.6))
  line("2", (10, 0))
  line("2", (10, 2))
  line("2", (11, 1))
  line("2", (11, 3))
  line("3", (13, 0.6))
  line("3", (13, 2.6))
  line("3", (12, -.4))
  line("3", (12, 1.6))
  line("4", (10, 2))
  line("4", (11, 3))
  line("4", (12, 1.6))
  line("4", (13, 2.6))
  line("5", (10, 0))
  line("5", (12, -.4))
  line("5", (13, .6))
  line("5", (11, 1))
  line("6", (11, 1))
  line("6", (11, 3))
  line("6", (13, 0.6))
  line("6", (13, 2.6))
})

#figure(scale(100%, latt-cubic), caption: [all cubic lattices])<cubic-latt>
*/

=== Crystal planes
The orientation of a crystal plane is given by Miller indices. This is done by finding the intercepts with the crystal axes in terms of $a_i$. Then we take the reciprocals and reduce to the smallest three integers having the same ratio, and if an intercept is a infinity we let that index be zero. The result is $(h k l)$, a negative index gets a minus sign like $(h overline(k) l)$. Many planes will be parallel e.g. $(001)$ and $(002)$. We denote equivalent planes by ${h k l}$, so the set of cube faces would be ${1 0 0}$. The indices $[u v w]$ indicates a direction and is basically a vector. In the cubic system $(h k l) perp [h k l]$. Here $expval(u v w)$ is a family of equivalent directions. In this course we only use $(h k l)$ and $[h k l]$.


In hexagonal crystals we introduce a fourth index $i = -(h + k)$ to get the Miller-Bravais indices $(h k i l)$, similarly we have $[U V T w]$ with $ U = 1/3 (2 u - v)", " V = 1/3 (2 v - u)", and" T = -(u + v) $


== Examples
=== Sodium chloride
The lattice is a fcc, with a basis consisting of one $"Na"^+$ and one $"Cl"^-$ seperated by one-half the body diagonal of the cube. Each unit cube contains four of these with $"Cl"^-$ at ${000; 1/2 1/2 0; 1/2 0 1/2; 0 1/2 1/2}$.

=== Cesium chloride
The lattice is a sc, with a basis consisting of one $"Cs"^+$ at $000$ and one $"Cl"^-$ at $1/2 1/2 1/2$. Either of these can however be treated as the center, therefore both have eight opposing neighbors.

=== Packings
There is an infinite number of ways to maximize the packing fraction for identical spheres. One is fcc (ABCABC$dots$), another is the hexagonal close-packed structure (ABABAB$dots$).

#pagebreak()
= Wave Diffraction and the Reciprocal Lattice
== Diffraction & Fourier
Crystal structure can be studied using diffraction, one way is using Bragg's law. Here we suppose incident waves are reflected from parallel planes of atoms in the crystal, these will then interfere destructively and constructively. If parallel lattice planes have spacing $d$, then the path difference is $2 d sin theta$ with $theta$ being measured with respect to the plane. If the path difference is $n lambda$ with $n in NN$ then we have constructive interference
$
  2 d sin theta = n lambda
$
this is a result of the lattice periodicity, note that if the planes were perfectly reflecting then we would only see radiation from the first plane and therefore a perfect reflection, instead what we actually see is the reflection from many, many planes. This also means that we don't get any information about the intensity or basis---the electron distribution---which we want.

We know the crystal structure is invariant under $bold(T)$, so any local property must also be invariant under $bold(T)$---we want the electron density $n(bold(r))$, so we require $n(bold(r)+bold(T)) = n(bold(r))$. This screams Fourier series---so we write $n(bold(r))$ as
$
  n(bold(r)) = sum_bold(G) n_bold(G) exp[i bold(G) dot bold(r)]
$
where $bold(G)$ is a set of vectors such that all $bold(T)$ leave the crystal invariant---reciprocal lattice vectors---$n_bold(G)$ determines the scattering amplitude. We can do the inverse Fourier transform to get the amplitudes
$
  n_bold(G) = V_c^(-1) integral_"cell" dd(V) n(bold(r)) exp[-i bold(G) dot bold(r)]
$

Now we just need to find $bold(G)$---we essentially proceed by construction. The reciprocal lattice axis vectors are given by
$
  bold(b)_1 = 2pi (bold(a)_2 times bold(a)_3)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)",  " bold(b)_2 = 2pi (bold(a)_3 times bold(a)_1)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)",  " bold(b)_3 = 2pi (bold(a)_1 times bold(a)_2)/(bold(a)_1 dot bold(a)_2 times bold(a)_3)
$
if the $bold(a)_i$ are primitive, then so are the $bold(b)_i$. Then
$
  bold(b)_i dot bold(a)_j = 2 pi delta_(i j)
$
and
$
  bold(G) = v_1 bold(b)_1 + v_2 bold(b)_2 + v_3 bold(b)_3
$
with $v_i in ZZ$---typically these are denoted $h k l$. With this definition we have the required invariance
$
  n(bold(r) + bold(T)) &= sum_bold(G) n_bold(G) exp[i bold(G) dot bold(r)] exp[i bold(G) dot bold(T)] \
  &= sum_bold(G) n_bold(G) exp[i bold(G) dot bold(r)]
$
since $exp[i bold(G) dot bold(T)]=1$ by construction. So every crystal has two lattices associated with it the crystal lattice (physical lattice) and the reciprocal lattice (diffraction lattice).

== Conditions
#thm[
  The set of reciprocal lattice vectors $bold(G)$ determines the possible x-ray reflections.
]

We consider two beams scattered from volume elements $bold(r)$ apart---the incoming and outgoing wavevectors are $bold(k)$ and $bold(k)'$ and the difference in phase factors will be $exp[i(bold(k)-bold(k)') dot bold(r)]$. We define the scattering amplitude as
$
  F equiv integral dd(V) n(bold(r)) exp[i(bold(k)-bold(k)') dot bold(r)] = integral dd(V) n(bold(r)) exp[- i Delta bold(k) dot bold(r)]
$
where $bold(k)+Delta bold(k) = bold(k)'$ is the scattering vector. Using the expansion for $n(bold(r))$ we obtain
$
  F = sum_bold(G) integral dd(V) n_bold(G) exp[i (bold(G)-Delta bold(k)) dot bold(r)]
$
so we get the diffraction condition $ Delta bold(k) = bold(G) => F = V n_G $ it can be shown that $F$ is very small when $Delta bold(k) eq.not bold(G)$---so $bold(G)$ determines the possible reflections. Another way to express $Delta bold(k) = bold(G)$ are the Laue equations
$
  bold(a)_i dot Delta bold(k) = 2 pi v_i
$
If the scattering is elastic $hbar bold(omega)$ is conserved so $k^2 = k'^2$, so $bold(k)+bold(G) = bold(k)'$ becomes
$
  2 bold(k) dot bold(G) + G^2 = 0 <=> 2 bold(k) dot bold(G) = G^2
$
this is equivalent to Bragg's law---one can show that the spacing $d(h k l)$ between parallel lattice planes normal to $bold(G) = h bold(b)_1 + k bold(b)_2 + l bold(b)_3$ is $d(h k l) = 2 pi \/abs(bold(G))$, with these the condition becomes
$
  2 (2pi)/lambda sin theta = (2 pi)/d(h k l)
$
or $2 d(h k l) sin theta = lambda => 2 d sin theta = n lambda$. This can also be interpreted in terms of Brillouin zones, which are Wigner-Seitz primitive cells in the reciprocal lattice---this can be done to higher orders, notably $bold(k)$ touching the _edge_ of the Brillouin zones give diffraction.

== Examples of Lattices
=== sc Lattice
We have the primitive translation vectors
$
  bold(a)_1 = a hat(bold(x))",  " bold(a)_2 = a hat(bold(y))",  " bold(a)_3 = a hat(bold(z))
$
the volume of a cell is $bold(a)_1 dot bold(a)_2 times bold(a)_3 = a^3$ so
$
  bold(b)_1 = (2 pi)/(a) hat(bold(x))",  " bold(b)_2 = (2 pi)/a hat(bold(y))",  " bold(b)_3 = (2 pi)/a hat(bold(z))
$

=== bcc Lattice
We have the primitive translation vectors
$
  bold(a)_1 & = a/2 (-hat(bold(x))+hat(bold(y))+hat(bold(z))) \
  bold(a)_2 & = a/2 (hat(bold(x))-hat(bold(y))+hat(bold(z))) \
  bold(a)_3 & = a/2 (hat(bold(x))+hat(bold(y))-hat(bold(z)))
$
the volume is $a^3\/2$, so
$
  bold(b)_1 = (2pi)/a (hat(bold(y))+hat(bold(z)))",  " bold(b)_2 = (2pi)/a (hat(bold(x))+hat(bold(z)))",  " bold(b)_3 = (2pi)/a (hat(bold(x))+hat(bold(y)))
$

=== fcc Lattice
We have the primitive translation vectors
$
  bold(a)_1 = a/2 (hat(bold(y))+hat(bold(z)))",  " bold(a)_2 = a/2 (hat(bold(x))+hat(bold(z)))",  " bold(a)_3 = a/2 (hat(bold(x))+hat(bold(y)))
$
the volume is $a^3\/4$, so
$
  bold(b)_1 & = (2pi)/a (- hat(bold(x))+hat(bold(y))+hat(bold(z))) \
  bold(b)_2 & = (2pi)/a (hat(bold(x))-hat(bold(y))+hat(bold(z))) \
  bold(b)_3 & = (2pi)/a (hat(bold(x))+hat(bold(y))-hat(bold(z)))
$

== Factors
When $Delta bold(k) = bold(G)$ we can write
$
  F_bold(G) = N integral_"cell" dd(V) n(bold(r)) exp[-i bold(G) dot bold(r)] = N S_bold(G)
$
where $S_bold(G)$ is the structure factor---note it's the integral over a single cell. We can also write
$
  n(bold(r)) = sum_(j-1)^s n_j (bold(r)-bold(r)_j)
$
where $n_j (bold(r)-bold(r)_j)$ is the contribution to $n(bold(r))$ from the $j^"th"$ atom (of $s$) in the basis. We obtain
$
  S_bold(G) &= sum_j integral dd(V) n_j (bold(r)-bold(r)_j) exp[-i bold(G) dot bold(r)] \
  &= sum_j exp[-i bold(G) dot bold(r)_j] integral dd(V) n_j (bold(rho)) exp[-i bold(G) dot bold(rho)] \
  &= sum_j f_j exp[- i bold(G) dot bold(r)_j]
$
with $bold(rho) = bold(r)-bold(r)_j$ and $f_j$ being the atomic form factor and $S_bold(G)$ now being called the structure factor of the basis---note $bold(G) dot bold(r)_j = 2 pi (v_1 x_j + v_2 y_j + v_3 z_j)$.

\*Oxford splits this chapter in two, and does some delta-function business but the end result should be the same.

#pagebreak()
= Crystal Binding and Elastic Constants
What holds a crystal together?---in solids the electrostatic interaction between electrons and nuclei is entirely responsible for cohesion. There are four principal types of crystalline binding---van der Waals, ionic, metallic and covalent.

== Inert gases
These form the simplest crystals and the $e^-$ distribution is similar to the free atoms. In this case the atoms pack as tight as possible, so most are fcc. Tiny distortions of the $e^-$ distribution leads to the van der Waals interaction---small fluctuations lead to induced dipoles.

To model this consider two identical harmonic oscillators separated by $R$. Each has charges $plus.minus e$ with separations $x_1$ and $x_2$---these oscillate along $x$, let $p_1$ and $p_2$ be their momenta. Then the Hamiltonian
$
  cal(H)_0 = 1/(2 m) p_1^2 + 1/2 C x_1^2 + 1/(2 m) p_2^2 + 1/2 C x_2^2
$
we assume each has $omega_0$, then $C = m omega_0^2$. Let $cal(H)_1$ be the Coulomb interaction then
$
  cal(H)_1 = e^2/R + e^2/(R+x_1-x_2) - e^2/(R+x_1) - e^2/(R-x_2)
$
in the point-dipole limit $abs(x_1), abs(x_2) << R$ we find to lowest order
$
  cal(H)_1 tilde.equiv - (2 e^2 x_1 x_2)/R^3
$
writing
$
  x_s equiv 1/sqrt(2) (x_1 + x_2)",  " x_a equiv 1/sqrt(2) (x_1 - x_2)
$
we can obtain
$
  x_1 = 1/sqrt(2) (x_s + x_a)",  " x_2 = 1/sqrt(2) (x_s - x_a)
$
and equivalent expressions for $p_1$ and $p_2$---subbing these in the total Hamiltonian is
$
  cal(H) = [1/(2 m) p_s^2 + 1/2 (C-(2 e^2)/R_3) x_s^2] + [1/(2 m) p_a^2 + 1/2 (C+ (2 e^2)/R^3) x_a^2]
$
so have two coupled modes $s$ and $a$, by inspection
$
  omega = [(C plus.minus (2 e^2)\/R^3)/m]^(1\/2) = omega_0 [1 plus.minus 1/2 ((2 e^2)/(C R^3)) - 1/8 ((2 e^2)/(C R^3))^2 + dots]
$
the zero-point energy is $hbar\/2 (omega_s + omega_a)$ which is lower than the uncoupled value $hbar omega_0$ by
$
  dd(U, d: Delta) = 1/2 hbar(dd(omega_s, d: Delta)+dd(omega_a, d: Delta)) = - hbar omega_0 1/8 ((2 e^2)/(C R^3))^2 = - A/R^6
$
this is the van der Waals interaction or London interaction or induced dipole-dipole interaction$dots$ It is the principle attractive interaction in crystals of inert gases---and is solely a product of the dipole-dipole coupling---one can approximate $A tilde hbar omega_0 alpha^2$.

As two atoms get attracted they will begin to overlap, this lead to a repulsive force due to the Pauli exclusion principle---$e^-$ will be forced to higher states increasing the energy and leading to a repulsive interaction $tilde B\/R^12$. Combining these effects one can write the potential energy of two atoms as
$
  U(R) = 4 epsilon [(sigma/R)^12 - (sigma/R)^6]
$
this is the Lennard-Jones potential---other exponential forms $lambda exp(-R\/rho)$ with $rho$ being a measure of the interaction range can also be used and are usually easier to handle analytically.

We can approximate the total energy of a crystal with $N$ atoms by summing over the Lennard-Jones potential
$
  U_"tot" = 1/2 N(4 epsilon) [sum_j ' (sigma/(p_(i j)R))^12 - sum_j ' (sigma/(p_(i j)R))^6]
$
with $p_(i j) R$ being the distance between atom $i$ and any other atom $j$ in terms of the nearest neighbor distance $R$---the sums are known. For fcc we can find the equilibrium distance
$
  dv(U_"tot", R) = 0 => R_0/sigma = 1.09
$
using this one can then obtain
$
  U_"tot" (R_0) = -(2.15)(4 N epsilon)
$
which should hold for all inert gases---this is the cohesive energy when the atoms are at rest---of course this is also very naive and ignores quantum effects which become more evident for smaller atoms.

== Ionic crystals
Ionic crystals consist of positive and negative ions---the ionic bond comes from the electrostatic interaction of oppositely charged ions---as with inert gas atoms we expect that the ions have symmetric charge distributions, but with distortions where they touch.

The long range interaction between ions with $plus.minus q$ is the Coulomb interaction $plus.minus q^2\/r$---the ions will try to balance this with the repulsive interaction between ion cores, this is similar to the inert gas case. The van der Waals interaction is very weak for ionic bonds, instead the main contribution to the binding energy is the Madelung energy. We define
$
  U_i = sum_j ' U_(i j)
$
where we sum over all $j eq.not i$. We write $U_(i j)$ as a central field repulsive potential (instead of the Pauli potential) and a Coulomb interaction
$
  U_(i j) = lambda exp(- r_(i j)/rho) plus.minus q^2/r_(i j)
$
neglecting surface effects we write $U_"tot" = N U_i$ with $N$ being the amount of molecules---$2 N$ ions. We again use $r_(i j) equiv p_(i j) R$, only including the repulsive interaction for nearest neighbors we find
$
  U_(i j) = cases(
    lambda exp(-R\/rho)- q^2\/R #h(10pt) & "nearest",
    plus.minus 1\/p_(i j) q^2\/R & "else"
  )
$
so
$
  U_"tot" = N U_i = N (z lambda e^(-R\/rho) - (alpha q^2)/R)
$
where $z$ is the number of nearest neighbors of any ion and
$
  alpha equiv sum_j ' ((plus.minus))/p_(i j) equiv "Madelung constant"
$
the equilibrium separation can be written as
$
  N dv(U_i, R) = - (N z lambda)/rho e^(-R\/rho) + (N alpha q^2)/R^2 = 0 => R_0^2 e^(-R_0\/rho) = (rho alpha q^2)/(z lambda)
$
using this we obtain
$
  U_"tot" = - (N alpha q^2)/R_0 (1- rho/R_0)
$
with the Madelung energy being
$
  - (N alpha q^2)/R_0
$
for a crystal to be stable we require that $alpha$ is positive. If we take our reference ion to be negative then the plus sign applies to all positive ions---an equivalent definition of $alpha$ is
$
  alpha/R = sum_j ' ((plus.minus))/r_j
$
with $r_j$ being the distance from the $j^"th"$ ion to the reference ion. As an example consider an infinite line of alternating ions---let $R$ be the distance between adjacent ions, then
$
  alpha/R = 2 [1/R - 1/(2 R) + 1/(3 R) - 1/(4 R) + dots] => alpha = 2 (1-1/2+1/3-dots) = 2 ln 2
$
this is obviously way harder for three-dimensional structures.

== Other bonds
Aside from the mentioned bonds we quickly mention covalent bonds, metallic bonds and hydrogen bonds.

Covalent bonds are formed from two $e^-$---one from each atom. These tend to be partly localized between the two atoms and their spins when bonded are antiparallel. Notably this bond has strong directional properties and doesn't fill space as tightly as other bonds---it only allows four nearest neighbors. The binding energy depends on the spin orientation of the two $e^-$ due to the Pauli exclusion principle which will modify the charge distribution in accordance with the spin distribution. This spin-dependent Coulomb interaction is called the exchange interaction. The big point is that orbitals will hybridize, and sometimes this is favorably and a bond will form (bonding)---other times this is not at all favorably and a bond will not form (antibonding).

Metals have a large number of free $e^-$ zooming around---conduction $e^-$---the valence $e^-$ of the atoms become the conduction $e^-$ of the metal. In a metallic bond the energy of these valence $e^-$ is lower as compared with the free atom.

== Elasticity
We briefly discuss elasticity since we'll eventually discuss elastic waves in crystals.

We treat a crystal as a continuous homogeneous medium---this approximation is valid for elastic waves of $lambda$ longer that $tilde 10^(-8) "m"$. We'll only consider infinitesimal strains such that Hooke's law applies. Take three axis defined by $hat(bold(x))_i$ then after a small uniform deformation these get deformed
$
  bold(x)' &= (1 + epsilon.alt_(x x)) hat(bold(x)) + epsilon.alt_(x y) hat(bold(y)) + epsilon.alt_(x z) hat(bold(z)) \
  bold(y)' &= epsilon.alt_(y x) hat(bold(x)) + (1 + epsilon.alt_(y y)) hat(bold(y)) + epsilon.alt_(y z) hat(bold(z)) \
  bold(z)' &= epsilon.alt_(z x) hat(bold(x)) + epsilon.alt_(z y) hat(bold(y)) + (1 + epsilon.alt_(z z)) hat(bold(z))
$
the coefficients $epsilon.alt_(alpha beta)$ define the deformation---even though the old axes were of unit length this is no longer guaranteed. We define the displacement $bold(R)$ of the deformation as
$
  bold(R) equiv bold(r)'-bold(r) = x (bold(x)'-hat(bold(x))) + dots
$
or more generally
$
  bold(R) (bold(r)) = u(bold(r)) hat(bold(x)) + v(bold(r)) hat(bold(y)) + w(bold(r)) hat(bold(z))
$
with
$
  epsilon.alt_(x x) tilde.equiv pdv(u, x)",  " epsilon.alt_(y x) tilde.equiv pdv(u, y)",  etc."
$
we define the strain components by this
$
  e_(x x) equiv epsilon.alt_(x x) = pdv(u, x)
$
similarly for $e_(y y)$ and $e_(z z)$. For the other components we define
$
  e_(x y) equiv bold(x)' dot bold(y)' tilde.equiv epsilon.alt_(y x) + epsilon.alt_(x y) = pdv(u, y) + pdv(v, x)
$

We define the dilation by $V' = bold(x)' dot bold(y)' times bold(z)'$, and we obtain
$
  delta equiv (V'-V)/V tilde.equiv e_(x x)+e_(y y) + e_(z z)
$
we basically ignore terms of order $epsilon.alt^2$ in all of these.
=== Constants
We define the stress as the force acting on a unit area in the solid---we have nine components $X_x,X_y,X_z, dots$, with the capital letter denoting the direction, and the subscript denoting the normal to the plane to which the force is applied---we say that $Y_z = Z_y$ etc. giving us six stress components $X_x, Y_y, Z_z, Y_z, Z_x, X_y$. Assuming Hooke's law holds we can write
$
  e_(x x) = S_(1 1) X_x + S_(1 2) Y_y + S_(1 3) Z_z + S_(1 4) Y_z + S_(1 5) Z_x + S_(1 6) X_y
$
and so on for $e_(y y), e_(z z), e_(y z), e_(z x)$ and $e_(x y)$. Inversely we can write
$
  X_x = C_(11) e_(x x) + C_(12) e_(y y) + C_(13) e_(z z) + C_(14) e_(y z) + C_(15) e_(z x) + C_(16) e_(x y)
$
etc. etc. All the $S_(11) dots$ are called elastic compliance constants, and the $C_(11) dots$ are called elastic stiffness constants or moduli of elasticity.

We can write the energy density as (Hooke's law)
$
  U = 1/2 sum_(lambda=1)^6 sum_(mu=1)^6 tilde(C)_(lambda mu) e_lambda e_mu
$
with $1 -> 6 = {x x, y y, z z, y z, z x, x y}$, with
$
  X_x = pdv(U, e_(x x)) equiv pdv(U, e_1) = tilde(C)_(1 1) e_1 + 1/2 sum_(beta = 2)^6 (tilde(C)_(1 beta) + tilde(C)_(beta 1)) e_beta
$
now we see that
$
  C_(alpha beta) = 1/2 (tilde(C)_(alpha beta) + tilde(C)_(beta alpha)) = C_(beta alpha)
$
i.e. they are symmetric meaning the $36$ constants reduce to $21$.

If our crystal has some symmetry this number reduces further. We claim that for a cubic crystal
$
  U = 1/2 C_11 (e_(x x)^2 + e_(y y)^2 + e_(z z)^2) + 1/2 C_44 (e_(y z)^2 + e_(z x)^2 + e_(x y)^2) + C_12 (e_(y y)e_(z z) + e_(z z) e_(x x) + e_(x x) e_(y y))
$
for a cubic structure we require four three-fold rotation axes, this happens according to
$
  x -> y -> z -> x & "  " -x->z->y->-x \
    x-> z-> -y-> x & "  " -x->y->z->-x
$
all terms are invariant under these in the posed $U$---if a term were odd we could find a rotation which would change the sign since $e_(x y) = - e_(x(-y))$. So we have found that for a cubic crystal we have merely three elastic stiffness constants.

To conclude we can summarize all of this with Hooke's law
$
  bold(sigma) = bold(C) bold(epsilon)
$
with $bold(sigma)$ being all the $X_x dots$ the stress components, and $bold(epsilon)$ being the strain components $e_(x x) dots$. $bold(C)$ is then a tensor containing all the elastic constants.

=== Bulk modulus
Consider the case $e_(x x)=e_(y y)=e_(z z) = 1/3 delta$ then for a cubic crystal
$
  U = 1/6 (C_11 + 2 C_12) delta^2
$
we define the bulk modulus $B$ by
$
  U = 1/2 B delta^2
$
which is equivalent to
$
  B = - V dv(p, V)
$
so for a cubic crystal
$
  B = 1/3 (C_11+2C_12)
$
the compressibility $K$ is defined as $K = B^(-1)$.

=== Elastic waves
For a volume element in a crystal we can obtain
$
  rho pdv(u, t, 2) = pdv(X_x, x)+pdv(X_y, y)+pdv(X_z, z)
$
or plugging stuff in
$
  rho pdv(u, t, 2) &= C_11 pdv(e_(x x), x) + C_12 (pdv(e_(y y), x)+pdv(e_(z z), x)) + C_44 (pdv(e_(x y), y)+pdv(e_(z x), z)) \
  &= C_11 pdv(u, x, 2) + C_44 (pdv(u, y, 2)+pdv(u, z, 2)) + (C_12+C_44) (pdv(v, x, y)+pdv(w, x, z))
$
by symmetry we have similar equations for $v$ and $w$. These equations are hell---we consider one simple solution
$
  u = u_0 exp[i(K x-omega t)]
$
with $u$ being the $x$-component of the particle displacement---with $K=2 pi\/lambda$ and $omega = 2pi nu$. Plugging this guy in gives
$
  omega^2 rho = C_11 K^2
$
so the velocity of a longitudinal wave in the $[1 0 0]$ direction is
$
  v_s = nu lambda = omega/K = (C_11/rho)^(1\/2)
$
take a transverse wave
$
  v = v_0 exp[i(K x- omega t)]
$
and we obtain the velocity of a transverse wave in the $[1 0 0]$ direction
$
  omega^2 rho = C_44 K^2 => v_s = (C_44/rho)^(1\/2)
$
similarly is obtained by $w$. So for $bold(K)$ parallel to $[1 0 0]$ the two shear waves have equal velocities---this is not generally true.

To find the velocities one can also solve a system of equations see the problems---further note that for longitudinal wave the motion is parallel to $bold(K)$ and for transverse motion it is perpendicular to $bold(K)$---in general there are three modes, see problems.


