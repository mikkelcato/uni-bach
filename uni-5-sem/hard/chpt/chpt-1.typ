//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Introduction
This course is split up in three-four parts with the first covering the mathematical description of crystal, in terms of the lattice and reciprocal lattice. The second part covers elastic waves and phonons in crystals. The third part covers electrons and energy bands. And the last part covers semiconductors and the Fermi surface.

#pagebreak()
#text(size: 35pt, strong("The Lattice"))

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

