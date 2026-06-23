
#import "temp.typ": *

#let d1 = diagram(
  node-stroke: .1em,
  edge-stroke: .1em,
  node-fill: green.lighten(70%),
  spacing: 0.6cm,
  node((-1, 0), radius: .0em, fill: white, stroke: 0em, name: <A>),
  edge("--"),
  node((0, 0), radius: .8em, name: <B>),
  edge("-"),
  node((1, 0), radius: .8em, name: <C>),
  edge("-"),
  node((2, 0), radius: .8em, name: <D>),
  edge("-"),
  node((3, 0), radius: .8em, name: <E>),
  edge("-"),
  node((4, 0), radius: .8em, name: <F>),
  edge("-"),
  node((5, 0), radius: .8em, name: <G>),
  edge("-"),
  node((6, 0), radius: .8em, name: <H>),
  edge("--"),
  node((7, 0), radius: 0em, fill: white, stroke: 0em, name: <I>),
  node((3.5, .5), radius: 0em, name: <a1>),
  node((3.5, -.5), radius: 0em, name: <a2>),
  edge(<a1>, <a2>, "..", stroke: red),
  /*node(enclose: ((0, 0), (6, 0)), fill: white, stroke: 0em, shape: bracket.with(
    dir: top,
    label: $A B$,
    sep: .3em,
  )),*/
  node(enclose: ((4, 0), (6, 0)), fill: white, stroke: 0em, shape: bracket.with(
    label: $B$,
  )),
  node(enclose: ((0, 0), (3, 0)), fill: white, stroke: 0em, shape: bracket.with(
    label: $A$,
  )),
) // example system

#let d2 = diagram(
  spacing: 1.69cm,
  edge-stroke: .1em,
  node((0, 0), name: <a>),
  edge("..", $phi.alt_2$, label-sep: -1.35em),
  node((2, 0), name: <b>),
  edge("-"),
  node((2, 1.5), name: <c>),
  edge("..", $phi.alt_1$, label-sep: .1em),
  node((0, 1.5), name: <d>),
  edge(<a>, <d>, "-"),
  edge(
    <b>,
    <c>,
    "<->",
    shift: .5em,
    $beta$,
    label-sep: -1.2em,
    stroke: .075em,
    mark-scale: 40%,
  ),
  node((-0.5, 0.75), $braket(phi.alt_2, e^(-beta h), phi.alt_1)=$),
) // example geometry for path integral

#let d3 = diagram(
  spacing: 1.69cm,
  edge-stroke: .1em,
  node((0, 0), name: <A>),
  edge("..", [open], label-sep: -1.35em),
  node((2, 0), name: <B>),
  edge("-"),
  node((2, 1.5), name: <C>),
  edge("..", $phi.alt_1$, label-sep: .1em),
  node((0, 1.5), name: <D>),
  edge(<A>, <D>, "-"),
  edge(
    <B>,
    <C>,
    "<->",
    shift: .5em,
    $beta$,
    label-sep: -1.2em,
    stroke: .075em,
    mark-scale: 40%,
  ),
  node((-0.3, 0.75), $ket(Psi)=$),
) // example geometry for state

#let d4 = diagram(
  spacing: 1.69cm,
  edge-stroke: .1em,
  node((0, 0), name: <A>),
  edge("..", [open], label-sep: -1.35em),
  node((2, 0), name: <B>),
  edge("-"),
  node((2, 1.5), name: <C>),
  edge("-", $oo$, label-sep: .1em),
  node((0, 1.5), name: <D>),
  edge(<A>, <D>, "-"),
  node((-0.4, 0.75), $ket(Omega)_"line"=$),
) // example geometry for state

#let double-pen = diagram(
  spacing: 3cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  edge(<A>, <B>, "-"),
  node((1, 1), name: <B>),
  edge(<A>, <C>, "-"),
  node((1, -1), name: <C>),
  node((2, 0), name: <D>),
  edge(<B>, <D>, "--"),
  edge(<C>, <D>, "--"),
  node((3, -1), name: <E>),
  node((3, 1), name: <F>),
  edge(<D>, <E>, "--"),
  edge(<D>, <F>, "--"),
  edge(<F>, <B>, "zigzag"),
  edge(<E>, <C>, "zigzag"),
  node((4, 0), name: <G>),
  edge(<G>, <F>, "-"),
  edge(<G>, <E>, "-"),
  node((3.2, -0.1), [R], radius: 0em),
  node((2.12, -0.7), [I], radius: 0em),
  node((1.1, -0.1), [L], radius: 0em),
  node((3.15, -1.2), $cal(i)^+$, radius: 0em),
  node((3.15, 1), $cal(i)^-$, radius: 0em),
  node((4.15, 0), $cal(i)^0$, radius: 0em),
  node((3.75, -.7), $cal(I)^+$, radius: 0em),
  node((3.75, +.5), $cal(I)^-$, radius: 0em),
  node((2, 1.3), []),
)

#let pen-bh = diagram(
  spacing: 1.69cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  edge(<A>, <B>, "zigzag"),
  node((1.5, 0), name: <B>),
  node((2.5, 1), name: <C>),
  edge(<B>, <C>, "-"),
  node((0, 3.5), name: <D>),
  node((0, 1.5), name: <E>),
  edge(<A>, <D>, "-"),
  edge(<D>, <C>, "-"),
  edge(<B>, <E>, "--"),
  edge(<A>, (1.75, 1.75), "-", stroke: orange),
)

#let mink-pen = diagram(
  spacing: 1.69cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  node((0, 3), name: <B>),
  node((1.5, 1.5), name: <C>),
  edge(<A>, <B>, "-"),
  edge(<A>, <C>, "-"),
  edge(<B>, <C>, "-"),
  node((0, -.2), $cal(i)^+$, radius: 0em),
  node((0, 3), $cal(i)^-$, radius: 0em),
  node((1.75, 1.5), $cal(i)^0$, radius: 0em),
  node((1.05, .5), $cal(I)^+$, radius: 0em),
  node((1, 2.25), $cal(I)^-$, radius: 0em),
  node((0, 3.2), []),
)

#let pen-ads = diagram(
  spacing: 1.69cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  node((0, 3), name: <B>),
  node((3, 0), name: <C>),
  node((3, 3), name: <D>),

  edge(<A>, <B>, "-", [L]),
  edge(<C>, <D>, "-", [R], label-sep: -1.5em),
  edge(<A>, <C>, "zigzag"),
  edge(<B>, <D>, "zigzag"),
  edge(<A>, <D>, "--"),
  edge(<B>, <C>, "--"),
)

#let pen-evap = diagram(
  spacing: 2.12cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  edge(<A>, <B>, "zigzag"),
  node((1.5, 0), name: <B>),
  node((2.8, .7), name: <C>),
  node((0, 3.5), name: <D>),
  node((0, 1.5), name: <E>),
  edge(<A>, <D>, "-"),
  edge(<D>, <C>, "-"),
  edge(<B>, <E>, "--"),
  edge(<A>, (1.75, 1.75), "-", stroke: orange),
  node((1.5, -0.5), name: <1A>),
  edge(<1A>, <C>),
  edge(<1A>, <B>),
  node((1.70, -0.20), name: <1B>),
  edge(
    <B>,
    <1B>,
    "~>",
    stroke: orange,
    shift: -0.2em,
  ),
)

#let kruskal = diagram(
  spacing: 2.12cm,
  edge-stroke: .08em,
  node((0, 0), name: <A>),
  node((0, 3), name: <B>),
  node((3, 0), name: <C>),
  node((3, 3), name: <D>),

  edge(<A>, <D>, "--"),
  edge(<B>, <C>, "--"),
  edge(
    <A>,
    <C>,
    [$r=0$],
    "zigzag",
    bend: -40deg,
    label-pos: 70%,
    label-angle: 15deg,
    label-size: 0.85em,
    label-sep: -1.1em,
  ),
  edge(<B>, <D>, "zigzag", bend: 40deg),

  edge(
    (0, 1.5),
    (3, 1.5),
    "->",
    [X],
    label-pos: 100%,
    label-sep: -1.4em,
    mark-scale: 50%,
  ),
  edge(
    (1.5, 0),
    (1.5, 3),
    "<-",
    [T],
    label-pos: 0%,
    label-sep: -1.4em,
    mark-scale: 50%,
  ),
)

#let semicircs(node, extrude) = {
  import cetz.draw: *
  let (w, h) = node.size
  let r = 0.5 * calc.max(w, h)
  group({
    arc((0, 0), start: 0deg, stop: -180deg, radius: r, anchor: "origin")
  })
}


#let hh-state-r = diagram(
  spacing: 1.69cm,
  edge-stroke: .1em,
  node((1, 0), name: <A>),
  node((0, 1), name: <B>),
  node((3, 0), name: <C>),
  node((2, 1), name: <D>),
  node((4, 1), name: <E>),

  edge(
    <A>,
    <B>,
  ),
  edge(<A>, <C>),
  edge(<A>, <D>, "--"),
  edge(<C>, <D>, "--"),
  edge(
    <C>,
    <E>,
  ),
  edge(<D>, <E>, $phi.alt_R$, "..", stroke: red),
  edge(<D>, <B>, $phi.alt_L$, "-", stroke: blue),
  node(
    enclose: ((0.125, 0.0), (3.875, 2.0)),
    name: <p0>,
    width: 1em,
    shape: semicircs.with(),
    stroke: black,
    fill: gray.lighten(70%),
    layer: 1,
  ),
  node((2.2, 1), name: <a1>, layer: 2),
  node((1.8, 1), name: <a2>, layer: 2),
  edge(
    <a1>,
    <a2>,
    "<-",
    bend: 90deg,
    stroke: black,
    label: $t_E = beta/2$,
    mark-scale: 50%,
    layer: 2,
  ),
)

#let entanglementsurf = diagram(
  node(
    (1.4, 0.5),
    $y$,
    stroke: 1pt,
    radius: 0.1em,
    fill: black,
  ),
  node((2.6, 0.4), $B$),
  node((1.75, 0.75), $x$, stroke: 1pt, radius: 0.1em, fill: black, layer: 1),
  node(
    (2, 0.8),
    $A$,
    stroke: black + 1pt,
    radius: 35pt,
    fill: green.lighten(85%),
  ),
  node(
    enclose: (
      (1.2, 0.2),
      (2.8, 1.8),
    ),
    stroke: (dash: "dashed"),
  ),
)

#let rindlergeometry = diagram(
  edge-stroke: .1em,
  node((-2.5, 0), layer: 1),
  edge("..", stroke: red, label: $phi.alt_L$),
  node((0, 0)),
  edge("..>", label: $x$, label-pos: 100%, stroke: blue, mark-scale: 50%),
  node((2.5, 0)),
  edge((0, 0), (2.5, 0), label: $phi.alt_R$, stroke: 0pt),
  node((0, 0)),
  edge("->", label: $t_E$, label-pos: 100%, stroke: black, mark-scale: 50%),
  node((0, -1.5)),
  node((.3, 0)),
  edge("<-", bend: 90deg, stroke: black, label: $eta_E = pi$, mark-scale: 50%),
  node((-.3, 0)),
  node(
    enclose: (
      (2.5, 0),
      (-2.5, 0),
      (0, +2),
    ),
    inset: -0.1pt,
    fill: gray.lighten(70%),
  ),
)

#let rindlerwedge = diagram(
  edge-stroke: .06em,
  edge((0, 0), (2, 0), "->", mark-scale: 50%, label: $x$, label-pos: 100%),
  edge((0, 0), (0, -2), "->", mark-scale: 50%, label: $t$, label-pos: 100%),
  edge(
    (0, 0),
    (2, -2),
    "--",
    label: $x = plus.minus t$,
    label-pos: 110%,
    label-side: center,
    label-fill: none,
  ),
  edge((0, 0), (2, 2), "--"),
  edge((2, -1.8), (2, 1.8), bend: -35deg, stroke: red),
)

#let rindlerwedge_alt = diagram(
  edge-stroke: .06em,
  edge((0, 0), (2, 0), "->", mark-scale: 50%, label: $x$, label-pos: 100%),
  edge((0, 0), (0, -2), "->", mark-scale: 50%, label: $t$, label-pos: 100%),
  edge(
    (0, 0),
    (2, -2),
    "--",
    label: $x = plus.minus t$,
    label-pos: 110%,
    label-side: center,
    label-fill: none,
  ),
  edge((0, 0), (2, 2), "--"),
)


#let AdSCFT = diagram(
  edge-stroke: 1pt,
  node(
    (0, 0),
    shape: cylinder.with(tilt: 12deg),
    radius: 5em,
    stroke: (paint: black, thickness: 1pt),
  ),
  node((0, .75), $G eq.not 0$),
  node(
    (2, 0),
    shape: cylinder.with(tilt: 12deg),
    radius: 5em,
    stroke: (paint: black, thickness: 1pt),
    fill: gray,
  ),
  node((2, .75), $G = 0$),
  edge((.85, 0), (1.15, 0), "<->", $"AdS"\/"CFT"$),
  edge((-.75, 0), (-.75, -.25), "->", $t$, label-side: left),

  node((0, .5), shape: ellipse, width: 10em, height: 2.25em, stroke: (
    paint: black,
    thickness: 1pt,
  )),
  node(
    (2, -.65),
    shape: ellipse,
    width: 10em,
    height: 2.1em,
    stroke: (
      paint: black,
      thickness: 1pt,
    ),
    fill: white,
  ),
  node((0, .5), layer: 2, name: <1>, radius: 2pt, fill: black),
  node((.62, 0), layer: 2, name: <2>),
  edge(<1>, <2>, "--"),
  node((0, -.62), layer: 2, name: <3>, radius: 2pt, fill: black),
  edge(<2>, <3>, "--"),
)

#let pagecurve = diagram(
  spacing: 1.69cm,
  edge-stroke: 0.1em,
  edge((0, 0), (0, -4), "->", label: $S_R$, label-pos: 100%, mark-scale: 50%),
  edge((0, 0), (4, 0), "->", label: $t$, label-pos: 100%, mark-scale: 50%),
  edge((0, 0), (2, -2), layer: 1),
  edge((2, -2), (4, -3.75), "--", stroke: red),
  edge((4, -3.75), (4.3, -3.75), "--", stroke: red),
  edge((2, -2), (3.75, 0), bend: 10deg),
  edge(
    (0, 0),
    (2, 0),
    "-|",
    label: $tilde t_"Page"$,
    label-side: right,
    label-pos: 100%,
    mark-scale: 70%,
  ),
  edge((0, -3.75), (2, -2), "--", stroke: blue),
)

#let AMPSfig = diagram(
  edge-stroke: black + .1em,
  node-stroke: 0pt,
  edge(
    (1, -1.2),
    (1, 1.2),
    "--",
    label: $r_s$,
    label-pos: 110%,
    label-side: center,
  ),
  node((0.8, 0), $C$),
  node((1.2, 0), $B$),
  node((3, 0), $A$),
  node(
    enclose: (
      (2.75, -1),
      (2.75, 1),
      (3.25, 0),
    ),
    fill: green.lighten(60%),
  ),
  node(
    enclose: (
      (-0.5, -1.25),
      (-0.5, 1.25),
      (1, 0),
    ),
    fill: gray.lighten(50%),
    inset: -.1pt,
  ),
  node(
    enclose: (
      (1.5, -1.25),
      (1.5, 1.25),
      (1, 0),
    ),
    fill: blue.lighten(60%),
    inset: -.1pt,
  ),
)

#let entanglementradiation = diagram(
  spacing: 1.69cm,
  node((-2, .5), name: <A>),
  node((4, -.5), name: <B>),
  node(
    enclose: (<A>, <B>),
    stroke: 1pt,
    shape: parallelogram.with(angle: 60deg),
    fill: green.lighten(60%),
  ),
  node((-1, .2), [BH], layer: 3),
  node((-.5, -.1), name: <A1>, layer: 3),
  node((-.75, .6), name: <A2>, layer: 3),
  edge(<A1>, <A2>, "-", stroke: 1pt),
  node((-1.5, -.1), name: <B1>, layer: 3),
  node((-1.25, .6), name: <B2>, layer: 3),
  edge(<B1>, <B2>, "-", stroke: 1pt),
  node((-.8, .6), name: <c1>, layer: 1),
  node((-1.2, .6), name: <c2>, layer: 1),
  edge(<c1>, <c2>, stroke: green.lighten(60%) + 4pt),
  node((2.25, 0), [radiation], layer: 3),

  node((2, .4), $bullet$, layer: 3),
  node((3, .4), $bullet$, layer: 3),
  node((2.5, .45), $bullet$, layer: 3),
  node((1.5, .45), $bullet$, layer: 3),

  edge((-.8, .6), (1.7, .2), stroke: green.lighten(60%) + 4pt, bend: -50deg),
  edge((-.925, .6), (2.0, .6), stroke: green.lighten(60%) + 4pt, bend: -50deg),
  edge((-1.075, .6), (2.7, .2), stroke: green.lighten(60%) + 4pt, bend: -50deg),
  edge((-1.2, .6), (3.0, .6), stroke: green.lighten(60%) + 4pt, bend: -50deg),
)

#let ERdiagram = diagram(
  spacing: 1.69cm,
  node-fill: gray,
  node-stroke: black + .15em,
  edge-stroke: black + .15em,
  node((0, 0), name: <L>, shape: ellipse, width: 2.5em, height: 12em),
  node((2.75, 0), name: <R>, shape: ellipse, width: 2.5em, height: 12em),
  edge((0, -.6), (2.75, -.6), bend: -55deg),
  edge((0, +.6), (2.75, +.6), bend: 55deg),
)

#let pen-AMPS = diagram(
  spacing: 1.69cm,
  edge-stroke: .1em,
  node((0, 0), name: <A>),
  edge(<A>, <B>, "zigzag"),
  node((1.5, 0), name: <B>),
  node((2.8, .7), name: <C>),
  node((0, 3.5), name: <D>),
  node((0, 1.5), name: <E>),
  edge(<A>, <D>, "-"),
  edge(<D>, <C>, "-"),
  edge(<B>, <E>, "--"),
  //edge(<A>, (1.75, 1.75), "-", stroke: orange),
  node((1.5, -0.5), name: <1A>),
  edge(<1A>, <C>),
  edge(<1A>, <B>),
  node((1.70, -0.20), name: <1B>),

  edge((0, .2), (.85, .2), bend: -2deg, stroke: 1.2pt + black),
  edge(
    (.85, .2),
    (1.4, .1),
    $C$,
    stroke: 1.2pt + red,
    bend: -8deg,
    label-pos: 20%,
    label-sep: 2pt,
  ),
  edge(
    (1.4, .1),
    (1.9, .1),
    $B$,
    stroke: 1.2pt + blue,
    bend: 30deg,
    label-sep: -15pt,
  ),
  edge(
    (1.9, .1),
    (2.8, .7),
    $A$,
    stroke: 1.2pt + green,
    bend: -4deg,
    label-pos: 50%,
    label-sep: +1pt,
  ),
)

#let sketch-AMPS = diagram(
  spacing: 3.38cm,
  edge-stroke: .1em,
  edge((0, 0), (2, 0), "->", mark-scale: 50%, label: [space], label-pos: 100%),
  edge((0, 0), (0, -2), "->", mark-scale: 50%, label: [time], label-pos: 100%),
  edge((0, 0), (1.75, -1.75), "--"),
  edge((0, -1.75), (1.75, -1.75), "zigzag"),
  edge((0, -.9), (1.75, -.9), "-", $t_"Page"$, label-pos: 12%, stroke: gray),
  node((1.3, -1.1), $B$, shape: rect, fill: blue.lighten(50%)),
  node((1.7, -1.1), $A$, shape: rect, fill: green.lighten(50%)),
  node((.9, -1.1), $C$, shape: rect, fill: red.lighten(50%)),
  edge((1.3, -1.1), (1.7, -1.1), "-", stroke: red),
  edge((1.3, -1.1), (.9, -1.1), "-", stroke: red),
)

#let EPRgedanken = diagram(
  spacing: 2.5cm,
  edge-stroke: .1em,
  node-fill: black,
  mark-scale: 50%,
  node((-.2, 0.4), $A$, layer: 3, shape: rect, fill: green.lighten(50%)),
  node((3.2, 0.4), $B$, layer: 3, shape: rect, fill: blue.lighten(50%)),

  node((1.4, 1), name: <1a1>, radius: 2pt),
  node((1.6, 1), name: <2a1>, radius: 2pt),
  edge(<1a1>, <2a1>, "-", stroke: red),

  node((0, .25), name: <1a1l>, radius: 2pt),
  edge(<1a1>, <1a1l>, "-->--"),
  node((3, .25), name: <2a1r>, radius: 2pt),
  edge(<2a1>, <2a1r>, "-->--"),

  node((1.4, 1.1), name: <1b1>, radius: 2pt),
  node((1.6, 1.1), name: <2b1>, radius: 2pt),
  edge(<1b1>, <2b1>, "-", stroke: red),

  node((0, .35), name: <1b1l>, radius: 2pt),
  edge(<1b1>, <1b1l>, "-->--"),
  node((3, .35), name: <2b1r>, radius: 2pt),
  edge(<2b1>, <2b1r>, "-->--"),

  node((1.4, 1.2), name: <1c1>, radius: 2pt),
  node((1.6, 1.2), name: <2c1>, radius: 2pt),
  edge(<1c1>, <2c1>, "-", stroke: red),

  node((0, .45), name: <1c1l>, radius: 2pt),
  edge(<1c1>, <1c1l>, "-->--"),
  node((3, .45), name: <2c1r>, radius: 2pt),
  edge(<2c1>, <2c1r>, "-->--"),

  node((1.4, 1.3), name: <1d1>, radius: 2pt),
  node((1.6, 1.3), name: <2d1>, radius: 2pt),
  edge(<1d1>, <2d1>, "-", stroke: red),

  node((0, .55), name: <1d1l>, radius: 2pt),
  edge(<1d1>, <1d1l>, "-->--"),
  node((3, .55), name: <2d1r>, radius: 2pt),
  edge(<2d1>, <2d1r>, "-->--"),
)

#let smoothness = diagram(
  spacing: 1cm,
  node((0, -1), fill: black.lighten(30%), radius: 20pt, layer: 2),
  node((0, 0.5), fill: black.lighten(30%), radius: 10pt, layer: 2),
  node((3, 0.5), fill: black.lighten(30%), radius: 10pt, layer: 2),
  node((0, 1.1), [Bob]),
  node((3, 1.1), [Alice]),
  node((0, -.4), [Bob]),
  node((2.5, -.7), [Alice]),
  node((2.75, -1), fill: orange, radius: 3pt, layer: 2),
  node((2.5, -1), fill: orange, radius: 3pt, layer: 2),
  node((2.25, -1), fill: orange, radius: 3pt, layer: 2),
  node((2, -1), fill: orange, radius: 3pt, layer: 2),
  node((3, -1), fill: orange, radius: 3pt, layer: 2),
  edge((-.75, -1), (-.75, .5), "-->", bend: -50deg, [After $t_"Page"$]),
)

