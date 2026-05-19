#import "@preview/equate:0.3.2": equate
#import "@preview/physica:0.9.8": *
#import "@preview/inknertia:0.1.0": feynman
#import feynman: *
#import "@preview/fletcher:0.5.8" as fletcher: cetz, diagram, edge, node
#import fletcher.shapes: (
  bracket, cylinder, ellipse, parallelogram, rect, trapezium, triangle,
)
#import "@preview/cetz:0.5.0": *
#import "@preview/subpar:0.2.2"

/*
template for main
*/
// note template
#let note(
  title: none,
  name: [],
  prof: [],
  doc,
) = {
  // page
  set page(
    fill: white,
    paper: "a4",
  )
  // general text
  set text(
    size: 12pt,
    fill: black,
    font: "New Computer Modern", // LaTeX font
  )

  v(100pt)
  // title
  set align(center)
  text(size: 30pt, fill: black, title)
  v(100pt)
  align(center)[#name]
  v(10pt)
  datetime.today().display("[day] [month repr:long] [year]")
  v(120pt)

  align(left)[Supervised by #prof]

  v(50pt)
  // abstract
  par(justify: true)[
        #text(size: 12pt, [
             #highlight[resumé.]
            ])
      ]
  pagebreak()

  // headings
  set heading(
    numbering: "1.",
  )
  show heading: set block(below: 1.2em)
  show heading.where(level: 1): set text(18pt, black)
  show heading.where(level: 2): set text(14pt, black)

  show outline.entry.where(
    level: 1,
  ): set block(above: 1.5em)
  show outline.entry.where(
    level: 2,
  ): set block(above: 1.25em)
  outline(depth: 2)
  set par(
    first-line-indent: 1.5em,
    leading: 1.5em,
    justify: true,
  )

  set page(numbering: "1")

  set scale(reflow: true)
  show figure.caption: set text(size: 11pt, style: "italic")
  show figure: set block(spacing: 2em)
  show ref.where(form: "normal"): set ref(supplement: it => {
    if it.func() == figure {
      "Figure"
    }
  })
  set enum(
    indent: 1.5em,
    body-indent: 0.75em,
  )

  show math.equation: set text(font: "New Computer Modern Math")

  show heading.where(level: 1): it => {
    counter(math.equation).update(0)
    it
  }

  set heading(
    numbering: (..x) => {
      let nums = x.pos()
      nums.at(0) = nums.at(0) - 1
      numbering("1.", ..nums)
    },
  )

  set math.equation(
    numbering: (..nums) => {
      let chapter = counter(heading.where(level: 1)).get().first() - 1
      numbering("(1.1)", chapter, nums.pos().first())
    },
  )

  show: equate.with(
    breakable: true,
    sub-numbering: false,
    number-mode: "label",
  )

  pagebreak()
  set align(left)
  doc
}
