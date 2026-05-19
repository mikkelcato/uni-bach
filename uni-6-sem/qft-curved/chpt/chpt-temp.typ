#import "@preview/equate:0.3.2": equate
#import "@preview/physica:0.9.8": *
#import "@preview/fletcher:0.5.8" as fletcher: cetz, diagram, edge, node

/*

template for chapters

*/

// note template
#let chpt-note(
  doc,
) = {
  // page
  set page(
    fill: white,
    paper: "us-letter",
    numbering: "1",
  )

  // general text
  set text(
    size: 12pt,
    fill: black,
    font: "New Computer Modern", // LaTeX font
  )

  show math.equation: set text(font: "New Computer Modern Math")

  // headings
  show heading: set block(below: 1.2em)
  show heading.where(level: 1): set text(16pt, black)
  show heading.where(level: 2): set text(14pt, black)
  set heading(
    numbering: "1.",
  )

  set par(
    first-line-indent: 1.5em,
    leading: 1.5em,
    justify: true,
  )
  show: equate.with(
    breakable: true,
    sub-numbering: true,
    number-mode: "label",
  )
  set math.equation(
    numbering: "(1.1)",
    //number-align: bottom,
  )

  set enum(
    indent: 1.5em,
    body-indent: 0.75em,
  )

  show figure.caption: set text(size: 11pt, style: "italic")

  set align(left)
  doc
}
