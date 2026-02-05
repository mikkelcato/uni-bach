#import "@preview/equate:0.3.2": equate
#import "@preview/physica:0.9.8": *
#import "@preview/theoretic:0.3.1"

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
    first-line-indent: 1em,
    spacing: 1.2em,
    leading: 0.65em,
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

  show ref: theoretic.show-ref

  set align(left)
  doc
}

//
#let theorem = theoretic.theorem.with(
  supplemenet: "Theorem",
  kind: "theorem",
  variant: "plain",
  options: (
    body-font: (style: "normal"),
    block-args: (inset: 8pt),
  ),
)

#let definition = theoretic.theorem.with(
  supplement: "Definition",
  kind: "definition",
  variant: "plain",
  options: (
    body-font: (style: "normal"),
    block-args: (inset: 8pt),
  ),
)

#let example = theoretic.theorem.with(
  supplement: "Example",
  kind: "example",
  variant: "plain",
  options: (
    body-font: (style: "normal"),
    block-args: (inset: 8pt),
  ),
)

#let proof = theoretic.proof.with(suffix: $script(square)$)
