//**** init-ting
#import "@preview/physica:0.9.5": *
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures
#import "temp.typ": *


#show: note.with(
  title: [
    Solid State Physics
  ],
  name: [Mikkel Kielstrup Holst],
  prof: [Line Jelver],
)

#show figure.caption: emph

#show ref: it => {
  let eq = math.equation
  let el = it.element
  if el != none and el.func() == eq {
    // Override equation references.
    link(el.location(), numbering(
      el.numbering,
      ..counter(eq).at(el.location()),
    ))
  } else {
    // Other references as usual.
    it
  }
}

#include "chpt/chpt-1.typ" // lattice/reciprocal lattice
#include "chpt/chpt-2.typ" // vibrations and phonons
#include "chpt/chpt-3.typ" // electron models
#include "chpt/chpt-4.typ" // fermi surface and semiconductors
