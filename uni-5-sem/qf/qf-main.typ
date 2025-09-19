//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *quantum physics*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on quantum physics taken during the SDU course. Based primarily on _Modern Quantum Mechanics_ by Sakurai and other notes.
  ],
)

#include "chpt/chpt-1.typ"
#include "chpt/chpt-2.typ"
