//**** init-ting
#import "@preview/physica:0.9.7": *
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
    Notes on quantum physics taken during the SDU course. Based primarily on _Modern Quantum Mechanics_ by Sakurai. We assume some familiarity with quantum mechanics. All errors are mine.
  ],
)

#include "chpt/chpt-1.typ" // formality
#include "chpt/chpt-2.typ" // dynamics
#include "chpt/chpt-3.typ" // angular momentum
#include "chpt/chpt-4.typ" // symmetries
#include "chpt/chpt-5.typ" // perturbative methods
#include "chpt/chpt-6.typ" // scattering
