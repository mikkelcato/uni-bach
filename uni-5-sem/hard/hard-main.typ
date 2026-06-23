//**** init-ting
#import "@preview/cetz:0.4.1" // drawings
#import "@preview/subpar:0.2.2" // subfigures
#import "../temp.typ": *


#show: lect-note.with(
  title: [
    Solid State Physics
  ],
  name: [Mikkel Kielstrup Holst],
  prof: [Line Jelver],
)

#include "chpt/chpt-1.typ" // lattice/reciprocal lattice
#include "chpt/chpt-2.typ" // vibrations and phonons
#include "chpt/chpt-3.typ" // electron models
#include "chpt/chpt-4.typ" // fermi surface and semiconductors
