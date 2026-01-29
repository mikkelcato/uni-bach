//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *astrophysics and cosmology*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on astrophysics and cosmology taken during the SDU course. Based on Ryden's _Introduction to Cosmology_ and _An Introduction to Modern Astrophysics_ by Carroll and Ostie. This course was very ragged as are these notes. All errors are mine. ],
)

//
#include "chpt/chpt-1.typ" // cosmology
#include "chpt/chpt-2.typ" // structure and formation
#include "chpt/chpt-3.typ" //stellar structure

