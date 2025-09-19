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
    Notes on astrophysics and cosmology taken during the SDU course---follows three part structure; from very big to big to small. Based on Ryden's _Introduction to Cosmology_ and _An Introduction to Modern Astrophysics_ by Carroll and Ostie---supplemented with notes taken during lecture.
  ],
)

//
#include "chpt/chpt-1.typ" // cosmology
#include "chpt/chpt-2.typ"

