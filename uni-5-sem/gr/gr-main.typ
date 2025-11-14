//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *general relativity and cosmology*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on relativity and cosmology taken during the SDU course---based primarily on lecture notes. The purpose of the course is to essentially follow Einstein's derivation of his field equations and then apply them to cosmology. For more of the underlying mathematics see my other notes.
  ],
)

#include "chpt/chpt-1.typ" // Relativity
#include "chpt/chpt-2.typ" // EFE
#include "chpt/chpt-3.typ" // Cosmology

