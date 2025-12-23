//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *introductory general relativity*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on relativity and cosmology taken during the SDU course---based primarily on lecture notes. All errors are likely mine.
  ],
)

#include "chpt/chpt-1.typ" // Relativity
#include "chpt/chpt-2.typ" // EFE
#include "chpt/chpt-3.typ" // Cosmology

