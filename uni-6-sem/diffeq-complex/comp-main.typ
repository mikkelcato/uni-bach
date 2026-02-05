//**** init-ting
#import "temp.typ": *


#show: note.with(
  title: [
    PDEs and Complex Analysis
  ],
  name: [Mikkel Kielstrup Holst],
  prof: [Michele Della Morte],
)

= Introduction
Everywhere in physics we describe phenomena by how something changes and typically relate this change to the quantity itself. This means differential equations pop up all over the place. Some examples include the Einstein field equations
$
  R_(mu nu) - 1/2 g_(mu nu) R = 8 pi G T_(mu nu)
$
and the Euler-Lagrange equations
$
  pdv(cal(L), x^mu) - dv(, t) pdv(cal(L), dot(x)^mu) = 0
$
both of which look deceptively simple.

Since differential equations are everywhere we need to be able to solve them. This typically requires computing some very difficult integrals and many are simply impossible. However using techniques from complex analysis some become very easy. These tools are therefore introduced in the latter half of these notes.

#include "chpt/chpt-1.typ" // ODEs
#include "chpt/chpt-2.typ" // PDEs
//#include "chpt/chpt-3.typ" // solution methods for PDEs
//#include "chpt/chpt-4.typ" // complex analysis

