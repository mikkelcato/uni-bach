//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Elastic soft matter
Elastic soft matter could be rubber or gel. The main property of these materials is that we can deform them, and after they'll recover they shape---these typically consist of polymers, in solution (cross-links) or not in solution.

== Basic definitions
We can define shear stress $sigma$ and shear strain $gamma$ as
$
  sigma = F/S
$
and
$
  gamma = d/h
$
for the setup on pg. 29. If the material is an ideal elastic then
$
  sigma = G gamma
$
with $G$ being the shear modulus. If the material is an ideal viscous fluid then
$
  sigma = eta dot(gamma)
$
with $eta$ being the viscocity.

If we compress a material isotropically using some pressure $Delta P$, and the response is $V arrow V + Delta V$ then we can write
$
  Delta P = - K (Delta V)/V
$
with $K$ being the bulk modulus.

If we stretch some material uniformly $L arrow L'$ then we define $lambda = L'\/L$ and $epsilon = lambda - 1$. The elongational stress $sigma$ can be written as
$
  sigma = E epsilon
$
with $E$ being Young's modulus. We can write
$
  E = (9 K G)/(3 K + G)
$

The bulk and shear moduli both characterize how strongly a material wants to preserve its shape. We'll assume that the bulk modulus $K$ is infinite, since softness is due to a low shear modulus $G$---in this case $E = 3 G$.

== Polymer chains

== Kuhn's theory

