//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *
#import "@preview/mannot:0.3.0": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#text(size: 30pt, strong("Cosmology"))
= Friedmann-Robertson-Walker metric
In the preceding part of the course we developed all the essential tools of relativity and applied them to the canonical example---the black hole. The second part of this course is focused on the application of relativity to cosmology.

== Cosmological standard model
At _large enough_ scales clusters of galaxies etc. can be treated like pointlike particles in a continuous cosmic fluid. This combined with the assumption that there is no preferred point in the Universe leads to the cosmological principle:

_The universe is homogeneous and isotropic._

We also assume there exists some cosmological time, implying the Universe has a rest frame---i.e. at some fixed time $t$ matter is at rest. The motion of matter (the cosmic gas) would then follow this cosmological time. This also implies we can define a set of co-moving coordinates.

These two assumptions: the cosmological principle and the existence of cosmological time, make up the cosmological standard model.
