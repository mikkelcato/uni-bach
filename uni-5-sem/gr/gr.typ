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
    Notes on relativity and cosmology taken during the SDU course. Based primarily on _Relativity, Gravitation and Cosmology_ by Cheng and other notes.
  ],
)

= Introduction
Relativity is in short the idea that only relative motion is measurable, this can be stated as a symmetry---namely that our equations should be unchanged under coordinate transformations. Special relativity is the symmetry with respect to inertial frames, and general relativity extends this to general frames. For this reason we use the language of tensors later on. For now we'll just review special relativity.

/*
== Symmetries and Transformations
We begin by looking at some familiar symmetries and transformations. Note we can write (using the usual summation convention) $ arrow(V) = V_i hat(e)_i = V'_i hat(e)'_i $ where we do a change of basis ${hat(e)_i} arrow {hat(e)'_i}$. Take a simple rotation
$
  vec(V'_1, V'_2, V'_3) = mat(cos theta, sin theta, 0; - sin theta, cos theta, 0; 0, 0, 1) vec(V_1, V_2, V_3)
$
or expressed more neatly as $V'_i = [bold(R)]_(i j) V_j$, where $[bold(R)]$ is the rotation matrix---note that any vector should transform like this. Using this we can immediately find that if $F_i = m a_i$ then $F'_i = m a'_i$, i.e. the physics is the same since everything transforms nicely.

The most basic transformation we have is the Gallilean transformation between two inertial frames $O arrow O'$, it is given by
$
  x_i arrow x'_i = [bold(R)(alpha,beta,gamma)]_(i j) x_j - v_i t
$
where $arrow(v)$ is the relative velocity and $alpha, beta, gamma$ specify the relative orientations between frames, and time is assumed to be absolute---note that typically $[bold(R)(0,0,0)]_(i j) = delta_(i j)$. Newtonian relativity states that our physical laws are unchanged under this transformation, so it is the symmetry under this transformation.

The problems with this begin when we do electrodynamics, since we know that $c$ is constant. This violates the previous symmetry. Either Maxwell's equations only hold in one inertial frame, or they actually do obey relativity but the transformation between frames is wrong. The second is correct, and Maxwell's equations turn out to be unchanged under the Lorentz transformation. Which for two frames with relative velocity $arrow(v) = v hat(x)$ is given by
$
  x' = gamma (x - v t), #h(15pt) y' = y, #h(15pt) z' = z, #h(15pt) t' = gamma(t - v/c^2 x)
$
where the Lorentz factor is $gamma = (1- beta^2)^(-1\/2)$ with $beta = v\/c$---note that the Lorentz transformation implies that $c$ is the same in all inertial frames, among other things.

Einstein's relativity made this concrete with his two postulates; the principle of relativity and the constancy of the speed of light. From which the Lorentz transformation naturally follows.

A direct result is the relativity of simultaneity, or if $Delta t = 0$ in one frame then we don't necessarily have $Delta t' = 0$. We also have time dilation and length contraction given by $Delta t = gamma Delta t'$ and $Delta x = Delta x' \/ gamma$, with $O'$ being the rest frame and $O$ moving with respect to this frame.

A quantity which is unchanged under the Lorentz transformation is the spacetime interval
$
  Delta s^2 = Delta x^2 - c^2 Delta t^2
$
which can easily be checked. This is directly related to the proper time ($Delta x = 0$) where we have $Delta s^2 = -c^2 Delta tau^2$. Since there only is one rest-frame, then its time interval must be unique and all observers should agree on the value.

*/

== Review of Special Relativity
In special relativity we work in a 4D Minkowski spacetime, where we treat time as a proper coordinate. The objects living in this space are four-vectors $x^mu = vecrow(c t, x, y, z)$. To set up a coordinate system we choose a set of four basis vectors ${bold(e)_mu}$, and we define the metric by
$ g_(mu nu) equiv bold(e)_mu dot bold(e)_nu $
Any four-vector can be written as $A = A^mu bold(e)_mu$, and the scalar product between any two four-vectors becomes
$ A dot B = g_(mu nu) A^mu B^nu $
this notation is nice because it is easy to see when we have scalars, and scalars are always invariant. The most important is
$ s^2 = g_(mu nu) x^mu x^nu $
which we call the metric equation. Another central quantity is the spacetime interval
$
  s^2 = - c^2 t^2 + x^2 + y^2 + z^2
$
this is exactly the metric equation for $ g_(mu nu) = eta_(mu nu) = diag(-1, 1, 1, 1) $ this metric is constant, so Minkowski spacetime is a flat space. In the general theory of relativity curved spacetime is essentially gravity, which would necessarily make the metric non-constant. For this reason special relativity is the no-gravity limit of general relativity.

Lorentz tranformations are the coordinate transformations between two frames moving with a constant velocity with respect to each other in Minkowski spacetime---we have $x^mu arrow x'^mu$ by $x'^mu = tensor(Lambda, mu, -nu) x^nu$. Note that any four-vector must transform like this.

Having defined the position four-vector we can derive other four-vectors. We know that proper time is a scalar since $s^2 = -c^2 tau^2$ is invariant. So we can define the four-velocity by $ U^mu = dv(x^mu, tau) $ which is a four-vector since $dd(x^mu)$ is. This can be related to the ordinary velocity $ U^mu = gamma dv(x^mu, t) = gamma vecrow(c, v^1, v^2, v^3) $ note that we have $eta_(mu nu) U^mu U^nu = - c^2$ for massive particles, and $eta_(mu nu) U^mu U^nu = 0$ for massless particles. To define four-momentum we naturally do $p^mu equiv m U^mu = gamma vecrow(m c, m v^i)$, or $p^mu = vecrow(E\/c, p^i)$ finding the related invariant gives us the famous mass-energy equivalence
$
  E^2 = m^2 c^4 + arrow(p)^2 c^2
$
and the relativistic momentum and energy $p^i = gamma m v^i$, $E = gamma m c^2$.

#pagebreak()
= General Relativity
== Equivalence principle
=== Newtonian field theory
Newton's theory of gravity can be summarized by
$
  arrow(F) = - G_N (m M)/r^2 hat(r) = m arrow(g)
$
where we introduce the gravitational field
$
  arrow(g) = - G_N M/r^2 hat(r)
$
this can be expressed as
$
  integral.cont_S arrow(g) dot dd(arrow(A)) = - 4 pi G_N M
$
where the LHS represents the gravitational flux through any closed surface $S$ and $M$ is the total mass enclosed. As is familiar this can be rewritten as follows
$
  integral div arrow(g) dd(V) = - 4 pi G_N integral rho dd(V)
$
this must hold for any volume, and if we define $arrow(g) equiv - grad Phi$, with $Phi$ the gravitational potential, we get the field equation
$
  laplacian Phi = 4 pi G_N rho
$
from Newton's second law $arrow(F) = m arrow(a)$ we also get
$
  dv(arrow(r), t, 2) = - grad Phi
$
Importantly the Newtonian field theory of gravitation is clearly incompatible with special relativity as time and space are not treated equally---in fact it is a static theory, this reflects the underlying physics admitting an action at a distance description, implying an infinite signal speed.

=== The principle and implications
The principle of equivalence is an empirical principle stating the equivalence of inertial mass and gravitational mass---so we have $m_I = m_G$ this was first mentioned by Galileo.

Another way of stating this is by considering a frame freely falling in some gravitational field. Since all objects experience the same acceleration, gravity would be absent in this frame---so physics in a free-falling frame is equivalent to physics in an intertial frame without gravity. Likewise physics in a nonaccelerating frame with gravity $arrow(g)$ is equivalent to physics in a frame without gravity but accelerating with $arrow(a)=-arrow(g)$. So according to the equivalence principle accelerating reference frames can be treated as frames with gravity. This gives us a nice definition of an inertial frame, which is now a frame with no gravity. This is one of Einstein's big realizations.

As mentioned the equivalence principle is just a restatement of $m_I = m_G$, but this is only the case within mechanics---weak equivalence. But we could also extend the equivalence between gravity and inertia to all other areas of physics---strong equivalence, which we care about in this course. Since the motion of a body in a gravitational field is independent of the properties of the body, Einstein believed that gravity must be a feature of spacetime itself---i.e. curved spacetime. The equivalence principle allows us to say that acceleration of a body can't be caused by gravity, so gravity is not a force, or bodies movie freely in spacetime with gravity.

The equivalence principle nicely shows some physical implications of gravitation. These come about by considering the same situation from within a spaceship in free-fall---no gravitational effects---and from outside watching the spaceship accelerate in a gravitational field. The physics in each case should be the same. One of the cool consequences is the gravitational bending of light.

A similar effect is gravitational redshift which happens when the gravitational field is parallel to the ray direction. In this case we have a receiver at a distance $h$ above the emitter in a downward pointing gravitational field. In the free-fall frame an observer would see no gravitational effects, so they would detect no frequency shift; $(Delta omega)_"ff" = (omega_"rec" - omega_"em")_"ff" = 0$. The observer outside sees the spaceship accelerating, we assume this acceleration starts when the ray is emitted. Since it takes a finite time $Delta t = h\/c$ for the ray to reach the receiver, the receiver will be moving with velocity $Delta u = g Delta t$. So we'd expect a Doppler shift
$ (Delta omega\/omega)_"Doppler" = Delta u \/c $
this is the low-velocity approximation. Since the ray is compressed this shift must be a blueshift, $(Delta omega)_"Doppler" > 0$. But from the observer in free-fall we know the frequency is not changed, meaning the blueshift must be cancelled. This would happen if the gravity present in the outside frame redshifted the light by
$ (Delta omega\/omega)_"gravity" = - Delta u\/c $
this can be rewritten
$
  (Delta omega)/omega =- (g Delta t)/c = - (g h)/c^2 = - (Delta Phi)/c^2 => (omega_"rec"-omega_"em")/omega_"em" = - (Phi_"rec"-Phi_"em")/c^2
$
this is the formula for gravitational frequency shift, note that the order of $"rec"$ and $"em"$ is not important in the leading order case. As an example take a spherical body where the formula takes the form
$
  (Delta omega)/omega = (G_N M)/(c^2 R)
$
the RHS is very small, so this effect is tiny as with many other gravitational effects.

This effect is weird since someone stationary with respect to the emitter, would recieve a different number of waves per unit time than the emitted rate. This is explained by time dilation, so the number of waves doesn't change, but the time unit itself changes under gravity---this is gravitational time dilation. If $omega tilde 1\/dd(tau)$ we can write the frequency shift as
$
  dd(tau_1) = (1 + (Phi_1 - Phi_2)/c^2) dd(tau_2)
$
if our gravitational field is static this can be integrated,
$
  (Delta tau)/tau = (Delta Phi)/c^2
$
this is very different from the usual time dilation which is a statement about relative motion, here even static clocks run at different rates dependent on the potential at their location. Note that the two types are compatible, and the dilation is given by
$
  tau_1 = (1 + (2 Delta Phi)/c^2 - u^2/c^2)^(1\/2) tau_2
$
Both gravitational frequency shift and time dilation have been proven by experiment.

Different clock rates will lead to different speed measurements, and even the speed of light can be different. So gravitational time dilation implies that the vacuum has an effective index of refraction when a gravitational field is present. At a given $r$ with potential $Phi(r)$ the ratio
$
  dv(r, tau) = c
$
is the light speed according to the local proper time, which is a universal constant. Due to gravitational time dilation an observer somewhere else would measure a different speed using a clock located somewhere else. Typically we pick $r_1 = r$ and $r_2 = oo$. With $r_2$ being our reference point where $Phi(oo) = 0$, and $t equiv tau(oo)$ gives the coordinate time. Then we get
$
  dd(tau) = (1 + Phi(r)/(c^2)) dd(t)
$
this implies that
$
  c(r) equiv dv(r, t) = (1 + Phi(r)/c^2) dv(r, tau) = (1 + phi(r)/c^2) c
$
so the speed of light as measured by the distant observer is reduced by gravity. We can define
$
  n(r) equiv c/c(r) tilde.eq 1 - Phi(r)/c^2
$
To reiterate---the velocity of light has not changed, it is still a universal constant, but an observer with time $t$ measured by some distant clock will see that light appears to move slower. This can be used to calculate the bending of light.\*

#pagebreak()
== The geodesic equation
=== Relativity as geometry
Now we begin general relativity proper. According to Einstein's theory all gravitational effects are caused by the underlying curved spacetime---this is motivated by the equivalence principle. General relativity states that matter and energy warps spacetime $g_(mu nu) eq.not eta_(mu nu)$.

According to the equivalence principle we found
$
  dd(tau) = (1 + Phi/c^2)dd(t)
$
how can this be interpreted geometrically? For a proper time interval we have $dd(s^2) = g_(0 0) dd(x^0)dd(x^0)$, and we know $dd(s^2) = -c^2 dd(tau^2)$ for a proper time interval. This gives
$
  (dd(tau))^2 = - g_(0 0) (dd(t))^2
$
so by comparison we find
$
  g_(0 0) = - (1 + Phi/c^2)^2 tilde.eq - (1 + (2 Phi)/c^2)
$
In other words; gravity makes the metric element $g_(0 0)$ deviate from the flat spacetime value of $eta_(0 0)=-1$. Since $g_(0 0)$ is also directly related to the Newtonian gravitational potential $Phi$, we can say that the elements $g_(mu nu) (x)$ are the relativistic gravitational potentials---and of course these are the objects of interest in relativity. Note that in a curved space any small local region can always be described apporoximately as a flat space. Therefore the equivalence principle always allows us to transform gravity away in a local region, and in this region special relativity will be valid---since gravity will be absent.

Any field theory is a two-step process; source particle $arrow^"field equation"$ field $arrow^"e.o.m"$ test particle. We've already mentioned both parts of the Newtonian field theory, and notably it is static, meaning it has no time evolution---in fact it is the static limit of some field theory. What Einstein wanted was to find the relativistic generalization of the Newtonian field theory, this would automatically make the theory dynamic as well---with time and space being equal. Einstein made the leap that the curved spacetime itself is the gravitational field---this was based on his study of the equivalence principle. The effect of the gravitational interaction between two particles can be described as a source giving rise to a curved spacetime which in turn influences the other particle---or put strongly the equivalence principle requires a metric structure of spacetime, and particles follow geodesics in such a curved spacetime. A geodesic is just the shortest possible curve in the curved space, as nature is typically lazy this is the most plausible path. So general relativity is structured like; source $arrow^"e.f.e"$ curved spacetime $arrow^"geodesic eq."$ test particle---space is not a thing, spacetime is just an expression of the relations among things.

What we want to do in relativity is use the Einstein field equations to find the metric $g_(mu nu)$, which then determines the geodesics our particles move along.

=== The equation
We want to show that the curve with extremum length---the geodesic---can be described in terms of the metric. Any curve can be represented by a set of coordinates $x^mu (lambda)$, where $lambda$ is some parameter. As we have seen the metric determines the infinitesimal length $ dd(s) = sqrt(g_(mu nu) dd(x^mu)dd(x^nu)) $ to get a length we integrate
$
  s = integral dd(s) = integral sqrt((dv(s, lambda))^2) dd(lambda) = integral L dd(lambda)
$
with the Lagrangian being
$
  L = sqrt(g_(mu nu) dv(x^mu, lambda) dv(x^nu, lambda)) = L (x, dot(x))
$
where we define $dot(x)^mu equiv dd(x)^mu\/dd(lambda)$. To find the extremum line we use the variational principle
$
  delta s = delta integral L(x,dot(x)) dd(lambda)=0
$
giving the Euler-Lagrange equation
$
  dv(, lambda) pdv(L, dot(x)^mu) - pdv(L, x^mu) = 0
$
it is a fact that using a Lagrangian of the form $ L(x, dot(x)) = g_(mu nu) dot(x)^mu dot(x)^nu $ gives the same equation of motion---this is just much easier to work with. The derivatives are relatively simple to calculate
$
  pdv(L, x^mu) = pdv(g_(sigma rho), x^mu) dot(x)^sigma dot(x)^rho
$
$
  pdv(L, dot(x)^mu) &= pdv(, dot(x)^mu) [g_(sigma nu) dot(x)^sigma dot(x)^nu] \
  &= g_(sigma nu) pdv(dot(x)^sigma, dot(x)^mu) dot(x)^nu + g_(sigma nu) dot(x)^sigma pdv(dot(x)^nu, dot(x)^mu) \
  &= g_(sigma nu) delta_mu^sigma dot(x)^nu + g_(sigma nu) dot(x)^sigma delta_mu^nu = 2 g_(mu nu) dot(x)^nu
$
plugging everything in we find
$
  dv(, lambda) (g_(mu nu) dot(x)^nu) - 1/2 pdv(g_(sigma rho), x^mu) dot(x)^sigma dot(x)^rho = 0
$
which is our equation of motion. We can rewrite this into a neater form
$
  g_(mu nu) dv(x^nu, lambda, 2) + pdv(g_(mu nu), x^sigma)dv(x^sigma, lambda) dv(x^nu, lambda) - 1/2 pdv(g_(sigma rho), x^mu) dv(x^sigma, lambda) dv(x^rho, lambda)=0
$
in the second term the factor $(dd(x^sigma)\/dd(lambda))(dd(x^nu)\/dd(lambda))$ is symmetric under $sigma arrow.l.r nu$, so only the symmetric part of the coefficient, namely
$
  1/2 (pdv(g_(mu nu), x^sigma) + pdv(g_(mu sigma), x^nu))
$
will contribute. This happens since every tensor can be decomposed in a symmetric and antisymmetric part---since
$
  A_(sigma nu) = 1/2 (A_(sigma nu) + A_(nu sigma)) + 1/2(A_(sigma nu)-A_(nu sigma))
$
whenever a tensor gets contracted with a symmetric product only the symmetric part survives, since the antisymmetric part cancels. So we can write
$
  dv(x^nu, lambda, 2) + Gamma^nu_(sigma rho) dv(x^sigma, lambda) dv(x^rho, lambda) = 0
$
where
$
  g_(mu nu) Gamma^nu_(sigma rho) = 1/2 [pdv(g_(sigma mu), x^rho) + pdv(g_(rho mu), x^sigma) - pdv(g_(sigma rho), x^mu)]
$
is the Christoffel symbol or affine connection.

Naturally this should reduce to the Newtonian case in the Newtonian limit of a particle moving with $v << c$ in a static and weak gravitational field. Nonrelativistic speed $dd(x^i)\/dd(t) << c$ implies $dd(x^i) << c dd(t)$ so
$
  dv(x^i, lambda) << c dv(t, lambda) = dv(x^0, lambda)
$
keeping only the dominant term in the double sum then gives
$
  dv(x^mu, lambda, 2) + Gamma_(0 0)^mu dv(x^0, lambda) dv(x^0, lambda) = 0
$
for a static field $dd(g_(mu nu), d: partial) \/dd(x^0, d: partial) = 0$ since all time derivatives should vanish, then the Christoffel symbol becomes
$
  g_(mu nu) Gamma_(0 0)^mu = - 1/2 pdv(g_(0 0), x^nu)
$
for a weak field we assume the metric is close to $eta_(mu nu)$ meaning
$
  g_(mu nu) = eta_(mu nu) + h_(mu nu)
$
with $h_(mu nu) << 1$ being a small correction field. Since $eta_(mu nu)$ is a constant metric we have $dd(g_(mu nu), d: partial) \/ dd(x^sigma, d: partial) = dd(h_(mu nu), d: partial)\/dd(x^sigma, d: partial)$, so to leading order
$
  eta_(mu nu) Gamma_(0 0)^mu = -1/2 pdv(h_(0 0), x^nu)
$
which has for static $h_(0 0)$ gives
$
  - Gamma_(0 0)^0 = -1/2 pdv(h_(0 0), x^0) = 0",  " Gamma_(0 0)^i = -1/2 pdv(h_(0 0), x^i)
$
plugging this back in we get for $mu = 0$
$
  dv(x^0, lambda) = "constant"
$
and for $mu = i$
$
  dv(x^i, lambda, 2) + Gamma_(0 0)^i dv(x^0, lambda) dv(x^0, lambda) = (dv(x^i, (c t), 2) + Gamma_(0 0)^i) (dv(x^0, lambda))^2 = 0
$
where we use $dd(x^i)\/dd(lambda) = (dd(x^i)\/dd(x^0))(dd(x^0)\/dd(lambda))$ and since $dd(x^0)\/dd(lambda)$ is constant we can say
$
  dv(x^i, lambda, 2) = (dv(x^i, x^0, 2))(dv(x^0, lambda))^2
$
this is just the product rule, and we end up with
$
  dv(x^i, (c t), 2) - 1/2 pdv(h_(0 0), x^i) = 0
$
by comparison with the Newtonian equation of motion $dd(arrow(r), 2)\/dd(t, 2) = - grad Phi$ we find $h_(0 0) = - 2 Phi\/c^2$ and we recover
$
  g_(0 0) = - (1 + (2 Phi)/c^2)
$
so this is quite nice. We also get a new way to characterize a weak field $abs(h_(0 0)) << abs(eta_(0 0)) => abs(Phi\/c^2) << 1$.

With this curved spacetime description we have just found
$
  dd(tau) = sqrt(- g_(0 0))dd(t)", with" g_(00) = - (1 + (2Phi)/c^2)
$
before moving on we want to show how frequency shift results from this. Consider two wavefronts emitted $dd(t_"em")$ apart. In a static field this time interval is maintained until they are received $dd(t_"rec")$. So we have $dd(t_"em") = dd(t_"rec")$. We also have $omega = 1\/dd(tau)$ so we get
$
  omega_"rec"/omega_"em" &= dv(tau_"em", tau_"rec") = (sqrt(-(g_00)_"em") dd(t_"em"))/(sqrt(-(g_00)_"rec") dd(t_"rec")) = ((1+2(Phi_"em"\/c^2) ) / (1+2(Phi_"rec"\/c^2)))^(1\/2) \
  &= 1 + (Phi_"em"-Phi_"rec")/c^2 + cal(O)(Phi^2\/c^4)
$
and we have recovered our previous result.

#pagebreak()
== Some calculus


#pagebreak()
== Covariance principle

#pagebreak()
= Cosmology
