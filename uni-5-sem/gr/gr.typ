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
Relativity is in short the idea that only relative motion is measurable, this can be stated as a symmetry---namely that our equations should be unchanged under coordinate transformations. Special relativity is the symmetry with respect to inertial frames, and general relativity extends this to general frames.

#pagebreak()
= Special Relativity
== Galileo and Lorentz
If we consider a system of point particles interacting gravitationally we can write
$
  m_N dv(arrow(r)_N, t, 2) = G sum_M (m_N m_M (arrow(r)_M-arrow(r)_N))/abs(arrow(r)_M-arrow(r)_N)^3
$
as the force on the particle $N$. This equation is invariant under the Galilean transformations
$
  arrow(r)' = R arrow(r) + arrow(v) t + arrow(d)",  " t' = t + t_0
$
where $R$ is a member of $O(3)$---a rotation group---the other symbols have the obvious meanings. The Galilean transformations evidently form a $10$-parameter group, called the Galilean group. Newton was confused as to why there are many other transformations which do not leave Newtonian physics invariant. Newton's laws only hold in frames related by the Galilean transformations---the inertial frames. Newton wanted to know what made these frames special---he believed there must be some absolute space and that inertial frames are at rest or move with uniform velocity within it.

The problems start for real with Maxwell's electromagnetism, which is not invariant under Galilean transformations---e.g. because the speed of light is a universal constant. Initially it was believed that light propagated in the ether with electromagnetism only being valid in the frame at rest with respect tothe ether---this was disproven by Michelson and Morley. Einstein solved this problem by replacing Galilean invariance with Lorentz invariance, which correspond to the Lorentz transformations
$
  X'^alpha = tensor(Lambda, +alpha, -beta) X^beta + a^alpha
$
$tensor(Lambda, +alpha, -beta)$ is defined such that $tensor(Lambda, +alpha, -gamma) tensor(Lambda, +beta, -delta) tensor(eta, -alpha, -beta) = tensor(eta, -gamma, -delta)$ with $eta_(alpha beta)$ being the usual Minkowski metric defined as $eta_(alpha beta) = g_(alpha beta) equiv e_alpha dot e_beta$---proceeding in this manner we derive everything from this symmetry, which is essentially the cleanest way to do it. Since $a^alpha$ is constant we find
$
  dd(x'^alpha) = tensor(Lambda, +alpha, -gamma) dd(x^gamma)
$
we define the proper time (with $c equiv 1$)
$
  dd(tau^2) equiv dd(t^2)-dd(arrow(x)^2) = - eta_(alpha beta) dd(x^alpha) dd(x^beta)
$
this is from $dd(s^2) = -c^2 dd(t^2)+dd(arrow(x)^2)$ with $dd(arrow(x))=0 => dd(t)=dd(tau)$ giving $dd(s^2) = -c^2 dd(tau^2)$. This is invariant
$
  dd(tau'^2) &= - eta_(alpha beta) dd(x'^alpha) dd(x'^beta) \
  &= - eta_(alpha beta) tensor(Lambda, +alpha, -gamma) tensor(Lambda, +beta, -delta) dd(x^gamma) dd(x^delta) \
  &= - eta_(gamma delta) dd(x^gamma) dd(x^delta) \
  &= dd(tau^2)
$
This gives us a universal speed of light, since we have for light
$
  abs(dv(arrow(x), t)) = c = 1 => dd(tau) = 0
$
so light travels along null-geodesics---as will be explained later. Since $dd(tau)$ is invariant we have $dd(tau')=0$ so $ abs(dv(arrow(x)', t')) = 1 $ i.e. the speed of light is the same in all coordinates.

These transformations form the Lorentz group, typically we are interested in the proper Lorentz groups, which is a subgroup satisfying
$
  tensor(Lambda, +0, -0) >= 1",  " det Lambda =+1
$
this excludes non-fundamental symmetries corresponding to non-physical transformations. Further taking $a^alpha = 0$ gives us the homogeneous proper Lorentz groups, this has a subgroup with
$
  tensor(Lambda, +i, -j) = R_(i j)",  " tensor(Lambda, +i, -0)=tensor(Lambda, +0, -i) = 0",  " tensor(Lambda, 0, -0) = 1
$
so the difference between this and the Galilean group are the boosts. Let's assume an observer $cal(O)$ sees a particle at rest, while another observer $cal(O)'$ sees it having velocity $arrow(v)$, then
$
  dd(x'^alpha)=tensor(Lambda, alpha, -beta) dd(x^beta)
$
but we know $dd(arrow(x))=0$ so
$
  dd(x'^i) = tensor(Lambda, i, -0) dd(t)",  " dd(t') = tensor(Lambda, 0, -0) dd(t)
$
combining these give
$
  v^i = dv(x'^i, t') => tensor(Lambda, i, -0) = v^i tensor(Lambda, 0, -0)
$
and from $eta_(gamma delta) = tensor(Lambda, alpha, -gamma) tensor(Lambda, beta, -delta) eta_(alpha beta)$ we get
$
  -1 = tensor(Lambda, alpha, -0) tensor(Lambda, beta, -0) eta_(alpha beta) = sum_i (tensor(Lambda, i, -0))^2 - (tensor(Lambda, 0, -0))^2
$
which upon plugging in the previous gives
$
  tensor(Lambda, 0, -0) = gamma",  " tensor(Lambda, i, -0) = gamma v^i
$
the other elements are not uniquely determined, but we can pick
$
  tensor(Lambda, i, -j) = delta_(i j) + v_i v_j (gamma - 1)/v^2",  " tensor(Lambda, 0, -j) = gamma v_j
$

== Force, energy and momentum
We can generalize Newton's second law with
$
  f^alpha = m dv(x^alpha, tau, 2)
$
if a particle is a rest then $dd(tau)=dd(t)$ and this reduces to the usual non-relativistic force, and since $dd(x'^alpha)=tensor(Lambda, alpha, -beta)dd(x^beta)$ with $dd(tau)=dd(tau')$ this implies that $f'^alpha = tensor(Lambda, alpha, -beta) f^beta$, so it is a four-vector---this is nice because we can construct the four-force from the non-relativistic force.

Similarly we can define four-momentum
$
  p^alpha = m dv(x^alpha, tau)
$
we can then write Newton's law with a net force as
$
  dv(p^alpha, tau) = f^alpha
$
this is nice because it is covariant. Using $dd(tau) = sqrt(dd(t^2)-dd(arrow(x)^2)) = sqrt(1-v^2) dd(t) = dd(t)\/gamma$ and $arrow(v) = dd(arrow(x))\/dd(t)$ we can find
$
  p^0 = m gamma equiv E",  " arrow(p)=m gamma arrow(v)
$
in the low-velocity limit these become
$
  arrow(p) = m arrow(v) + cal(O)(v^3)",  " E = m + 1/2 m v^2 + cal(O)(v^4)
$
By linearity of the Lorentz transformation, if $p^alpha$ is conserved then $p'^alpha$ is conserved.

== Four-vectors and Tensors
Any four-vector transforming as
$
  V^alpha arrow V'^alpha = tensor(Lambda, alpha, -beta) V^beta", for" x^alpha arrow x'^alpha = tensor(Lambda, alpha, -beta)x^beta
$
is contravariant, and any four-vector transforming inversely
$
  U_alpha arrow U'_alpha = tensor(Lambda, -alpha, beta) U_beta
$
with $tensor(Lambda, -alpha, beta) equiv eta_(alpha gamma) eta^(beta delta) tensor(Lambda, gamma, -delta)$ and $eta^(beta delta)=eta_(beta delta)$ with $eta^(beta delta) eta_(alpha delta) = tensor(delta, beta, -alpha)$ is covariant. These are inverses, meaning
$
  tensor(Lambda, -alpha, gamma) tensor(Lambda, alpha, -beta) = eta_(alpha delta) eta^(gamma epsilon) tensor(Lambda, delta, -epsilon) tensor(Lambda, alpha, -beta) = eta_(epsilon beta) eta^(gamma epsilon) = tensor(delta, -beta, gamma)
$
this means that contractions become invariant
$
  U'_alpha V'^alpha = tensor(Lambda, -alpha, gamma) tensor(Lambda, +alpha, -beta) U_gamma V^beta = tensor(delta, -beta, gamma) U_gamma V^beta = U_beta V^beta
$
which is very nice. Every contravariant four-vector has a covariant buddy and vice versa given by
$
  V_alpha equiv eta_(alpha beta) V^beta",  " U^alpha equiv eta^(alpha beta) U_beta
$
so we can raise and lower indices,
$
  eta^(alpha beta) V_beta = eta^(alpha beta) eta_(beta gamma) V^gamma = V^alpha
$
we can check that $V_alpha$ is actually covariant,
$
  V'_alpha = eta_(alpha beta) V'^beta = eta_(alpha beta) tensor(Lambda, beta, -gamma)V^gamma = eta_(alpha beta) eta^(gamma delta) tensor(Lambda, beta, -gamma) V_delta = tensor(Lambda, -alpha, delta)V_delta
$
similarly we can check that $U^alpha$ is contravariant,
$
  U'^alpha = eta^(alpha beta) U'_beta = eta^(alpha beta) tensor(Lambda, -beta, gamma) U_gamma = eta^(alpha beta) eta_(gamma delta) tensor(Lambda, -beta, gamma) U^delta = tensor(Lambda, alpha, -delta) U^delta
$

Tensors are essentially generalized four-vectors, so they have more indices, where contra- and covariant indices transform as one would expect,
$
  tensor(T, gamma, -alpha beta) arrow tensor(T', gamma, -alpha beta) = tensor(Lambda, gamma, -delta) tensor(Lambda, -alpha, epsilon) tensor(Lambda, -beta, rho) tensor(T, delta, -epsilon rho)
$
indices within a tensor can also be contracted $tensor(T, alpha gamma) equiv tensor(T, alpha, -beta, gamma beta)$,
$
  tensor(T', alpha gamma) = tensor(T', alpha, -beta, gamma beta) &= tensor(Lambda, alpha, -delta) tensor(Lambda, -beta, epsilon) tensor(Lambda, gamma, -rho) tensor(Lambda, beta, -kappa) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(delta, -kappa, epsilon) tensor(T, delta, -epsilon, rho kappa) \
  &= tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, -epsilon, rho epsilon) = tensor(Lambda, alpha, -delta) tensor(Lambda, gamma, -rho) tensor(T, delta, rho)
$
we can also make many other construction which are also tensors, e.g. by taking derivatives
$
  tensor(T, -alpha, beta gamma) equiv pdv(, x^alpha) tensor(T, beta gamma)
$
is a tensor when $tensor(T, beta gamma)$ is a tensor---this is more obvious if we write,
$
  partial_alpha equiv pdv(, x^alpha)
$
which is covariant. The direct product of two tensors is a tensor,
$
  tensor(T, alpha, -beta, gamma) equiv tensor(A, alpha, -beta) tensor(B, gamma)
$
the linear combination of tensors is a tensor,
$
  tensor(T, alpha, -beta) equiv a tensor(R, alpha, -beta) + b tensor(S, alpha, -beta)
$
the Minkowski metric is by definition a tensor, implying that $tensor(delta, alpha, -beta)$ is a tensor---this also implies that raising or lowering indices $tensor(T, -alpha, delta, -gamma) equiv eta^(delta beta) tensor(T, -alpha beta gamma)$ conserves tensor-ness. The Levi-Civita tensor $epsilon^(alpha beta gamma delta)=-epsilon_(alpha beta gamma delta)$ is a tensor. The zero tensor $tensor(T, alpha, -beta)=0$ is always zero, $tensor(T', alpha, -beta)=0$, since $tensor(T, alpha, -beta)=tensor(S, alpha, -beta) => tensor(T', alpha, -beta)=tensor(S', alpha, -beta)$.

#pagebreak()
= The equivalence principle
== Newtonian field theory
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

== The principle and implications
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
= Tensors in relativity
In special relativity the metric was the Minkowski metric $eta_(alpha beta)$, in general we have some other metric $g_(alpha beta) = e_alpha dot e_beta$ which has an inverse $g^(alpha beta) = e^alpha dot e^beta$---we also have $e_alpha e^beta = tensor(delta, -alpha, beta)$, implying $g_(alpha beta) g^(beta gamma) = tensor(delta, -alpha, gamma)$.

Before we had $dd(s^2) = eta_(alpha beta) dd(x^alpha) dd(x^beta)$, now this is just $ dd(s^2) = g_(alpha beta) dd(x^alpha) dd(x^beta) $ which is sometimes used as the definition of the metric---relativity requires that our equations are covariant under any transformation that leaves this invariant. All the previous stuff regarding contractions is of course still valid, now we just use the metric $g_(alpha beta)$ instead of $eta_(alpha beta)$, and we allow arbitrary coordinate transformations.

As opposed to special relativity we now deal with position dependent local transformations---this is because our basis vectors now depend on where we are in space, so $g_(alpha beta) equiv g_(alpha beta) (x)$, this implies that our transformations must also be coordinate dependent. Since space is locally flat it is fair to demand that a transformation $dd(x'^alpha) arrow dd(x^alpha)$ leaves $dd(s^2)$ invariant. From the chain-rule we can write
$
  dd(x'^alpha) = pdv(x'^alpha, x^beta) dd(x^beta)
$
this immediately tells us how contravariant four-vectors transform
$
  x^alpha -> x'^alpha = pdv(x'^alpha, x^beta) x^beta
$
similarly one can look at
$
  pdv(, x'^alpha) = pdv(x^beta, x'^alpha) pdv(, x^beta)
$
to get how covariant four-vectors transform
$
  x_alpha -> x'_alpha = pdv(x^beta, x'^alpha) x_beta
$
as expected these are inverses
$
  pdv(x^beta, x'^alpha) pdv(x'^alpha, x^gamma) = tensor(delta, beta, -gamma)
$
now we can write how the metric should transform
$
  g'_(alpha beta) = pdv(x^gamma, x'^alpha) pdv(x^delta, x'^beta) g_(gamma delta)
$

We'd like to be able to take derivatives of tensors, but it's easy to see that $partial_alpha x^beta$ doesn't transform correctly, and is thus not a tensor. To solve this we differentiate the definition of $x^beta$ to get
$
  partial'_beta x'^alpha &= pdv(, x'^beta) (pdv(x'^alpha, x^gamma) x^gamma) \
  &= pdv(x^delta, x'^beta) partial_delta (pdv(x'^alpha, x^gamma) x^gamma) \
  &= pdv(x^delta, x'^beta) pdv(x'^alpha, x^gamma) partial_delta x^gamma + pdv(x^delta, x'^beta) x^gamma partial_delta pdv(x'^alpha, x^gamma) \
  &= pdv(x^delta, x'^beta) pdv(x'^alpha, x^gamma) partial_delta x^gamma + pdv(x^delta, x'^beta) pdv(x'^alpha, x^delta, x^gamma) x^gamma \
  &= pdv(x^delta, x'^beta) pdv(x'^alpha, x^gamma) partial_delta x^gamma + pdv(x'^alpha, x'^beta, x^gamma) x^gamma
$
since $ pdv(, x'^beta) = pdv(x^delta, x'^beta) pdv(, x^delta) $ the second term is what makes this a non-tensor. For this reason we seek a covariant derivative $D_alpha$, with the following property
$
  D_beta x^alpha -> D'_beta x'^alpha = pdv(x^gamma, x'^beta) pdv(x'^alpha, x^delta) D_gamma x^delta
$
consider the derivatives of the four-vector itself $partial_alpha bold(x)$, this transforms like
$
  partial_alpha bold(x) -> partial'_alpha bold(x) = pdv(x^gamma, x'^alpha) partial_gamma bold(x)
$
since the four-vector is coordinate independent. We can act with $bold(e)'^beta$---the inverse basis vectors---on both sides giving
$
  bold(e)'^beta dot partial'_alpha bold(x) = pdv(x^gamma, x'^alpha) pdv(x'^beta, x^delta) bold(e)^delta dot partial_gamma bold(x)
$
so we have just found $D_alpha x^beta = bold(e)^beta dot partial_alpha bold(x)$---making the covariant derivative the expansion coefficients of $partial_alpha bold(x)$. This gives us a nice relation because we can write
$
  partial_beta x^alpha = bold(e)^alpha dot (partial_beta bold(x)) + bold(x) dot (partial_beta bold(e)^alpha)
$
since $x^alpha = bold(e)^alpha dot bold(x)$. We now define $partial_beta bold(e)^alpha = - tensor(Gamma, alpha, -gamma beta) bold(e)^gamma$, or $bold(x) dot (partial_beta bold(e)^alpha) = - tensor(Gamma, alpha, -gamma beta) x^gamma$, giving
$
  D_beta x^alpha = partial_beta A^alpha + tensor(Gamma, alpha, -gamma beta) A^gamma
$
similarly
$
  D_beta x_alpha = partial_beta x_alpha - tensor(Gamma, gamma, -beta alpha) A_gamma
$
a mixed tensor $tensor(T, alpha, -beta) = A^alpha B_beta$ will have a covariant derivative
$
  D_gamma tensor(T, alpha, -beta) = partial_gamma tensor(T, alpha, -beta) - tensor(Gamma, delta, -gamma beta) tensor(T, alpha, -delta) + tensor(Gamma, alpha, -gamma epsilon) tensor(T, epsilon, -beta)
$
so each index gets its own Christoffel symbol, with a negative for covariant indices. A big example is $g_(mu nu)$,
$
  D_lambda g_(mu nu) = partial_lambda g_(mu nu) - tensor(Gamma, rho, -lambda mu) g_(rho nu) - tensor(Gamma, rho, -lambda nu) g_(mu rho)
$
in how we've defined things we have
$
  tensor(Gamma, mu, -nu lambda) = tensor(Gamma, mu, -lambda nu)
$
even though the metric is position dependent, it is covariantly constant $D_lambda g_(mu nu) = 0$. To see this consider
$
  partial_lambda (bold(e)_mu dot bold(e)_nu) &= (partial_lambda bold(e)_mu) dot bold(e)_nu + bold(e)_mu dot (partial_lambda bold(e)_nu) \
  &= tensor(Gamma, rho, -lambda mu) bold(e)_rho dot bold(e)_nu + tensor(Gamma, rho, -lambda nu) bold(e)_mu dot bold(e)_rho
$
so in terms of $g_(mu nu)$:
$
  partial_lambda g_(mu nu) - tensor(Gamma, rho, -lambda mu) g_(rho nu) - tensor(Gamma, rho, -lambda nu) g_(mu rho) = D_lambda g_(mu nu) = 0
$
this is an essential property. Now we derive a representation of the Christoffel symbols in terms of the metric. We start with three equations
$
  D_lambda g_(mu nu) &= partial_lambda g_(mu nu) - tensor(Gamma, rho, -lambda mu) g_(rho nu) - tensor(Gamma, rho, -lambda nu) g_(mu rho) = 0 \
  D_nu g_(lambda mu) &= partial_nu g_(lambda mu) - tensor(Gamma, rho, -nu lambda) g_(rho mu) - tensor(Gamma, rho, -nu mu) g_(lambda rho) = 0 \
  - D_mu g_(nu lambda) &= - partial_mu g_(nu lambda) + tensor(Gamma, rho, -mu nu) g_(rho lambda) + tensor(Gamma, rho, -mu lambda) g_(nu rho) = 0
$
so we rotate the indices cyclically. Summing these and using symmetry gives
$
  partial_lambda g_(mu nu) + partial_nu g_(lambda mu) - partial_mu g_(nu lambda) - 2 tensor(Gamma, rho, -lambda nu) g_(mu rho) = 0
$
or equivalently
$
  tensor(Gamma, lambda, -mu nu) = 1/2 g^(lambda rho) (partial_nu g_(mu rho) + partial_mu g_(nu rho) - partial_rho g_(mu nu))
$
this is the fundamental theorem of Riemannian geometry.

#pagebreak()
= The geodesic equation
== Relativity as geometry
According to Einstein's theory all gravitational effects are caused by the underlying curved spacetime---this is motivated by the equivalence principle. General relativity states that matter and energy warps spacetime $g_(mu nu) eq.not eta_(mu nu)$.

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

== The equation
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
  dv(x^nu, lambda, 2) + tensor(Gamma, nu, -sigma rho) dv(x^sigma, lambda) dv(x^rho, lambda) = 0
$

=== Newtonian limit
Naturally this should reduce to the Newtonian case in the Newtonian limit of a particle moving with $v << c$ in a static and weak gravitational field. Nonrelativistic speed $dd(x^i)\/dd(t) << c$ implies $dd(x^i) << c dd(t)$ so
$
  dv(x^i, lambda) << c dv(t, lambda) = dv(x^0, lambda)
$
keeping only the dominant term in the double sum then gives
$
  dv(x^mu, lambda, 2) + tensor(Gamma, mu, -00) dv(x^0, lambda) dv(x^0, lambda) = 0
$
for a static field $dd(g_(mu nu), d: partial) \/dd(x^0, d: partial) = 0$ since all time derivatives should vanish, then the Christoffel symbol becomes
$
  g_(mu nu) tensor(Gamma, mu, -0 0) = - 1/2 pdv(g_(0 0), x^nu)
$
for a weak field we assume the metric is close to $eta_(mu nu)$ meaning
$
  g_(mu nu) = eta_(mu nu) + h_(mu nu)
$
with $h_(mu nu) << 1$ being a small correction field. Since $eta_(mu nu)$ is a constant metric we have $dd(g_(mu nu), d: partial) \/ dd(x^sigma, d: partial) = dd(h_(mu nu), d: partial)\/dd(x^sigma, d: partial)$, so to leading order
$
  eta_(mu nu) tensor(Gamma, mu, -0 0) = -1/2 pdv(h_(0 0), x^nu)
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
= The covariance principle
The covariance principle states that; our equations must be covariant under the general coordinate transformations which leave $dd(s^2)$ invariant, and our equations should reduce to the correct special relativistic form in local inertial frames---additionally our equations reduce to Newtonian equations in the Newtonian limit. It can be treated as a generalization of the equivalence principle.

This gives us a clear path to convert our equations from special relativity to general relativity, since the tensor formalism developed for the two only differ in their derivatives. We simply replace $partial$ by $D$ to convert our equations,
$
  partial -> D " "(= partial + Gamma)
$
this is minimal---gravitational---coupling. Christoffel symbols are derivatives of the metric---our potentials---so they bring the gravitational field strength into our equations naturally.

As an example consider Maxwell's equations and the Lorentz force law from special relativity,
$
  dv(U^mu, tau) = q/c F^(mu nu) U_nu ->^"minimal coupling" (D U^mu)/(dd(tau)) = q/c F^(mu nu) U_nu
$
$
  partial_mu F^mu nu = - 1/c j^nu & ->^"m.c" D_mu F^(mu nu) =-1/c j^nu \
  partial_mu tilde(F)^(mu nu) = 0 & ->^"m.c" D_mu tilde(F)^(mu nu)=0
$
by taking $F^(mu nu) = 0$ in the Lorentz force law we can rederive the geodesic equation,
$
  (D U^mu)/dd(tau) = 0
$
just using the definition of $D$ and $U^mu = dd(x^mu)\/dd(tau)$ gives
$
  dv(x^mu, tau, 2) + tensor(Gamma, mu, -nu lambda) dv(x^nu, tau) dv(x^lambda, tau) = 0
$
which is just the geodesic equation.
