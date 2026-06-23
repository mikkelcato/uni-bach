#import "../../temp.typ": *
#show: chpt-note.with()

= Advanced Methods for PDEs
== Separation of variables
We assume coordinates $vecrow(x, y, z, t)$ and take our PDE to be defined on a set $X subset.eq RR^4$. We only consider homogeneous equations.

A solution is _separable_ if we can factorise it as
$
  u(x,y,z,t) = X(x) Y(y) Z(z) T(t)
$
We may refer to a solution as _partially separable_ with the meaning being obvious.

#example(name: "The heat equation")[We consider the equation
  $
    laplacian u (bold(r)) = 1/c^2 pdv(u(bold(r)), t, 2)
  $
  We make the ansatz $u = X Y Z T$ giving
  $
    X''/X + Y''/Y + Z''/Z = 1/c^2 T''/T
  $
  This can only be true for all ${x,y,z,t}$ if each term is constant. We end up with four separable ODEs
  $
    X''/X & = -l^2";  " Y''/Y = -m^2";  " Z''/Z = -n^2";  " T''/T = - c^2 mu^2
  $
  along with the requirement
  $
    - (l^2+m^2+n^2) = -mu^2
  $
  we refer to the constant as _separation constants_. All four ODEs have similar solutions with one being
  $
    X = A e^(i l x) + B e^(-i l x)
  $
  where $A,B,dots$ are specified using boundary conditions.
]

Generally a separable solution with $n$ independent variables has $n-1$ independent separation constants.

Assuming the PDE is linear then the _superpositions principle_ holds. Then we can form general solutions by constructing superpositions of solutions corresponding to different values of separation constants. Consider $X subset.eq RR^2$ with $u_(lambda_1) = X_(lambda_1) Y_(lambda_1)$ being a solution where the separation constant takes the value $lambda_1$. Then
$
  u = sum_i a_i X_(lambda_i) Y_(lambda_i)
$
is also a solution for any $a_i$ provided the $lambda_i$ work.

== Integral transforms
We show how the _Fourier transform_#footnote[Details in exercises.] can be used to solve PDEs by again considering the heat equation
$
  alpha pdv(u, x, 2) = pdv(u, t)",   for" t>0 "and" u(x,0) = f(x) "with" x in RR
$
To make life simple we assume $u -> 0$ and $partial u -> 0$ as $abs(x) -> oo$. Then we can take the Fourier transform of the PDE
$
  alpha/sqrt(2 pi) integral_(-oo)^oo pdv(u, x, 2) e^(-i k x) dd(x) = 1/sqrt(2 pi) pdv(, t) integral_(-oo)^oo u e^(-i k x) dd(x)
$
which by definition becomes
$
  - alpha k^2 tilde(u) (k,t) = pdv(tilde(u) (k,t), t)
$
with
$
  tilde(u) = 1/sqrt(2pi) integral_(-oo)^oo u(x,t) e^(-i k x) dd(x) tilde "Fourier transform"
$
The solution is
$
  tilde(u) (k,t) = tilde(u) (k,0) e^(-alpha k^2 t)
$
where
$
  tilde(u) (k,0) = 1/sqrt(2 pi) integral_(-oo)^oo u(x,0) e^(-i k x) dd(x) = tilde(f) (k)
$
We find
$
  tilde(u) (k,t) = tilde(f) (k) e^(-alpha k^2 t) = sqrt(2 pi) tilde(f) (k) tilde(G) (k,t)
$
where
$
  tilde(G) (k,t) = 1/sqrt(2 pi) e^(-alpha k^2 t)
$
We see $tilde(u)$ is a product of two Fourier transforms. We can then apply the _convolution theorem_ to find
$
  u(x,t) = integral_(-oo)^oo G(x-x',t) f(x') dd(x')
$
We need $G(x,t)$#footnote[This is the _Green's function_ as defined below.]
$
  G(x,t) &=^"inverse Fourier" 1/(2 pi) integral_(-oo)^oo e^(-alpha k^2 t) e^(i k x) dd(k) \
  &= 1/sqrt(4 pi alpha t) e^(-x^2/(4 alpha t))
$
We have found the solution
$
  u = 1/sqrt(4 pi alpha t) integral_(-oo)^oo exp(- (x-x')^2/(4 alpha t)) f(x') dd(x')
$
Which becomes the previously discussed solution for $f(x) = delta(x)$.

== Green's functions
Let $cal(P)$ denote a linear partial differential operator. Then a linear inhomogeneous PDE has the form
$
  cal(P) u(bold(r)) = rho(bold(r))
$
where $bold(r) in RR^n$.

=== ODEs and Green's functions
Let $x,b in CC^n$ and $A in "GL"(CC,n)$. We want to solve the linear system
$
  A x = b
$
By definition $A$ is _invertible_ so the solution is
$
  x = A^(-1) b equiv M b
$
with $A_(i j) M_(j k) = delta_(i k)$. We will see that finding a _Green's function_ amounts to being the continuous equivalent of the matrix inverse.

Before defining the Green's function we define the _Dirac $delta$-function_ by
$
  delta(x) = lim_(epsilon -> 0^+) delta_epsilon (x)
$
with
$
  delta_epsilon (x) = 1/(epsilon sqrt(pi)) e^(-x^2\/epsilon^2)
$
or more simply as
$
  delta(x) = cases(0 #h(2em)&"for" x eq.not 0, oo &"for" x = 0)
$
We can also think of it as the derivative of the _step function_. We care about the $delta$-function due to the property
$
  integral dd(x) f(x) delta(x - x_0) = f(x_0)
$
i.e. it picks out the value of $f$ at $x_0$.

Consider a linear ODE of order $n$
$
  cal(P) y(x) = f(x)
$
where
$
  cal(P) = a_n (x) dv(, x, n) + a_(n-1) (x) dv(, x, n-1) + dots + a_1 (x) dv(, x) + a_0 (x)
$
We would like to invert $cal(P)$ and write $y = cal(P)^(-1) f$. We denote this object as the _Green's function_. Generalising the discrete $x_i = M_(i j) b_j$ we have
$
  y(x) = integral_a^b G(x,z) f(z) dd(z)
$
Then
$
  cal(P) y(x) = integral_a^b cal(P) G(x,z) f(z) dd(z) =^! f(x)
$
implying
$
  cal(P) G(x,z) = delta(x-z)
$
with $a <= x <= b$. Then if we can find $G(x,z)$ we have found our solution!

#example[
  Consider the second order ODE
  $
    dv(y, x, 2) = delta(x-1)
  $
  with
  $
    evaluated(dv(y, x))_(x=0) = 0";  " y(0) = 0
  $
  We consider two regions $x < 1$ and $x > 1$.

  For $x < 1$ we have
  $
    dv(y, x, 2) = 0 => y = A x + B
  $
  with the boundary conditions we immediately have $A = B = 0$.

  For $x > 1$ we have
  $
    y = C x + D
  $
  However, now the boundary conditions are of no use. To proceed we show $y$ is continuous at $x = 1$. By contradiction if it were discontinuous then $dv(y, x) tilde delta(x-1)$ and
  $
    dv(y, x, 2) tilde dv(, x) delta(x-1)
  $
  which is a contradiction. We can also show that $dv(y, x)$ is discontinuous as $x = 1$ since
  $
    dv(y, x) = Theta(x-1) + "const" tilde "discontinuous"
  $
  By continuity we have
  $
    0 = evaluated(C x + D)_(x=1) => C = -D
  $
  By discontinuity we have
  $
    underbracket(lim_(x->1^-) dv(y, x), 0) + 1 = underbracket(lim_(x->1^+) dv(y, x), C)
  $
  implying $C = 1$ and $D = -1$. Then
  $
    y(x) & = cases(x-1 #h(1em)&"for" x>=1, 0 &"for" x <= 1) \
         & = (x-1) Theta(x-1)
  $
]

We now consider more general properties of Green's functions for ODEs. We have
$
  sum_(m=0)^n a_m dv(, x, m) G(x,z) = 0";  for" x<z "and" x>z
$
Then
$
  lim_(epsilon->0) integral_(z-epsilon)^(z+epsilon) cal(P) G(x,z) dd(x) &= lim_(epsilon->0) sum_(m=0)^n integral_(z-epsilon)^(z+epsilon) a_m (x) dv(, x, m) G(x,z) dd(x) \
  &= lim_(epsilon->0) integral_(z-epsilon)^(z+epsilon) delta(x-z) dd(x) \
  &= 1
$
We can simplify the LHS since all
$
  dv(, x, m) G(x,z) "for" 0<=m< n-1 tilde "continuous"
$
by an argument similar to that in the previous example. These terms then vanish in the limit $epsilon -> 0$. The $(n-1)$th derivative has a finite discontinuity at $x = z$ and also vanishes. We are left with
$
  1 &= lim_(epsilon->0) integral_(z-epsilon)^(z+epsilon) a_n (x) dv(, x, n) G(x,z) dd(x) \
  &=^"ibp" lim_(epsilon->0) [a_n (x) dv(, x, n-1) G(x,z)]_(z-epsilon)^(z+epsilon)
$
implying that $dv(, x, n-1) G(x,z)$ has a discontinuity of size $a_n^(-1)$ at $x = z$. We have found $n$ constraints being that $G$ and all derivative up to order $n-2$ are continuous along with the above finite discontinuity.

=== PDEs and Green's functions
We consider the _Poisson equation_
$
  laplacian u(bold(r)) = rho(bold(r))
$
and will use Green's functions to solve boundary value problems in a given volume $V subset.eq RR^n$ with inhomogeneous boundary conditions.#footnote[Dirichlet or Neumann. With both given the problem is _overdetermined_.]

We recall _Green's second identity_. Let $phi.alt (bold(r))$ and $psi(bold(r))$ be scalar functions defined in some volume $V subset.eq RR^n$ with boundary $cal(S) = partial V$. Then#footnote[This follows from the _divergence theorem_ applied to $bold(F) = phi.alt grad psi - psi grad phi.alt$.]
$
  integral_V (phi.alt laplacian psi - psi laplacian phi.alt) dd(V) = integral_cal(S) (phi.alt grad psi - psi grad phi.alt) dot hat(n) dd(cal(S))
$
We need the $G$ satisfying
$
  laplacian G(bold(r),bold(r)_0) = delta(bold(r)-bold(r)_0)
$
with $bold(r)_0 in V$. We define $phi.alt (bold(r)) = u(bold(r))$ and $psi(bold(r)) = G(bold(r),bold(r)_0)$ and obtain
$
  integral_V [u(bold(r)) underbracket(laplacian G(bold(r),bold(r)_0), delta(bold(r)-bold(r)_0)) - G(bold(r),bold(r)_0) underbracket(laplacian u(bold(r)), rho(bold(r)))] dd(V_bold(r)) = integral_cal(S) [u(bold(r)) pdv(G(bold(r),bold(r)_0), n) - G(bold(r),bold(r)_0) pdv(u(bold(r)), n)] dd(cal(S))_bold(r)
$
where we integrate over $bold(r)$. We find
$
  u(bold(r)_0) = integral_V G(bold(r),bold(r)_0) rho(bold(r)) dd(V)_bold(r) + integral_(cal(S)) [u(bold(r)) pdv(G(bold(r),bold(r)_0), n) - G(bold(r),bold(r)_0) pdv(u(bold(r)), n)] dd(cal(S))_bold(r)
$
which is our solution.

We assume Dirichlet boundary conditions meaning $u(bold(r))$ is given on $cal(S)$. Let $f(bold(r)) : RR^n -> RR$ be a function with $u(bold(r)) = f(bold(r))$ for $bold(r) in cal(S)$. We are free to choose the boundary condition satisfied by $G(bold(r),bold(r)_0)$ so we pick $G(bold(r),bold(r)_0) = 0$ for $bold(r) in cal(S)$. This defines the _Dirichlet Green's function_. Then the third term on the RHS vanishes
$
  u(bold(r)_0) = integral_V G(bold(r),bold(r)_0) rho(bold(r)) dd(V)_bold(r) + integral_(cal(S)) f(bold(r)) pdv(G(bold(r),bold(r)_0), n) dd(cal(S))_bold(r)
$
We now need to find $G$ satisfying $G(bold(r),bold(r)_0) = 0$ for $bold(r) in cal(S)$. We make the ansatz
$
  G(bold(r),bold(r)_0) = F(bold(r),bold(r)_0) + H(bold(r),bold(r)_0)
$
and require
$
  laplacian F = delta(bold(r)-bold(r)_0)";  " laplacian H = 0
$
within the volume $V$ and $F + H = 0$ on $cal(S)$. We refer to $F$ as the _fundamental solution_.

#example[
  We consider the Poisson equation in three dimensions. We assume $F -> 0$ as $abs(bold(r)) -> oo$. We consider a sphere $cal(S)$ of radius $R$ centered at $bold(r)_0$ and integrate over the enclosed volume $V$
  $
    integral_V laplacian F dd(V)_bold(r) = underbracket(integral_V delta(bold(r)-bold(r)_0) dd(V)_bold(r), 1)
  $
  by Gauss'
  $
    integral_V laplacian F dd(V)_bold(r) = integral_cal(S) grad F dot hat(n) dd(cal(S))
  $
  and by symmetry $F(bold(r),bold(r)_0) eq F(abs(bold(r)-bold(r)_0)) eq F(r)$ implying $F$ has the same value everywhere on $cal(S)$. Writing the RHS in spherical coordinates we find#footnote[Using $dd(cal(S)) = r^2 sin theta dd(theta, phi.alt)$ and $grad f = pdv(f, r) hat(r)$.]
  $
    dv(F, r) = 1/(4 pi r^2) => F(r) = - 1/(4 pi r) + K
  $
  since $F -> 0$ as $r -> oo$ we have $K = 0$ and we obtain
  $
    F(bold(r),bold(r)_0) = - 1/(4 pi abs(bold(r)-bold(r)_0))
  $
]

#example[
  We can now solve
  $
    laplacian u(bold(r)) = - rho(bold(r))/epsilon_0 tilde "Gauss' law"
  $
  in $V = RR^3$ with $u(bold(r)) -> 0$ as $r -> oo$.

  The $F$ found above fully specifies $G$ since it already satisfies the boundary condition.#footnote[We set $H = 0$.] The surface integral in the solution we found vanishes aswell since $f(bold(r)) = 0$. The solution is then
  $
    u(bold(r)_0) = integral_V (rho(bold(r)) dd(V)_bold(r))/(4 pi epsilon_0 abs(bold(r)-bold(r)_0))
  $
  which should be familiar! This example was made trivial since $V = RR^3$ meaning we did not need to find $H$.
]

We now consider more general $V subset RR^3$ where $F eq.not 0$ on $cal(S)$. We will construct the corresponding $H$ by the _method of images_. The idea is to add solutions $H$ of the homogeneous equation to $F$ as to cancel the contributions of $F$ on $cal(S)$ meaning $G(bold(r),bold(r)_0) = 0$ on $cal(S)$. We do this by adding _copies_ of $F$ outside the volume $V$#footnote[We call these _image sources_.] to $F$ itself.

We have a single source $delta(bold(r)-bold(r)_0)$ inside the volume $V$ and add $N$ image sources outside the volume $V$
$
  sum_(l)^N underbracket(q_l, "strength") delta(bold(r)-bold(r)_l)
$
with $bold(r)_l$ outside the volume $V$. With an image source outside the volume $V$ its corresponding fundamental solution satisfies the homogeneous equation inside the volume $V$.#footnote[Since $delta(bold(r)-bold(r)_l) = 0$ everywhere inside the volume $V$.] We can write
$
  G(bold(r),bold(r)_0) &= F(bold(r),bold(r)_0) + underbracket(sum_l^N q_l F(bold(r),bold(r)_l), H(bold(r),bold(r)_0)) \
  &= "real" + "images"
$
We can then adjust $bold(r)_l$ and $q_l$ until $G(bold(r),bold(r)_0) = 0$ on $cal(S)$. With $G(bold(r),bold(r)_0)$ we have the solution to the Dirichlet boundary value problem.

#example[
  We consider
  $
    laplacian u (bold(r)) = 0 "for" z>0
  $
  with $u(bold(r)) = f(x,y)$ for $z = 0$. The volume is $V = {(x,y,z) in RR^3: z > 0}$ and the boundary is $partial V = RR^2$ alongside the surface at infinity. The Green's function must satisfy $G(bold(r),bold(r)_0) = 0$ in the plane $z=0$ and $G(bold(r),bold(r)_0) -> 0$ as $abs(bold(r))->oo$. This is achieved by an image source at $bold(r)_1$ being the reflection of $bold(r)_0$ about $z = 0$ with strength $q = -1$
  $
    G(bold(r),bold(r)_0) = - 1/(4 pi abs(bold(r)-bold(r)_0)) + 1/(4 pi abs(bold(r)-bold(r)_1))
  $
  where $bold(r)_0 = vecrow(x_0, y_0, z_0)$ and $bold(r)_1 = vecrow(x_0, y_0, -z_0)$.

  With $rho(bold(r)) = 0$ we find
  $
    u(bold(r)_0) = integral_cal(S) f(bold(r)) pdv(G(bold(r),bold(r)_0), n) dd(cal(S))_bold(r)
  $
  where the normal derivative in the $x y$-plane is
  $
    pdv(G, n) = - grad G dot hat(z) = - pdv(G, z)
  $
]

We now assume Neumann boundary conditions#footnote[This uniquely determines a solution up to an additive constant as shown.]
$
  f(bold(r)) = pdv(u(bold(r)), n) "on" cal(S)
$
These satisfy a _self-consistency_ relation
$
  integral_cal(S) f(bold(r)) dd(cal(S)) = integral_cal(S) grad u (bold(r)) dot hat(n) dd(cal(S)) =^"Gauss" integral_V laplacian u(bold(r)) dd(V) =^"Poisson" integral_V rho(bold(r)) dd(V)
$
Analogously to before we would be tempted to pick
$
  pdv(G(bold(r),bold(r)_0), n) = 0
$
This is not allowed due to the self-consistency relation
$
  integral_cal(S) pdv(G, n) dd(cal(S)) = integral_cal(S) grad G dot hat(n) dd(cal(S)) = integral_V laplacian G dd(V) = 1
$
The simplest choice is instead
$
  pdv(G(bold(r),bold(r)_0), n) = 1/A_cal(S)
$
with $bold(r)$ on $cal(S)$ and $A_cal(S)$ being the area of $cal(S)$. This defines the _Neumann Green's function_. Then we find
$
  u(bold(r)_0) = integral_V G(bold(r),bold(r)_0) rho(bold(r)) dd(V) + underbracket(1/A integral_cal(S) u(bold(r)) dd(cal(S)), expval(u)_cal(S) " " tilde "constant") - integral_cal(S) G(bold(r),bold(r)_0) f(bold(r)) dd(cal(S))
$
Where $expval(u)_cal(S)$ is a _freely specifiable_ constant. Assuming $cal(S)$ is at infinity and $u(bold(r)) -> 0$ as $abs(bold(r)) -> oo$ we force $expval(u)_cal(S) = 0$. The construction of $G(bold(r),bold(r)_0)$ proceeds as before by the method of images.#footnote[Assuming $G(bold(r),bold(r)_0) = F(bold(r),bold(r)_0) + H(bold(r),bold(r)_0)$ is motivated from $laplacian G(bold(r),bold(r)_0) = delta(bold(r)-bold(r)_0)$ which defines $G(bold(r),bold(r)_0)$. ] The _tuning_ of $q_l$ and $bold(r)_l$ is now done to satisfy the _area condition_.
