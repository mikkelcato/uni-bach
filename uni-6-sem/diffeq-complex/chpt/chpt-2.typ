#import "chpt-temp.typ": *
#show: chpt-note.with()

= Partial Differential Equations
Consider a higher dimensional space $X$ with $X subset.eq RR^n$. A _partial differential equation_ relates a function $u : X -> RR$ or $CC$ which depends on $n$ independent variables $x_1, dots, x_n$ to its partial derivatives. The computation of $u$ typically involves a boundary value problem leading to a dependence on the geometry of $X$ as well as the values of $u$ and its derivatives on the boundary $dd(X, d: partial)$.

#example[
  The following are all examples of _homogeneous_ and _linear_ partial differential equations we find in physics:
  $
      "The transport equation:" #h(1em) a pdv(u, x) + b pdv(u, y) & = 0 \
      "The Laplace equation:" #h(1em) pdv(u, x, 2) + pdv(u, y, 2) & = 0 \
    "The Schrödinger equation:" #h(1em) pdv(u, t) +i pdv(u, x, 2) & = 0 \
          "The heat equation:" #h(1em) pdv(u, t) - a pdv(u, x, 2) & = 0 \
         "The wave equation:" #h(1em) pdv(u, t, 2) - pdv(u, x, 2) & = 0
  $
]

Unless specified we assume $n = 2$.

== The heat equation
We briefly discuss the heat equation before discussing solution methods. We write the three-dimensional heat equation as

$
  pdv(u, t) = D dd(u, d: laplace)
$
This equation describes the evolution of some temperature distribution $u equiv u(bold(r),t)$. We call $D$ the diffusion constant and have defined the Laplacian $laplace equiv laplacian$.

#proof[
  We consider heat diffusing through a material with conductivity $k$, specific heat capacity $s$, and density $rho$. Let the heat inside a volume $V$ be $Q(t)$. The heat flux through $V$ is then
  $
    - dv(Q, t) = integral_(dd(V, d: partial)) (-k grad u) dot bold(n) dd(A)
  $
  which follows by definition. By the divergence theorem we can rewrite the RHS
  $
    - dv(Q, t) = integral_V div (-k grad u) dd(V)
  $
  We can also write $Q$ as
  $
    Q = integral s rho u dd(V)
  $
  implying
  $
    dv(Q, t) = integral s rho pdv(u, t) dd(V)
  $
  Comparing the above immediately gives the heat equation.
]

To make life simple we consider the one-dimensional heat equation
$
  pdv(u, t) - pdv(u, x, 2) = 0
$
We immediately see some trivial solutions
$
  u(x,t) & = x";  " u(x,t) = alpha x^2 + beta x + gamma + 2 alpha t
$
For non-trivial solutions we need $t > 0$ so $u$ is defined on $RR times (0,oo)$. An important such solution is
$
  u(x,t) = 1/sqrt(4 pi t) e^(-x^2\/4t)
$
At a fixed time $t_0$ we define $alpha_0 equiv sqrt(2 t_0)$. The solution then takes the form
$
  u(x,t_0) = 1/(sqrt(2 pi) alpha_0) e^(-x^2\/2alpha_0^2)
$
which is Gaussian! And since $alpha(t) = sqrt(2 t)$ the width of this Gaussian increases with time while the magnitude decreases.

== First order PDEs
=== Linear homogeneous PDEs
We define the first-order linear _partial differential operator_ $cal(P)$ by
$
  cal(P) & = A(x,y) pdv(, x) + B(x,y) pdv(, y) + C(x,y) \
  cal(P) & : C^1 (X) -> C^0 (X)
$
with $C^1 (X)$ being the space of differential functions on $X subset RR^2$. This guy acts on functions $u in C^1 (X)$ as
$
  u(x,y) |-> (cal(P) u) (x,y) = A(x,y) pdv(u, x) + B(x,y) pdv(u, y) + C(x,y) u
$
We say a partial differential equation is homogeneous if
$
  cal(P) u = 0
$
and it is linear if $cal(P)$ is linear
$
  cal(P) (alpha u + beta v) = alpha cal(P) u + beta cal(P) v
$
Then by definition to solve a homogeneous partial differential equation we need to determine the kernel of $cal(P)$.

=== Using equipotential lines
We assume $C = 0$. We will attempt to solve the problem by rewriting it in terms of an ordinary differential equation. When $n = 2$ we try to find equipotential lines defined by $p = p(x,y)$ along which $u(x,y)$ is constant meaning we can write $u(x,y) = f(p)$.

With this assumption we have
$
  pdv(u, x) = dv(f, p) pdv(p, x)";  " pdv(u, y) = dv(f, p) pdv(p, y)
$
implying
$
  underbracket((A pdv(p, x) + B pdv(p, y)), =^! 0) dv(f, p) = 0
$
By definition of $p$ we have
$
  dd(p) = pdv(p, x) dd(x) + pdv(p, y) dd(y) = 0
$
implying
$
  dd(x)/A = dd(y)/B tilde "parametrisation of lines"
$

#example[
  Consider
  $
    x pdv(u, x) - 2 y pdv(u, y) = 0
  $
  Then
  $
    dd(x)/x = - dd(y)/(2y)
  $
  implies
  $
    x^2 y = "constant" =^! p
  $
  Then the solution is
  $
    u(x,y) = f(x^2 y)
  $
  with $f$ being arbitrary.

  Additional boundary conditions restrict our solution. Assuming $ u(x=1,y) = 2 y + 1 $ we immediately find the particular solution
  $
    u(x,y) = 2 x^2 y + 1
  $
  Assuming
  $
    u(x=1, y=1) = 4
  $
  we find many more solutions
  $
    u(x,y) & = 4 x^2 y \
           & = 3 x^2 y + 1 \
           & = 4
  $
  This happens since we are much less restrictive with this particular boundary condition.
]

When $C eq.not 0$ the above does not work. We make the ansatz $u(x,y) = h(x,y) f(p)$ with $f(p)$ as above and $h(x,y)$ being arbitrary.

#example[
  Consider
  $
    x pdv(u, x) + 2 pdv(u, y) - 2 u = 0
  $
  Using the ansatz we find
  $
    (x pdv(h, x) + 2 pdv(h, y) - 2 h) f + (x pdv(p, x) + 2 pdv(p, y)) h dv(f, p) = 0
  $
  We assume we know some particular solution $h(x,y)$ to the original equation. Then the first term vanishes and we have
  $
    x pdv(p, x) + 2 pdv(p, y) = 0
  $
  which we can solve as before! We find
  $
    p = x e^(-y\/2)
  $
  and the general solution becomes
  $
    u(x,y) = h(x,y) f(x e^(-y\/2))
  $
  where $h(x,y)$ is any particular solution. An example in the above could be $h(x,y) = x^2$ which is easily seen by inspection.
]
=== Linear inhomogeneous PDEs
A partial differential equation is homogenoeus if given a solution $u(x,y)$ then $lambda u(x,y)$ with $lambda in CC$ is also a solution. Clearly this is equivalent to the linear case where we had
$
  cal(P) u = 0
$
A partial differential boundary value problem is homogeneous if the partial differential equation is homogeneous and given that if $u(x,y)$ satisfies the boundary conditions then $lambda u(x,y)$ satisfies the boundary conditions.

#example[
  An example of an inhomogeneous boundary value problem is
  $
    u(x,y) = 2 x^2 y + 1
  $
  with $u(x=1,y) = 2 y + 1$ which we derived above.
]

With $v in C^0 (X)$ with $X subset RR^2$ then an inhomogeneous partial differential equation can be written as
$
  cal(P) u = v
$
with $u in C^1 (X)$. This can also be written as
$
  A(x,y) pdv(u, x) + B(x,y) pdv(u, y) + C(x,y) u + R(x,y) = 0
$
The solution of these equations are related to their homogeneous counterpart since the general solution can be written as
$
  "particular solution to inhomogeneous" + "general solution to homogeneous"
$

#example[
  Consider
  $
    y pdv(u, x) - x pdv(u, y) = 3 x
  $
  with $u(x,y=0) = x^2$. We can find the solution to the homogeneous equation by equipotential lines
  $
    dd(x)/y = - dd(y)/x => p = x^2 + y^2
  $
  Then $u_"homogeneous" (x,y) = f(x^2+y^2)$. We now need some particular solution to the inhomogeneous equation. By inspectiong it is easy to find $u_"particular" (x,y) = -3 y$. The general solution is then
  $
    u(x,y) = f(x^2+y^2) - 3 y
  $
  Applying the boundary condition $u(x,y=0) = f(x^2) =^! x^2$ meaning $f(p) = p$ so
  $
    u(x,y) = x^2+y^2 - 3 y
  $
]

=== Characteristics for first order PDEs
Consider the general first-order partial differential equation
$
  A(x,y) pdv(U, x) + B(x,y) pdv(u, y) = F(x,y,u)
$
with the boundary condition $u(x,y) = phi.alt(s)$ which is defined along some curve $Gamma subset RR^2$. Let the curve be parametrised by $(x,y) = (x(s),y(s))$. Then
$
  dv(u, s) = pdv(u, x) dv(x, s) + pdv(u, y) dv(y, s)
$
These can be written as a system of equations in the following way
$
  underbracket(mat(dv(x, s), dv(y, s); A(x,y), B(x,y)), M) vec(pdv(u, x), pdv(u, y)) = vec(dv(phi.alt, s), F(x,y,u))
$
Given $det M eq.not 0$ this can be solved to determine $grad u(x,y)$.

The condition $det M = 0$ determines a set of curves we call the _characteristics_. This determinant leads to
$
  dv(y, x) = (B(x,y))/(A(x,y))
$
which the curves must satisfy. This condition is the same which determined the equipotential lines.

When $vecrow(dv(x, s), dv(y, s))$ depends linearly on $vecrow(A, B)$ we see that $det M = 0$. This is equivalent to saying that the curve on which the boundary condition is specified coincides with the characteristics. Since these correspond to equipotential lines the function $u(x,y)$ is constant along the characteristics. Given this is not the case we can solve to find $grad u$ and given $u(x,y) = phi.alt (s)$ on $Gamma$ we can specify the solution locally. This is just done by a Taylor expansion.

#example[
  Consider
  $
    x pdv(u, x) - 2 y pdv(u, y) = 0
  $
  with $u(x,y) = 2 y + 1$ on the line $x = 1$, $y in[0,1]$. We have already found the general solution to be
  $
    u(x,y) = f(x^2 y)
  $
  with the characteristics $c = x^2 y$. The particular solution obeying the boundary condition is
  $
    u(x,y) = 2 x^2 y + 1
  $
  with $y in [0,1]$.

  This solution is determined everywhere along a given characteristic since $f(x^2 y)$ is constant. The more general solution is
  $
    u(x,y) = 2 x^2 y + 1 + g(x^2 y)
  $
  with $g(z) = 0$ for $0 <= z <= 1$.
]

Generally if $Gamma$ is itself a characteristic curve then the solution is unspecified everywhere expect on $Gamma$. While if $Gamma$ crosses some characteristic multiple times then this can lead to overdetermination meaning no solution may exist. The line $Gamma$ should only cross each characteristic once and then the solution is specified along these.

== Second order PDEs
The general second order partial differential equation is of the form
$
  A(x,y) pdv(u, x, 2) + B(x,y) pdv(u, y, x) + C(x,y) pdv(u, y, 2) + D(x,y) pdv(u, x) + E(x,y) pdv(u, y) + F(x,y) u = R(x,y)
$
We classify these equations in the following way:
$
  "hyperbolic" & <--> B^2 > 4 A C \
   "parabolic" & <--> B^2 = 4 A C \
    "elliptic" & <--> B^2 < 4 A C
$
This classification is local since the behaviour might change within $X subset.eq RR^2$.

=== Homogeneous second order PDEs with constant coefficients
We assume constant coefficients with $D = E = F = R= 0$. Then
$
  A pdv(u, x, 2) + B pdv(u, y, x) + C pdv(u, y, 2) = 0
$
We make the ansatz $u(x,y) = f(p)$ and obtain
$
  A [dv(f, p, 2) (pdv(p, x))^2 + dv(f, p) pdv(p, x, 2)] + B [dv(f, p, 2) pdv(p, y) pdv(p, x)+ dv(f, p) pdv(p, y, x) ] + C [dv(f, p, 2) (pdv(p, y))^2 + dv(f, p) pdv(p, y, 2)] = 0
$
We assume $p$ is linear in $x$ and $y$. Then
$
  p = a x + b y
$
and we find
$
  underbracket((A a^2 + B a b + C b^2), =^! 0) dv(f, p, 2) = 0
$
This is a quadratic
$
  lambda_(1,2) = b/a = (-B plus.minus sqrt(B^2-4 A C))/(2 C)
$
Then any functions of
$
  p_1 = x + lambda_1 y";  " p_2 = x + lambda_2 y
$
satisfy the original equation. Then the general solution can be written as
$
  u(x,y) = f(x+lambda_1 y) + g(x+lambda_2 y)
$
#example[
  Consider
  $
    pdv(u, x, 2) + pdv(u, y, 2) = 0
  $
  Then $A = C = 1$ and $B = 0$ meaning the equation is elliptic. We find $lambda_(1,2) = plus.minus i$.
]

#example[
  Consider $A, C in RR$ and $B = 0$. The solutions have arguments
  $
    "hyperbolic" & -> x plus.minus alpha y \
      "elliptic" & -> x plus.minus i beta y
  $
  with $alpha, beta in RR$.
]

#example[
  Consider $B^2 = 4 A C$. Then the solution is two-fold degenerate with
  $
    lambda = - B/(2 C)
  $
  yielding the solution
  $
    u(x,y) = f(x - B/(2 C) y)
  $
  which is not completely general!
]

We would like to generalise our solution. We make the ansatz
$
  u(x,y) = h(x,y) g(x- B/(2 C) y)
$
and eventually obtain
$
  (A pdv(h, x, 2) + B pdv(h, y, x) + C pdv(h, y, 2)) g = 0
$
Then $h$ must be any solution to the original partial differential equation. The simplest such solution is just $h(x,y) = x$ giving the general solution for the parabolic case
$
  u(x,y) = f(x- B/(2 C) y) + x g (x- B/(2 C) y)
$
which is quite nice.

A complementary derivation relies on using _variable transformations_ since these can trivialise the partial differential equation. We start with the hyperbolic and elliptic cases where $B^2 eq.not 4 A C$. We define
$
  xi equiv x + lambda_1 y";  " eta equiv x + lambda_2 y
$
implying
$
  pdv(u, x) = (pdv(, xi) + pdv(, eta)) u";  "pdv(u, y) = (lambda_1 pdv(, xi) + lambda_2 pdv(, eta)) u
$
Then the partial differential equation becomes
$
  [A (pdv(, xi) + pdv(, eta))^2 + B (pdv(, xi) + pdv(, eta)) (lambda_1 pdv(, xi) + lambda_2 pdv(, eta)) + C (lambda_1 pdv(, xi) + lambda_2 pdv(, eta))^2] u(xi,eta) = 0
$
We define $lambda_(1,2)$ by
$
  A + B lambda_i + C lambda_i^2 = 0
$
meaning they are the same as before. Then the partial differential equation becomes
$
  [2 A + B (lambda_1 + lambda_2) + 2 C lambda_1 lambda_2] pdv(u, xi, eta) = 0
$
Since $B^2 eq.not 4 A C$ the term in the brackets never vanishes so
$
  pdv(u, xi, eta) =^! 0
$
This has the solution
$
  pdv(u, eta) = F(eta)";  " pdv(u, xi) = G(xi)
$
implying
$
  u(xi,eta) & = f(xi) + g(eta) \
            & = f(x+lambda_1 y) + g(x+lambda_2 y)
$
as we found before.

For the parabolic case we assume $B^2 = 4 A C$ and define
$
  lambda equiv - B/(2 C)
$
Then we do
$
  xi = x+ lambda y";  " eta = x
$
Which gives
$
  A pdv(u, eta, 2) = 0
$
implying
$
  u(x,y) = x g(x+lambda y) + f(x+lambda y)
$
These transformations are typically not easy to find.

=== Homogeneous mixed first and second order PDEs with constant coefficients
We consider partial differential equations of the form
$
  alpha pdv(u, x, 2) = pdv(u, t)
$
with $u equiv u(x,t)$ and $alpha$ constant. This is simply the heat equation.

The previous strategy of $u = f(p)$ will not work since we will not be able to factor $f(p)$. We could instead try separation of variables which we discuss later or try setting both sides equal to a constant. We instead proceed by dimensional analysis and define
$
  eta equiv x^2/(alpha t)
$
which is dimensionless. We then try $u = f(eta)$. Then
$
     pdv(u, x) & = (2 x)/(alpha t) dv(f, eta) \
  pdv(u, x, 2) & = 2/(alpha t) dv(f, eta) + ((2 x)/(alpha t))^2 dv(f, eta, 2) \
     pdv(u, t) & = - x^2/(alpha t^2) dv(f, eta)
$
meaning the partial differential equation becomes
$
  4 eta dv(f, eta, 2) + (2+ eta) dv(f, eta) = 0
$
which is an ordinary differential equation! We can solve this to find $dv(f, eta)$ by separation of variables which can then be integrated to give
$
  f(eta) & prop integral_(eta_0)^eta mu^(-1\/2) e^(-mu\/4) dd(mu) \
         & prop integral_(zeta_0)^(zeta = x\/2 sqrt(alpha t)) e^(-nu^2) dd(nu)
$
where we do $zeta equiv eta^(1\/2)\/2$. Taking $zeta_0 = 0$ this is the _error function_
$
  u(x,t) prop "erf"(x/(2 sqrt(alpha t)))
$
note that for $zeta >= zeta_0$ we require $x >= 0$ and $t > 0$.

=== Characteristics for second order PDEs
We consider a general linear second order partial differential equation of the form
$
  A(x,y) pdv(u, x, 2) + B(x,y) pdv(u, x, y) + C(x,y) pdv(u, y, 2) = F(x,y,u,pdv(u, x),pdv(u, y))
$
Generally the boundary conditions should be such that $u$ and its first partial derivatives are specified along a suitable set of boundaries bordering or enclosing the region of interest.

We classify boundary conditions as follows:

1. _Dirichlet_ boundary conditions specify the value of $u$ at each point of the boundary.

2. _Neumann_ boundary conditions specify the value of $dd(u, d: partial)\/dd(n, d: partial)$ being the _normal derivative_ at each point of the boundary. We define $ pdv(u, n) equiv grad u dot hat(n) $ with $hat(n)$ being normal to the boundary at each point.

3. _Cauchy_ boundary conditions specify $u$ and $dd(u, d: partial)\/dd(n, d: partial)$ at each point of the boundary.

We consider the third case in depth. Assume there is some curve $Gamma subset.eq RR^2$ on which the boundary conditions are defined. We parametrise this cuve as $(x,y) = (x(s),y(s))$ where $s$ is the _curve parameter_. We denote the boundary conditions along $Gamma$ by $u(x,y) = phi.alt (s)$ and $dd(u, d: partial)\/dd(n, d: partial) = psi (s)$. At each point of $Gamma$ we have
$
   dd(bold(r)) & = dd(x) hat(x) + dd(y) hat(y) tilde "tangent" \
  hat(n) dd(s) & = dd(y) hat(x) - dd(x) hat(y) tilde "normal"
$
Then along $Gamma$ we have
$
  pdv(u, s) &equiv grad u dot dv(bold(r), s) = pdv(u, x) dv(x, s) + pdv(u, y) dv(y, s) = dv(phi.alt (s), s) \
  pdv(u, n) &equiv grad u dot hat(n) = pdv(u, x) dv(y, s) - pdv(u, y) dv(x, s) = psi(s)
$
We can write this as
$
  mat(dv(x, s), dv(y, s); dv(y, s), - dv(x, s)) vec(pdv(u, x), pdv(u, y)) = vec(dv(phi.alt, s), psi)
$
These can be solved to find $pdv(u, x)$ and $pdv(u, y)$ along $Gamma$ which provide information about the local neighbourhood of $Gamma$ and can be used to determine $u$.

Consider
$
  dv(, s) = dv(x, s) pdv(, x) + dv(y, s) pdv(, y)
$
giving
$
  dv(, s) pdv(u, x) & = dv(x, s) pdv(u, x, 2) + dv(y, s) pdv(u, x, y) \
  dv(, s) pdv(u, y) & = dv(x, s) pdv(u, y, x) + dv(y, s) pdv(u, y, 2)
$
Then we have
$
  mat(A, B, C; dv(x, s), dv(y, s), 0; 0, dv(x, s), dv(y, s)) vec(pdv(u, x, 2), pdv(u, x, y), pdv(u, y, 2)) = vec(F, dv(, s) pdv(u, x), dv(, s) pdv(u, y))
$
The characteristics are again determined by
$
  matrixdet(A, B, C; dv(x, s), dv(y, s), 0; 0, dv(x, s), dv(y, s)) = 0 => A (dv(y, x))^2 - B dv(y, x) + C = 0
$
This is an ordinary differential equation for the curves in the $x y$-plane along which the second partial derivatives of $u$ cannot be found. These curves have
$
  dv(y, x) = (B plus.minus sqrt(B^2-4 A C))/(2 A)
$
assuming $A eq.not 0$. These can be classified as before $tilde$ hyperbolic, parabolic, elliptic.

Assuming $A, B, C$ are constants then the characteristics are linear $x + lambda y = "const"$.

#example[
  Consider
  $
    pdv(u, x, 2) - 1/c^2 pdv(u, t, 2) = 0
  $
  We have $A = 1$, $C = -1\/c^2$ and $B = 0$ meaning the equation is hyperbolic. Then
  $
    dv(x, t) = plus.minus c
  $
  We find two families of characteristics
  $
    x minus.plus c t = "const"
  $
  The general solution is then given by
  $
    u(x,y) = underbracket(f(x-c t), "forward") + underbracket(g(x+ c t), "backward")
  $
  We specify Cauchy boundary conditions $u$ and $dd(u, d: partial)\/dd(n, d: partial)$ along a segment on the $x$-axis i.e. $t= 0$ and $x in [0,L]$ with $L > 0$. Then the solution is specified on every straight line
  in the family corresponding to a slope $+c$ and on every straight line in the family corresponding to a slope $-c$ which intersect on the line segment. Along the forward characteristics $p = x - c t tilde "constant"$ so $f(x- c t)$ is uniquely defined for the characteristics that intersect the line segment. Likewise for the backward characteristics $p = x + c t$. The solution is then _partially defined_ when either $f(x- c t)$ or $g(x + c t)$ are fully known and fully defined when both are known.
]

Certain boundary conditions are appropriate if we want well-defined solutions for the various classes of equation:
$
  "hyperbolic with open boundary surface" & tilde "Cauchy" \
   "parabolic with open boundary surface" & tilde "Dirichlet or Neumann" \
     "elliptic bounded by closed surface" & tilde "Dirichlet or Neumann"
$

=== Uniqueness
We prove the uniqueness of the Poisson equation when using either Dirichlet or Neumann boundary conditions.

#theorem[
  Let $u$ be real and let its first and second partial derivatives be continuous in a region $V subset.eq RR^n$ as well as on its boundary $S equiv dd(V, d: partial)$. Furthermore, let
  $
    laplacian u (bold(r)) = rho (bold(r)) "for" bold(r) in V
  $
  and either $u = f$ or $dd(u, d: partial)\/dd(n, d: partial) = g$ on $S$. Then the solution $u$ is unique.#footnote[Up to some constant.]
]

#proof[
  Let $u_1 (bold(r))$ and $u_2 (bold(r))$ be distinct solutions satisfying everything mentioned. We want to show that these differ by at most a constant. We define $omega(bold(r)) = u_1 (bold(r)) - u_2 (bold(r))$. Within $V$ this function satisfies the Laplace equation
  $
    laplacian omega = 0
  $
  and on $S$ we have
  $
    u_1 = u_2 = f " or " pdv(u_1, n) = pdv(u_2, n) = g
  $
  depending on which boundary conditions are given. Then we have
  $
    omega = 0 " or " pdv(omega, n) = 0
  $
  on $S$. We now consider
  $
    div (omega grad omega) = omega laplacian omega + grad omega dot grad omega
  $
  integrating
  $
    integral_V (underbracket(omega laplacian omega, 0) + grad omega dot grad omega) dd(V) &= integral_S omega (grad omega dot hat(n)) dd(S) \
    &= integral_S omega pdv(omega, n) dd(S) \
    &=^"by assumption" 0
  $
  We are left with
  $
    integral_V abs(grad omega)^2 dd(V) = 0 => grad omega = 0
  $
  implying $omega = u_1 - u_2 tilde "constant"$.
]

As a corollary it follows that for Dirichlet boundary conditions with $u_1 = u_2$ on $S$ then $u_1 = u_2$ everywhere in $V$. While for Neumann boundary conditions $u_1$ and $u_2$ can differ by a constant in $V$.
