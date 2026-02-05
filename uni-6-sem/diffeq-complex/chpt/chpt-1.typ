#import "chpt-temp.typ": *
#show: chpt-note.with()

= Ordinary Differential Equations
An _ordinary differential equation_ of order $n$ is an equation of the form
$
  f(x,y,dv(y, x),dots,dv(y, x, n)) = 0
$
with $y = y (x)$. We typically care about _linear_ differential equations in physics. These are of the form
$
  dv(y, x, n) = sum_(i=0)^(n-1) a_i (x) dv(y, x, i) + underbracket(r(x), "source")
$
with $a_i (x)$ and $r(x)$ being continuous. When $r(x) = 0$ we call the equation _homogeneous_, these always have the trivial solution $y(x) = 0$.

== Separation of variables
Consider an equation of the form
$
  dv(y, x) = F(x,y)
$
Assuming $F(x,y) = f(x) g(y)$ we can integrate the above giving
$
  integral f(x) dd(x) = integral 1/g(y) dd(y)
$
This is called _separation of variables_.


#example[Consider
  $
    dv(y, x) = x (1+y)
  $
  We find
  $
    integral x dd(x) & = integral dd(y)/(1+y) \
                     & => y = A e^(x^2\/2) - 1
  $
  with $A$ being an arbitrary constant.
]

== Exact equations
We can write
$
  F(x,y) = - A(x,y)/B(x,y)
$
Then
$
  A dd(x) + B dd(y) = 0
$
This differential is _exact_ if
$
  A = pdv(U, x)";  " B = pdv(U, y)
$
Meaning we have
$
  dd(U) equiv pdv(U, x) dd(x) + pdv(U, y) dd(y) =^! A dd(x) + B dd(y)
$
For some function $U(x,y)$. This is the case if
$
  pdv(A, y) = pdv(B, x)
$
Then $dd(U) = 0$ implying $U(x,y) = k$. We integrate the equation for $A$ to find
$
  U = integral A dd(x) + F(y)
$
To determine $F(y)$ we use the equation for $B$ giving
$
  B & = pdv(U, y) \
    & = integral pdv(A, y) dd(x) + dv(F, y)
$
which is quite nice.

#example[Consider
  $
    x dv(y, x) +3 x + y = 0
  $
  We find
  $
    A = 3 x +y";  " B = x
  $
  We have
  $
    pdv(A, y) = pdv(B, x)
  $
  so it is exact. Then
  $
    U & = integral 3x + y dd(x) + F(y) \
      & = 3/2 x^2 + x y + k_1 + F(y) =^! k_2
  $
  or
  $
    3/2 x^2 + x y + F(y) = k
  $
  We also have
  $
    underbracket(B, x) & = integral dd(x) + dv(F, y) \
                       & = x + k_1 + dv(F, y) \
                       & =^! pdv(U, y) = x + underbracket(dv(F, y), =^! 0)
  $
  so $k_1 = 0$ and
  $
    F(y) = k_2
  $
  Then
  $
    3/2 x^2 + x y = k
  $
  and we are done.
]

== Inexact equations
We now assume the equation is _inexact_ meaning
$
  A dd(x) + B dd(y) = 0";  " pdv(A, y) eq.not pdv(B, x)
$
We will try to make this exact by finding an _integrating factor_ $mu(x, y)$ such that
$
  pdv((mu A), y) = pdv((mu B), x)
$
This works since multiplying the differential by $mu$ does not change it. Assume $mu = mu(x)$ then
$
  mu pdv(A, y) = mu pdv(B, x) + B dv(mu, x)
$
or
$
  dd(mu)/mu = f(x) dd(x)
$
with
$
  f(x) equiv 1/B (pdv(A, y)-pdv(B, x))
$
We can integrate this to obtain
$
  mu(x) = exp[integral f(x) dd(x)]
$
Alternatively with $mu = mu(y)$ we find
$
  mu(y) = exp[integral g(y) dd(y)]
$
with
$
  g(y) equiv 1/A (pdv(B, x)-pdv(A, y))
$
Any arbitrary constant in the above has been thrown away since we just need some $mu$ that works.


#example[Consider
  $
    dv(y, x) + 2/y + (3y)/(2 x) = 0
  $
  where
  $
    A = 4 x+ 3 y^2";  " B=2 x y
  $
  This is clearly inexact. We find
  $
    f(x) = 2/x
  $
  giving
  $
    mu(x) = x^2
  $
  Then we multiply the differential by $x^2$ giving
  $
    (4 x^3 + 3 x^2 y^2) dd(x) + 2 x^3 y dd(y) = 0
  $
  which is exact by construction. We proceed as before
  $
    U & = integral A dd(x) + F(y) \
      & = x^4 + x^3 y^2 + k_1 + F(y) =^! k_2
  $
  so
  $
    x^4 + x^3 y^2 + F(y) = k
  $
  And
  $
    underbracket(B, 2 x^3 y) & = integral pdv(A, y) dd(x) + dv(F, y) \
                             & = integral 6 x^2 y dd(x) + dv(F, y) \
                             & = 2 x^3 y + k_1 + dv(F, y) \
                             & =^! pdv(U, y) = 2 x^3 y + underbracket(dv(F, y), 0)
  $
  so $k_1 = 0$ and
  $
    F(y) = k_2
  $
  Then
  $
    x^4 + x^3 y^2 = k
  $
  and we are done.
]
We now consider a special type of inexact ordinary differential equation. Namely those of the form
$
  dv(y, x) + P(x) y = Q (x)
$
which is a first order linear ordinary differential equation. We have
$
  A = P(x) y - Q(x)";  " B = 1
$
We find
$
  mu(x) = exp[integral P(x) dd(x)]
$
After making the differential exact we find
$
  Q(x) mu(x) & = mu(x) dv(y, x) + mu(x) P(x) y \
             & = dv(, x) (mu(x) y)
$
Since
$
  dv(, x) (mu(x) y) = mu(x) dv(y, x) + P(x) mu(x) y
$
We can then integrate to find
$
  y = 1/mu(x) integral Q(x) mu(x) dd(x)
$

#example[Consider
  $
    dv(y, x) + 2 x y = 4 x
  $
  We find
  $
    mu(x) = e^(x^2)
  $
  Then
  $
    y & = e^(-x^2) integral e^(x^2) 4 x dd(x) \
      & = 2 e^(-x^2) integral e^u dd(u) \
      & = 2 e^(-x^2) (e^(x^2) + k_1) \
      & = k e^(-x^2) + 2
  $
  and we are done.]
