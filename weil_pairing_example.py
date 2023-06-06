from sage.all_cmdline import *
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.schemes.elliptic_curves.ell_point import *
from sage.structure.factory import *

## BLS EXAMPLE
F2 = GF(2**28, 'b')
b = F2.gen()
E2 = EllipticCurve([F2(0), F2(0), F2(1), F2(1), F2(1)])
m = E2.order()
print(m)
n = 113
P = int(m/n**2)*E2.random_point()
Q = int(m/n**2)*E2.random_point()

print(m / n)
print(m / (n * n))
print(n * P)
print(n * Q)
print(P.order())
print(Q.order())

# The result is not 1
x = P.weil_pairing(Q, n)
print(x)

T = E2.random_point()
U = E2.random_point()

# The result is 1
y = T.weil_pairing(U, m)
print(y)
