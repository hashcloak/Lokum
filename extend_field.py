from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_point import *

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_cp6_782 import *
from test_pairing import *
    
q = 570227033427643938477351787526785557771164366858588273239167105706995178753794255646220989581954976296644199541099698255998850583874501725806067165976078609861
Fq = GF(q)
Fqz = Fq['z']; (z,) = Fqz._first_ngens(1)
print(Fqz)
Fq8_abs = Fq.extension(z**8 - 2, names=('w',)); (w,) = Fq8_abs._first_ngens(1);
print(Fq8_abs)
Fq2 = Fq.extension(z**2 - 2, names=('j',)); (j,) = Fq2._first_ngens(1);
print(Fq2)
Fq2Z = Fq2['Z']; (Z,) = Fq2Z._first_ngens(1);
print(Fq2Z)
Fq4 = Fq2.extension(Z**2 - j, names=('s',)); (s,) = Fq4._first_ngens(1)
print(Fq4)
Fq4W = Fq4['W']; (W,) = Fq4W._first_ngens(1)
print(Fq4W)
Fq8 = Fq4.extension(W**2 - s, names=('r',)); (r,) = Fq8._first_ngens(1)
print("x")
print(Fq8)

