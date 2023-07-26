"""
author: Aurore Guillevic <aurore.guillevic@inria.fr>
File on generating Cocks-Pinch curves from a prime r,
parameters from Mikerah <mikerah@hashcloak.com>
date: March 8, 2023
"""
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.number_field.number_field import NumberField

#from tnfs.gen.generate_curve_utils import str_py_binary_form

q = ZZ(19958386796366484580525069329849932218129718221411371521433091928786090572515874127845723606834618339413069824836929295092920273020066744418556009330309419)
r = ZZ(3618502788666131213697322783095070105623107215331596699973092056135872020481)
k = 5
tr = ZZ(162523082385094613493926538014996591663598424186107575522851399605278493878876)
D = 7
j = hilbert_class_polynomial(-7).roots()[0][0]

# what are the parameters that led to q? that is, find q = ((tr0 + ht*r)^2 + 7*(y0+hy*r)^2)/4
t0 = (tr % r)
zeta_5 = t0-1
assert (zeta_5^5) % r == 1
ht = tr//r
assert tr == t0 + ht*r
y = sqrt((4*q-tr^2)//7)
y0 = y % r
hy = y // r
assert y == y0 + hy*r
assert (tr^2 + D*y^2) % 4 == 0
assert q == (tr^2 + D*y^2)//4

#str_py_binary_form(r)
r == 2**251+2**196+2**192+1
assert r.is_prime()
Fr = GF(r)
QQx.<x> = QQ[]
rx = x^251 + x^196 + x^192 + 1
assert rx.is_irreducible()
K.<a> = NumberField(rx)
(K(-7)).is_square() # -7 is not a square so there will be no "nice" formula for sqrt(-D) mod r

sqrt_7 = ZZ(sqrt(Fr(D)))
sqrt_7_ = r-sqrt_7
L.<z> = NumberField(x^2+x+2)
L.discriminant()

# check that rx does not factor over L but r factors (into conjugate prime ideals above r)
rx.roots(L)
# equivalently:
LX.<X> = L[]
rx(X).factor()

L(q).factor()
L(r).factor()

t0 = Fr.random_element()
t0 = t0^((r-1)//5)
while t0 == Fr(1):
    t0 = Fr.random_element()
    t0 = t0^((r-1)//5)
print("t0 = {} is a 5th primitive root of unity".format(t0))
print("sqrt(-7) mod r = {}".format(sqrt_7))
max_h = 100
for (ti,ei) in [(t0,1), (t0^2,2), (t0^3,3), (t0^4,4)]:
    t_ = ZZ(ti)
    y_ = ZZ((ti-2)*sqrt_7)
    print("t0 = {}\ny0 = {}".format(t_, y_))
    ht = -max_h
    for t1 in range(t_-max_h*r, t_+max_h*r, r):
        hy = -max_h
        # there are two choices of sqrt_7: sqrt_7 and r-sqrt_7
        # y_+hy*r and r-y_+h_y*r = -(y_-(h_y+1)*r)
        # y_-hy*r is not exactly the same as r-y_+hy*r
        for y1 in range(y_-max_h*r, y_+max_h*r, r):
            p4 = (t1^2 + D*y1^2)
            if (p4 % 4) == 0:
                p = p4 // 4
                if p.nbits() <= 512:
                    if p.is_prime():
                        print("found prime q of {} bits at\nt= t0^{}{:+d}*r, \nt={}\ny=y0{:+d}*r\ny={}\np={}\n".format(p.nbits(), ei, ht, t1, hy, y1, p))
            hy = hy+1
        ht = ht+1
# many alternative Cocks-Pinch curves with q smaller than 512 bits.
"""
t0 = 3308959683784840091244335558813507016181706711517320724035349135300124977711 is a 5th primitive root of unity
sqrt(-7) mod r = 209880913220488156313834499058596184255382000503803956865172539718143124843
t0 = 3308959683784840091244335558813507016181706711517320724035349135300124977711
y0 = 972226531867112912424057242609662492407766450545671205315101752060443862612

found prime q of 509 bits at
t= t0^1+17*r, 
t=64823507091109070724098822871429698811774529372154464623577914089609949325888
y=y0+2*r
y=8209232109199375339818702808799802703653980881208864605261285864332187903574
p=1168456878587509537335431496700422980516633904232772247889202989338172690739525702652736357629558044514068566944158319258680469138981368241571362370200719

found prime q of 509 bits at
t= t0^3+17*r, 
t=62518265236988912029814048515112361934391985679315883489943096793625832170434
y=y0+0*r
y=725387255838288000852399593166548692797099371264050037953575361856639912248
p=978054198734756119518149330319389679799322530504493849515846609856692151993582595018913029551254115170884403801923256897302100747548522004514989801660721

"""
