from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

if __name__ == "__main__":
    arithmetic(False)

    QQx = QQ['x']; (x,) = QQx._first_ngens(1)

    r = ZZ(3618502788666131213697322783095070105623107215331596699973092056135872020481)
    q = ZZ(570227033427643938477351787526785557771164366858588273239167105706995178753794255646220989581954976296644199541099698255998850583874501725806067165976078609861)
    c_possible = ZZ(157586456811283429841422804972005719747057336671136400349968118675472488492656787210)
    c = q // r

    r_cp6_782 = ZZ(0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001)
    q_cp6_782 = ZZ(0x3848c4d2263babf8941fe959283d8f526663bc5d176b746af0266a7223ee72023d07830c728d80f9d78bab3596c8617c579252a3fb77c79c13201ad533049cfe6a399c2f764a12c4024bee135c065f4d26b7545d85c16dfd424adace79b57b942ae9)
    c_cp6_782 = q_cp6_782 // r_cp6_782
    c_res = ZZ(86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951561028)
    print(c_cp6_782)
    print(ZZ(c_cp6_782))
    print(c_res)