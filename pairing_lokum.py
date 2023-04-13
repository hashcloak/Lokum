from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_cp6_782 import *
from test_pairing import *


def ate_pairing_cp6_782(Q, P, a, tr):
    T = tr-1
    m,S1 = miller_function_ate(Q, P, a, T)
    f = final_exp_cp6_782(m)
    return f


def final_exp_hard_cp6_782(m):
    # returns m^(-W0 + q*W1) = (m^(-1))^W0 * (m^q)^W1
    W1_2naf=[0,1,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,-1,0,-1,0,0,1,0,0,-1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,1,0,0,0,0,-1,0,1,0,0,1,0,-1,0,-1,0,1,0,0,1,0,0,0,0,-1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,-1,0,0,1,0,0,0,1,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,0,-1,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,0,0,1]
    W0_2naf=[1,0,0,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,-1,0,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,-1,0,-1,0,0,-1,0,-1,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1,0,-1,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,-1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,-1,0,0,1,0,0,0,0,-1,0,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,-1,0,-1,0,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,-1,0,1,0,-1,0,-1,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,1,0,0,0,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,0,1,0,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,-1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,0,0,-1,0,-1,0,0,-1,0,-1,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,1,0,1,0,0,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,0,0,1,0,0,-1,0,1,0,0,-1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,-1,0,0,-1,0,0,0,-1,0,-1,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,-1,0,-1,0,1,0,0,0,0,0,1,0,-1,0,-1,0,0,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,0,0,0,-1,0,-1,0,-1,0,1,0,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,0,1,0,0,1,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,1,0,1]
    f0 = m
    f0_inv = f0.frobenius(3) # f0^(q^3) = 1/f0
    f1 = m.frobenius()
    f1_inv = f1.frobenius(3) # f1^(q^3) = 1/f1
    f,S,M = multi_exp_2naf(f0_inv,f0,f1,f1_inv,W0_2naf,W1_2naf,verbose=False)
    return f

def final_exp_cp6_782(m):
    f = final_exp_easy_k6(m)
    f = final_exp_hard_cp6_782(f)
    return f

def final_exp_easy_k6(m):
    """cost: f3 + i6 + 2*m6 + f
    """
    mp3 = m.frobenius(3) # assuming the extension is in one layer
    im = 1/m
    f = mp3*im # m^(q^3-1)
    f = f * f.frobenius() # m^((q^3-1)*(q+1))
    return f

def test_ate_pairing_cp6_782(E, E6_abs, r, c, E2, c2, tr, D_twist):
    q = E.base_field().characteristic()    
    exponent = ((q**6-1) // r)
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f = ate_pairing_cp6_782(Q, P, E.a4(), tr)
    
    f1_ate, S1  = miller_function_ate(Q,P,E.a4(), tr-1)
    f_ate = f1_ate**exponent
    f_ate_ = final_exp_cp6_782(f1_ate)
    ok = f_ate == f_ate_ == f
    # check it is bilinear
    #ok = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f_ate2 = ate_pairing_cp6_782(aa*Q, bb*P, E.a4(), tr)
            ok = ok and (f_ate2 == f**(aa*bb))
            bb += 1
        aa += 1
    print("check that ate pairing is bilinear: {}".format(ok))
    return ok    

def test_miller_function_ate_cp6_782(E, E6_abs, r, c, E2, c2, tr, D_twist=False):
    # ate pairing: impossible to check, impossible to compute ate with Magma or Sage
    q = E.base_field().characteristic()    
    exponent = ((q**6-1) // r)
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f1_ate, S1  = miller_function_ate(Q, P, E.a4(), tr-1)
    f_ate = f1_ate**exponent
    f_ate_ = final_exp_cp6_782(f1_ate)
    ok = f_ate == f_ate_
    # check it is bilinear
    #ok = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f_ate2, S1  =  miller_function_ate(aa*Q, bb*P, E.a4(), tr-1)
            f_ate2 = final_exp_cp6_782(f_ate2) # f_ate2^((q^6-1) // r)
            ok = ok and (f_ate2 == f_ate**(aa*bb))
            bb += 1
        aa += 1
    print("check that miller_function is bilinear: {}".format(ok))
    return ok

def test_curve(u, a, b, q, r, tr, c, y, D, xi):

    print("###################################################################")
    print("test curve E: y^2 = x^3 + a*x + b")
    print("a = {} b = {:#x}".format(a, b))

    assert q == (tr**2 + D*y**2)//4
    print("q = {:#x} # {} bits".format(q, q.nbits()))
    Fq = GF(q)
    E = EllipticCurve([Fq(a),Fq(b)])
    res = check_curve_order(E,r*c)
    print("check_curve_order(E,r*c): {}".format(res))
    print("ate pairing, Miller loop scalar is tr-1 =\nT = {}\nT = {:#x}".format(tr-1, tr-1))
    
    # absolute extension
    #Fqz.<z> = Fq[]
    #Fq6_abs.<w> = Fq.extension(z**6 - xi)
    #Fq3.<j> = Fq.extension(z**3 - xi)
    #Fq3Z.<Z> = Fq3[]
    #Fq6.<s> = Fq3.extension(Z**2 - j)
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1);
    Fq6_abs = Fq.extension(z**6 - xi, names=('w',)); (w,) = Fq6_abs._first_ngens(1);
    Fq3 = Fq.extension(z**3 - xi, names=('j',)); (j,) = Fq3._first_ngens(1);
    Fq3Z = Fq3['Z']; (Z,) = Fq3Z._first_ngens(1);
    Fq6 = Fq3.extension(Z**2 - j, names=('s',)); (s,) = Fq6._first_ngens(1)

    #E6 = EllipticCurve([Fq6(a),Fq6(b)])
    E6_abs = EllipticCurve([Fq6_abs(a),Fq6_abs(b)])
    
    # G2
    # M-twist
    E_M = EllipticCurve([(Fq3(a))*j**2, (Fq3(b))*j**3])
    print("M-twist\na=a*w^2={}\nb=b*w^3={}".format(E_M.a4(),E_M.a6()))
    E_D = EllipticCurve([(Fq3(a))/j**2, (Fq3(b))/j**3])
    print("D-twist\na=a/w^2={}\nb=b/w^3={}".format(E_D.a4(),E_D.a6()))
    # Sage does not know how to compute efficiently the curve order over that extension
    #order_EM = E_M.order()
    #order_ED = E_D.order()
    #assert (order_EM % r) == 0
    #assert (order_ED % r) == 0
    #assert order_EM == order_ED

    #tr1 = tr
    #tr2 = tr^2 - 2*q
    #tr3 = tr*tr2 - q*tr1
    tr3 = tr**3 - 3*q*tr # order E/Fq3 = (q^3 + 1 - tr3) = (q+1-tr) * (q^2-q+1 + t^2 + q*t + t)
    order_twist = q**3 + 1 + tr3
    assert order_twist % r == 0
    c2 = order_twist // r
    
    print("check_curve_order(E_M, order_twist): {}".format(check_curve_order(E_M, order_twist)))
    print("check_curve_order(E_D, order_twist): {}".format(check_curve_order(E_D, order_twist)))
    
    exponent = (q**2 - q + 1)//r
    e1 = exponent // q
    e0 = exponent % q
    
    e1_2naf = bits_2naf(e1)
    e0_2naf = bits_2naf(e0)    

# The above functions are used in the cp6_782 pairing impl. I need to transform it into lokum
# Use miller_function_ate for ate_pairing


# To use this function we need a point in Q, P, a, tr
def ate_pairing_lokum(Q, P, a, tr):
    T = tr-1
    m, S1 = miller_function_ate(Q, P, a, T)
    f = final_exp_lokum(m)
    return f

def final_exp_lokum(m):
    f = final_exp_easy_k8(m)
    f = final_exp_hard_lokum(f)
    return f

def final_exp_easy_k8(m):
    # The below parameters needs to be changed
    mp3 = m.frobenius(3)


def test_ate_pairing_lokum(E, E6_abs, r, c, E2, c2, tr, D_twist):
    q = E.base_field().characteristic()    
    exponent = ((q**6-1) // r) ## CHANGE q**8-1
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f = ate_pairing_cp6_782(Q, P, E.a4(), tr)
    
    f1_ate, S1  = miller_function_ate(Q,P,E.a4(), tr-1)
    f_ate = f1_ate**exponent
    f_ate_ = final_exp_cp6_782(f1_ate)
    ok = f_ate == f_ate_ == f
    # check it is bilinear
    #ok = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f_ate2 = ate_pairing_cp6_782(aa*Q, bb*P, E.a4(), tr)
            ok = ok and (f_ate2 == f**(aa*bb))
            bb += 1
        aa += 1
    print("check that ate pairing is bilinear: {}".format(ok))
    return ok    

def test_miller_function_ate_cp6_782(E, E6_abs, r, c, E2, c2, tr, D_twist=False):
    # ate pairing: impossible to check, impossible to compute ate with Magma or Sage
    q = E.base_field().characteristic()    
    exponent = ((q**6-1) // r)
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f1_ate, S1  = miller_function_ate(Q, P, E.a4(), tr-1)
    f_ate = f1_ate**exponent
    f_ate_ = final_exp_cp6_782(f1_ate)
    ok = f_ate == f_ate_
    # check it is bilinear
    #ok = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f_ate2, S1  =  miller_function_ate(aa*Q, bb*P, E.a4(), tr-1)
            f_ate2 = final_exp_cp6_782(f_ate2) # f_ate2^((q^6-1) // r)
            ok = ok and (f_ate2 == f_ate**(aa*bb))
            bb += 1
        aa += 1
    print("check that miller_function is bilinear: {}".format(ok))
    return ok
def check_curve_order(E, order):
    # randomized probabilistic check, usefull for large parameters
    ok = True
    i=0
    while ok and i < 10:
        P = E.random_element()
        ok = order*P == E(0)
        i += 1
    return ok

# def test_miller_function_ate_cp6_782(E, E6_abs, r, c, E2, c2, tr, D_twist=False):
def test_miller_function_ate_lokum(E, E8_abs, r, c, E2, c2, tr, D_twist=False):
    q = E.base_field().characteristic()    
    exponent = ((q**8-1) // r)
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f1_ate, S1  = miller_function_ate(Q, P, E.a4(), tr-1)
    f_ate = f1_ate**exponent
    f_ate_ = final_exp_cp6_782(f1_ate)
    ok = f_ate == f_ate_
    # check it is bilinear
    #ok = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f_ate2, S1  =  miller_function_ate(aa*Q, bb*P, E.a4(), tr-1)
            f_ate2 = final_exp_cp6_782(f_ate2) # f_ate2^((q^6-1) // r)
            ok = ok and (f_ate2 == f_ate**(aa*bb))
            bb += 1
        aa += 1
    print("check that miller_function is bilinear: {}".format(ok))
    return ok
# Originally -> def test_lokum(u, a, b, q, r, tr, c, y, D, xi):
def test_lokum(a, b, q, r, tr, c, D):
    print("###################################################################")
    print("test curve E: y^2 = x^3 + a*x + b")
    print("a = {} b = {:#x}".format(a, b))

    #assert q == (tr**2 + D*y**2)//4
    print("q = {:#x} # {} bits".format(q, q.nbits()))
    Fq = GF(q)
    E = EllipticCurve([Fq(a),Fq(b)])

    print("The order of the curve is:", c)
    res = check_curve_order(E,r*c)
    print("check_curve_order(E,r*c): {}".format(res))
    print("ate pairing, Miller loop scalar is tr-1 =\nT = {}\nT = {:#x}".format(tr-1, tr-1))

    # absolute extension
    #Fqz.<z> = Fq[]
    #Fq6_abs.<w> = Fq.extension(z**6 - xi)
    #Fq3.<j> = Fq.extension(z**3 - xi)
    #Fq3Z.<Z> = Fq3[]
    #Fq6.<s> = Fq3.extension(Z**2 - j)
    # ABOVE ARE COMMENTED IN THE ORIGINAL

    # COMMENTIN BELOW FOR CONSTRUCTING OWN
    #Fqz = Fq['z']; (z,) = Fqz._first_ngens(1);
    #Fq6_abs = Fq.extension(z**6 - xi, names=('w',)); (w,) = Fq6_abs._first_ngens(1);
    #Fq3 = Fq.extension(z**3 - xi, names=('j',)); (j,) = Fq3._first_ngens(1);
    #Fq3Z = Fq3['Z']; (Z,) = Fq3Z._first_ngens(1);
    #Fq6 = Fq3.extension(Z**2 - j, names=('s',)); (s,) = Fq6._first_ngens(1)

    ## THIS IS FOR EXTENDING LOKUM FIELDS 1-4-8
    xi = ZZ(2)
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1);
    Fq8_abs = Fq.extension(z**8 - xi, names=('j',)); (j,) = Fq8_abs._first_ngens(1)
    Fq4 = Fq.extension(z**4 - xi, names=('j',)); (j,) = Fq4._first_ngens(1);
    Fq4Z = Fq4['Z']; (Z,) = Fq4Z._first_ngens(1);
    Fq8 = Fq4.extension(Z**2 - j, names=('s',)); (s,) = Fq8._first_ngens(1)

    E8_abs = EllipticCurve([Fq8_abs(a),Fq8_abs(b)])
    print(E8_abs)

    E_D = EllipticCurve([(Fq4(a))/j**2, (Fq4(b))/j**3])
    print("D-twist\na=a/w^2={}\nb=b/w^3={}".format(E_D.a4(),E_D.a6()))
    print("PRINT FQ4")
    print(Fq4)
    N = [x for x in range(1000) if kronecker(x, q) == -1]; N
    print(N)
    print("The order of the twist is calculated...")
    order_twist = ZZ(105728290513165610986047842809751318841280324636548185919618761670766498489126319374080364643814474159560267765206481614470062820665803111802542081483890984079007411741658964681748397228590555002124666865154706881326247670437929905649005732882915063211059733358699029466871349150330576303159118955009483248411687688459258653412137668989908623154784117709434748451700616886945074033208009205147510709322838869285502070178802449805606834287959659872025482054893432653565073089295339866221304319849957430798002705359030400272401451744984691726540651157819242361547135100554273343028342409843532297688541312108594686251215425527313523033316)
    print(order_twist)
    assert order_twist % r == 0
    c2 = order_twist // r
    print("c2 is calculated..." )
    print(c2)
    print("check_curve_order(E_D, order_twist): {}".format(check_curve_order(E_D, order_twist)))
    
    exponent = (q**2 - q + 1)//r ## Is this correct for Lokum?
    e1 = exponent // q
    e0 = exponent % q
    
    e1_2naf = bits_2naf(e1)
    e0_2naf = bits_2naf(e0)
    print(e1_2naf)
    print(e0_2naf)
    
    ##

if __name__ == "__main__":
    print("Hello world!")

    arithmetic(False)
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BLS12-377 seed
    # CHANGE: u0
    u0 = ZZ(0x8508C00000000001)
    # Lokum parameters
    r = ZZ(3618502788666131213697322783095070105623107215331596699973092056135872020481)
    q = ZZ(570227033427643938477351787526785557771164366858588273239167105706995178753794255646220989581954976296644199541099698255998850583874501725806067165976078609861)
    a = ZZ(408825162639581376425354251544709546275292709793540929634396921208238440519285149826351823244598728486772571705638177322608711025414928328755239670937630810045)
    b = ZZ(349808739615329235353531183724887241300117306968432494640082668760849877387756346906668495855253497480229722131408639028906859040081069711950364714123876057189)
    # cofactor of the curve
    c = ZZ(157586456811283429841422804972005719747057336671136400349968118675472488492656786156)
    tr= ZZ(3816706629471577761299121574721828830967739654627469433190085390763284409348826)
    D = ZZ(-143)

    Miller_loop_scalar = ZZ(3816706629471577761299121574721828830967739654627469433190085390763284409348825)
    assert Miller_loop_scalar == (tr-1)
    ## COMMENT FROM CP6 T = 506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621
    ## COMMENT FROM CP6 T = 0xcda91b419cd3f63a81dba4eb8fc81231d84f4491c1e1a4dfaf393e5c14beec4affbc93fba900003f8fbe3c000000007a5
    # absolute extension
    ## CHANGE? What's xi? I think it's the generator in Fq file? But, what's the generator?
    #xi = 13
    
    ## CHANGE? What's tr3
    #tr3=ZZ(-0x87219fb9ab6223dea08125f7b21a9d72bfd2eca93424662b01c6ba77d0efc9ada28e36672abc6ddb9ff90435350b0907a055422d2327d4b7662461c32b1a122e8cd4fdcbdcb496d5395acf3575580d5ef3b2859d2fb197730dcd15e9e306306baa548e380d51e48cf9f029557488a5d8e51a4d13319fa7b0ecf8d1613fd2fca09e981bf2d5f42d62a7ef0e68ce62bce5455ea)
    # we have r^2 divides q^6+1-t6 = (q3+1-t3)*(q^3+1+t3),
    # r divides (q+1-t) and (q+1-t) divides #E(Fq3)=(q3+1-t3)
    # r divides also the quadratic twist of E(Fq3) of order q^3+1+t3
    #order_twist = q^3 + 1 + tr3
    ## CHANGE? How to find the order_twist
    #order_twist = ZZ(0x2b87fda17134c96a73372daf9a1407c1bfe240a89db104c580df9d215037525cbc4ba24cf9ce08b93a5652a8e16a97fc618f319fb7fa85b5a9e180cd3f8dabebe38a09b261651d60f8f9f6fcd973093d8719c862c2bcd8dc3b1986ef18419cd0d1d6be027e0248d61e0c45d27076027ec67b23ca21268d5c2076b44f0a32c9a8abe9b50297df4c5efe04c0cec1019fb862c0494c868105fb739b5347526f005a82c4cb19e83b1a26fa18378c42661feec7f593c522696a885ca18780d64da536e647fec41f299502df523124f0b29870a703e5bde1bc49259646f3bbd388ab5642f80a80a89c72a4f78017ce6e9e8ccb996212d7771fa8a4fc5500b2352b94d58a035f2c4635f33ff5435f393055e21a030ae458054a1372030868e0e0540a3fca151d50b90)
    #assert (order_twist % r) == 0
    
    #c2=ZZ(0x19e70d3618ca3a1ad2ad61580ef78d7f8e3ef002559614ed25f7ebe6c44b3be7ce44dd084e05d4a86d8ebb7fdf52ad09033887d09b4d2c4b99b0d9b077d21a3b8a42c4e1e6f5a9eea59baa2aeb4a1fee44392814554748b589c61531dbce7d1d80da227dd39c1d8f9eda95d6abf3ef8c39cbf0df977800f25b044f4e9c9327ecaaf050acb91f3a597504275948226574956a099ace2c3c541018fe4d10cc7bedd37f629e8e6348eb48c92998f3358bb74eb9b6630846b7a5e400154f633cf4216f0595cee44d02cea7dc01aae2140556e918748f067931d6f388069c063c11f0c9a51593693f88c98a12bd486d2fb4b77fca151d50b90)
    #assert c2*r == order_twist
    # Final exponentiation
    # hard part: (q^2-q+1)/r = W0 + W1*q
    #W0 = ZZ(7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033)
    #W1 = ZZ(86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986)

    # easy final exp costs 4 + 34 + inversion + 18 = 56 + inversion ~= 81 m
    #W1_2naf=[0,1,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,-1,0,-1,0,0,1,0,0,-1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,1,0,0,0,0,-1,0,1,0,0,1,0,-1,0,-1,0,1,0,0,1,0,0,0,0,-1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,-1,0,0,1,0,0,0,1,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,0,-1,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,0,0,1]
    #W0_2naf=[1,0,0,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,-1,0,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,-1,0,-1,0,0,-1,0,-1,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1,0,-1,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,-1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,-1,0,0,1,0,0,0,0,-1,0,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,-1,0,-1,0,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,-1,0,1,0,-1,0,-1,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,1,0,0,0,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,0,1,0,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,-1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,0,0,-1,0,-1,0,0,-1,0,-1,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,1,0,1,0,0,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,0,0,1,0,0,-1,0,1,0,0,-1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,-1,0,0,-1,0,0,0,-1,0,-1,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,-1,0,-1,0,1,0,0,0,0,0,1,0,-1,0,-1,0,0,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,0,0,0,-1,0,-1,0,-1,0,1,0,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,0,1,0,0,1,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,1,0,1]

    #exponent = (q**2 - q + 1)//r
    #e1 = exponent // q
    #e0 = exponent % q
    #e1 += 1
    #e0 -= q
    #print("e1 = {}\ne0 = {}\nW1 = {}\nW0 = {}".format(e1, e0, W1, W0))
    #assert -W0 + W1*q == exponent
    #log2_e = float(log(abs(e1), float(2.0)) + log(abs(e0), float(2.0)))
    #log2_W = float(log(abs(W1), float(2.0)) + log(abs(W0), float(2.0)))
    #print("log_2(e1) + log_2(e0) = {}".format(log2_e))
    #print("log_2(W1) + log_2(W0) = {}".format(log2_W))
    #M = Matrix(ZZ, 2, 2, [exponent, 0, -q, 1])
    #R = M.LLL()
    #print(R)
    #for i in range(2):
    #    print("R[{}][0] + q*R[{}][1] = ({}) * ((q^2-q+1)//r) + ({})".format(i,i, (R[i][0] + R[i][1]*q) // exponent, (R[i][0] + R[i][1]*q) % exponent ))
    #log2_R0 = float(log(abs(R[0][0]), float(2.0)) + log(abs(R[0][1]), float(2.0)))
    #print("log_2(R00) + log_2(R01) = {}".format(log2_R0))
    #log2_R1 = float(log(abs(R[1][0]), float(2.0)) + log(abs(R[1][1]), float(2.0)))
    #print("log_2(R10) + log_2(R11) = {}".format(log2_R1))
    #assert W0 == -e0
    #assert W1 == e1
    test_lokum(a, b, q, r, tr, c, D)