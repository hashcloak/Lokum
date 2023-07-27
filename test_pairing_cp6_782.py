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

def check_curve_order(E, order):
    # randomized probabilistic check, usefull for large parameters
    ok = True
    i=0
    while ok and i < 10:
        P = E.random_element()
        ok = order*P == E(0)
        i += 1
    return ok

# here maps from Fq6 to Fq6_abs manually otherwise pb with Sage, no frobenius power available
def map_psi_M_D_twist(Q2, E6_abs, D_twist=False):
    Fq6_abs_ext = E6_abs.base_field()
    w = Fq6_abs_ext.gen(0)
    Q2cx = (Q2[0]).polynomial().coefficients() # Q2[0] in Fq3 -> u0 + u1*x + u2*x^2
    Q2x = Fq6_abs_ext([Q2cx[0], 0, Q2cx[1], 0, Q2cx[2], 0])
    Q2cy = (Q2[1]).polynomial().coefficients()
    Q2y = Fq6_abs_ext([Q2cy[0], 0, Q2cy[1], 0, Q2cy[2], 0])
    if not D_twist:
        Q = psi_sextic_m_twist([Q2x,Q2y], w)
    else:
        Q = psi_sextic_d_twist([Q2x,Q2y], w)
    return E6_abs(Q)

def test_final_exp_basic_cp6_782(W0,W1, Fq6_abs):
    i = 0
    ok1 = True
    exponent = q**3-1
    while ok1 and i < 10:
        m = Fq6_abs.random_element()
        mp3 = m.frobenius(3)
        im = 1/m
        f = mp3*im # m^(p^3-1)
        ok1 = f == m**exponent # ok1 = f == m**(q**3-1)
        i = i+1
    print("check easy part of final exp: f^(q^3-1): {}".format(ok1))
    i = 0
    ok2 = True
    exponent = (q**3-1)*(q+1)
    while ok2 and i < 10:
        m = Fq6_abs.random_element()
        mp3 = m.frobenius(3)
        im = 1/m
        f0 = mp3*im # m^(p^3-1)
        f = f0 * f0.frobenius() # m^((q^3-1)*(q+1))
        #ok2 = f == (m**(q**3-1))**(q+1)
        ok2 = (f == f0**(q+1)) and (f == m**exponent)
        i = i+1
    print("check easy part of final exp: f^((q^3-1)*(q+1)) {}".format(ok2))
    i = 0
    ok3 = True
    exponent = (q**3-1)*(q+1)
    while ok3 and i < 10:
        m = Fq6_abs.random_element()
        mp3 = m.frobenius(3)
        im = 1/m
        f0 = mp3*im # m^(p^3-1)
        f = f0 * f0.frobenius() # m^((q^3-1)*(q+1))
        # 1 = f^(p^6-1) = f^((p-1)*(1+p+p^2+p^3+p^4+p^5))
        f_1 = f.frobenius(3) # now 1/f == f^(p^3)
        ok3 = f_1 == 1/f
        i = i+1
    print("check f^p^3 = 1/f: {}".format(ok3))
    i = 0
    ok4 = True
    exponent = ((q**6-1) // r)
    while ok4 and i < 10:
        m = Fq6_abs.random_element()
        mp3 = m.frobenius(3)
        im = 1/m
        f0 = mp3*im # m^(p^3-1)
        f = f0 * f0.frobenius() # m^((q^3-1)*(q+1))
        # 1 = f^(p^6-1) = f^((p-1)*(1+p+p^2+p^3+p^4+p^5))
        f_1 = f.frobenius(3) # now 1/f == f^(p^3)
        f = f_1**W0 * f.frobenius()**W1
        ok4 = f == m**exponent
        i = i+1
    print("check final exp: {}".format(ok4))
    ok = ok1 and ok2 and ok3 and ok4
    print("check basic final exp: {}".format(ok))
    return ok

def test_multi_exp(Fq6_abs, W0, W1):
    ok = True
    i = 0
    while ok and i < 10:
        f0 = Fq6_abs.random_element()
        f1 = Fq6_abs.random_element()
        ff = multi_exp(f0, W0, f1, W1)
        ft = f0**W0 * f1**W1
        ok = ft == ff
        i += 1
    print("check multi-exp: {}".format(ok))
    return ok

def test_multi_exp_2naf(Fq6_abs, W0_2naf, W1_2naf):
    ok = True
    i = 0
    while ok and i < 10:
        f0 = Fq6_abs.random_element()
        f1 = Fq6_abs.random_element()
        ff,M,S = multi_exp_2naf(f0, 1/f0, f1, 1/f1, W0_2naf, W1_2naf)
        ft = f0**W0 * f1**W1
        ok = ft == ff
        i += 1
    print("check multi-expo 2naf: {}".format(ok))
    return ok

def test_final_exp_cp6_782(Fq6_abs, r):
    q = Fq6_abs.characteristic()
    exponent = (q**6-1) // r
    ok =True
    i = 0
    while ok and i < 10:
        f = Fq6_abs.random_element()
        g = final_exp_cp6_782(f)
        h = f**exponent
        ok = g == h
        i += 1
    print("check final exp: {}".format(ok))
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

def test_tate_pairing_cp6_782(E, E6_abs, r, c, E2, c2, tr, D_twist=False):
    # swap P and Q of ate function to get a Tate Miller loop
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    #assert r*Q == E6_abs(0)
    P = c*E.random_element() # now P has order r
    #assert r*P == E(0)
    f1_tate, S2 = miller_function_tate(P, Q, E.a4(), r)
    #f_tate = f_tate^((q^6-1) // r)
    f_tate = final_exp_cp6_782(f1_tate)
    #print("f_tate = {}".format(f_tate))

    # check it is bilinear and non-degenerate
    Fq1 = f_tate.parent()(1)
    ok = True
    ok_non_degenerate = True
    aa = 1
    while ok and aa < 6:
        bb = 1
        while ok and bb < 6:
            f1_tate2, S1  = miller_function_tate(aa*P, bb*Q, E.a4(), r)
            f_tate2 = final_exp_cp6_782(f1_tate2) # f1_tate2^((q^6-1) // r)
            ok_non_degenerate = ok_non_degenerate and f_tate2 != Fq1
            ok = f_tate2 == f_tate**(aa*bb)
            bb = bb + 1
        aa = aa + 1
    print("check that Tate pairing is bilinear: {}".format(ok))
    print("check that Tate pairing is non-degenerate: {}".format(ok_non_degenerate))
    return ok

def test_ate_csb(E, E6_abs, E2, r, c, c2, tr, D_twist=False):
    i=0
    ok = True
    while ok and i < 10:
        P = c*E.random_element() # now P has order r
        Q2 = c2*E2.random_element()
        Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
        f_ate = ate_pairing_cp6_782(Q, P, E.a4(), tr)
        f1_ate_csb, S1  = miller_function_ate_csb(Q, P, E.a4(), tr-1)
        f_ate_csb = final_exp_cp6_782(f1_ate_csb)
        #print("f_ate_csb == f_ate: {}".format(f_ate_csb == f_ate))
        ok1 = f_ate_csb == f_ate
        #t_1Q = (tr-1)*Q
        t_1Q2 = (tr-1)*Q2
        t_1Q = map_psi_M_D_twist(t_1Q2, E6_abs, D_twist)
        #print("S1 == (tr-1)*Q: {}".format(E6_abs((S1[0]/S1[2]^2, S1[1]/S1[2]^3)) == t_1Q))
        ok2 = S1[0]/S1[2]**2 == t_1Q[0] and S1[1]/S1[2]**3 == t_1Q[1]
        ok = ok1 and ok2
        i += 1
    print("check miller_function_ate_csb: {}".format(ok))
    return ok
    
def test_tate_csb(E, E6_abs, E2, r, c, c2, D_twist=False):
    i=0
    ok = True
    while ok and i < 10:
        P = c*E.random_element() # now P has order r
        Q2 = c2*E2.random_element()
        Q = map_psi_M_D_twist(Q2, E6_abs,D_twist)
        f_tate = tate_pairing_cp6_782(P, Q, E.a4(), r)
        f1_tate_csb, S1  = miller_function_tate_csb(P, Q, E.a4(), r)
        f_tate_csb = final_exp_cp6_782(f1_tate_csb)
        ok1 = f_tate_csb == f_tate
        #print("f_tate_csb == f_tate: {}".format(ok1))
        #print("S1 == O: {}".format(S1[2] == 0))
        ok2 = S1[2] == 0
        ok = ok1 and ok2
        i += 1
    print("check miller_function_tate_csb: {}".format(ok))
    return ok

def test_g2_frobenius_eigenvalue_cp6(E2,E6_abs,c2,D_twist=False):
    ok = True
    i=0
    while ok and i < 10:
        Q2 = c2*E2.random_element()
        Q = map_psi_M_D_twist(Q2,E6_abs,D_twist)
        ok = Q[0].frobenius() == Q[0]**q
        Qp = [(Q[0]).frobenius(), Q[1].frobenius()]
        Qp = E6_abs(Qp)
        ok = ok and q*Q == Qp and (tr-1)*Q == Qp
        i += 1
    print("check that Frobenius(x) = x^q, pi(Q) == q*Q == (tr-1)*Q where r*Q == E(0): {}".format(ok))
    return ok

def test_miller_loop_opt_ate_cp6_782(E, E6_abs, E2, r, c, c2, u0, D_twist):
    P = c*E.random_element() # now P has order r
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    f1 = miller_loop_opt_ate_cp6_782(Q, P, E.a4(), u0)
    g1 = final_exp_cp6_782(f1)
    ok = True
    aa=1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            f2 = miller_loop_opt_ate_cp6_782(bb*Q, aa*P, E.a4(), u0)
            g2 = final_exp_cp6_782(f2)
            ok = ok and g2 == g1**(aa*bb)
            bb += 1
        aa += 1
    print("check miller_loop_opt_ate_cp6_782 is bilinear: {}".format(ok))
    return ok

def test_miller_loop_opt_ate_cp6_782_2naf(E, E6_abs, E2, r, c, c2, u0, D_twist):
    P = c*E.random_element() # now P has order r
    Q2 = c2*E2.random_element()
    Q = map_psi_M_D_twist(Q2, E6_abs, D_twist)
    f1 = miller_loop_opt_ate_cp6_782_2naf(Q, P, E.a4(), u0)
    g1 = final_exp_cp6_782(f1)
    ok = True
    aa=1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            f2 = miller_loop_opt_ate_cp6_782_2naf(bb*Q, aa*P, E.a4(), u0)
            g2 = final_exp_cp6_782(f2)
            ok = ok and g2 == g1**(aa*bb)
            bb += 1
        aa += 1
    print("check miller_loop_opt_ate_cp6_782_2naf_test is bilinear: {}".format(ok))
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
    print("--------------- FIELD EXTENSIONS -------------")
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1)
    print(Fqz)
    Fq6_abs = Fq.extension(z**6 - xi, names=('w',)); (w,) = Fq6_abs._first_ngens(1);
    print(Fq6_abs)
    Fq3 = Fq.extension(z**3 - xi, names=('j',)); (j,) = Fq3._first_ngens(1);
    print(Fq3)
    Fq3Z = Fq3['Z']; (Z,) = Fq3Z._first_ngens(1);
    print(Fq3Z)
    Fq6 = Fq3.extension(Z**2 - j, names=('s',)); (s,) = Fq6._first_ngens(1)
    print(Fq6)

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
    #print(order_ED)
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
    
    test_final_exp_basic_cp6_782(W0,W1, Fq6_abs)
    test_multi_exp(Fq6_abs,W0,W1)
    test_multi_exp_2naf(Fq6_abs,W0_2naf,W1_2naf)
    test_final_exp_cp6_782(Fq6_abs,r)

    print("test Miller function, ate and Tate pairing, M-twist")
    test_miller_function_ate_cp6_782(E,E6_abs,r,c,E_M,c2,tr,D_twist=False)
    test_ate_pairing_cp6_782(E,E6_abs,r,c,E_M,c2,tr,D_twist=False)
    test_tate_pairing_cp6_782(E,E6_abs,r,c,E_M,c2,tr,D_twist=False)

    print("test Miller function, ate and Tate pairing, D-twist")
    test_miller_function_ate_cp6_782(E,E6_abs,r,c,E_D,c2,tr,D_twist=True)
    test_ate_pairing_cp6_782(E,E6_abs,r,c,E_D,c2,tr,D_twist=True)
    test_tate_pairing_cp6_782(E,E6_abs,r,c,E_D,c2,tr,D_twist=True)

    # test CSB formulas
    print("test formulas from CSB'04, M-twist")
    test_ate_csb(E,E6_abs,E_M,r,c,c2,tr,D_twist=False)
    test_tate_csb(E,E6_abs,E_M,r,c,c2,D_twist=False)

    print("test formulas from CSB'04, D-twist")
    test_ate_csb(E,E6_abs,E_D,r,c,c2,tr,D_twist=True)
    test_tate_csb(E,E6_abs,E_D,r,c,c2,D_twist=True)

    ## test optimal ate pairing on CP6_782
    print("Test optimal ate pairing on CP6_782 curve, same Miller loop formula as for BW6_761")
    assert ((u0**3-u0**2-u0) + q*(u0+1)) % r == 0
    test_g2_frobenius_eigenvalue_cp6(E_M,E6_abs,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_cp6(E_D,E6_abs,c2,D_twist=True)

    test_miller_loop_opt_ate_cp6_782(E,E6_abs,E_M,r,c,c2,u0,D_twist=False)
    test_miller_loop_opt_ate_cp6_782_2naf(E,E6_abs,E_M,r,c,c2,u0,D_twist=False)

    test_miller_loop_opt_ate_cp6_782(E,E6_abs,E_D,r,c,c2,u0,D_twist=True)
    test_miller_loop_opt_ate_cp6_782_2naf(E,E6_abs,E_D,r,c,c2,u0,D_twist=True)

if __name__ == "__main__":
    arithmetic(False)
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BLS12-377 seed
    u0 = ZZ(0x8508C00000000001)
    # CP6_782 parameters
    r = ZZ(0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001)
    q = ZZ(0x3848c4d2263babf8941fe959283d8f526663bc5d176b746af0266a7223ee72023d07830c728d80f9d78bab3596c8617c579252a3fb77c79c13201ad533049cfe6a399c2f764a12c4024bee135c065f4d26b7545d85c16dfd424adace79b57b942ae9)
    #https://github.com/scipr-lab/zexe/blob/master/algebra/src/curves/cp6_782/g1.rs
    a = ZZ(5)
    b = ZZ(17764315118651679038286329069295091506801468118146712649886336045535808055361274148466772191243305528312843236347777260247138934336850548243151534538734724191505953341403463040067571652261229308333392040104884438208594329793895206056414)
    # cofactor of the curve
    c = ZZ(86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951561028)
    tr= ZZ(0xcda91b419cd3f63a81dba4eb8fc81231d84f4491c1e1a4dfaf393e5c14beec4affbc93fba900003f8fbe3c000000007a6)
    D = ZZ(339)
    y = ZZ(0xd0530d574d39bcf8401f2ad25ea64d565eaea5bd8a7d14ce7685a7532d5e468ffb1c299c166f9ca31cc4fed54e4b2ba60)
    Miller_loop_scalar = ZZ(506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621)
    assert Miller_loop_scalar == (tr-1)
    #T = 506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621
    #T = 0xcda91b419cd3f63a81dba4eb8fc81231d84f4491c1e1a4dfaf393e5c14beec4affbc93fba900003f8fbe3c000000007a5
    # absolute extension
    xi = 13
    
    tr3=ZZ(-0x87219fb9ab6223dea08125f7b21a9d72bfd2eca93424662b01c6ba77d0efc9ada28e36672abc6ddb9ff90435350b0907a055422d2327d4b7662461c32b1a122e8cd4fdcbdcb496d5395acf3575580d5ef3b2859d2fb197730dcd15e9e306306baa548e380d51e48cf9f029557488a5d8e51a4d13319fa7b0ecf8d1613fd2fca09e981bf2d5f42d62a7ef0e68ce62bce5455ea)
    # we have r^2 divides q^6+1-t6 = (q3+1-t3)*(q^3+1+t3),
    # r divides (q+1-t) and (q+1-t) divides #E(Fq3)=(q3+1-t3)
    # r divides also the quadratic twist of E(Fq3) of order q^3+1+t3
    #order_twist = q^3 + 1 + tr3
    order_twist = ZZ(0x2b87fda17134c96a73372daf9a1407c1bfe240a89db104c580df9d215037525cbc4ba24cf9ce08b93a5652a8e16a97fc618f319fb7fa85b5a9e180cd3f8dabebe38a09b261651d60f8f9f6fcd973093d8719c862c2bcd8dc3b1986ef18419cd0d1d6be027e0248d61e0c45d27076027ec67b23ca21268d5c2076b44f0a32c9a8abe9b50297df4c5efe04c0cec1019fb862c0494c868105fb739b5347526f005a82c4cb19e83b1a26fa18378c42661feec7f593c522696a885ca18780d64da536e647fec41f299502df523124f0b29870a703e5bde1bc49259646f3bbd388ab5642f80a80a89c72a4f78017ce6e9e8ccb996212d7771fa8a4fc5500b2352b94d58a035f2c4635f33ff5435f393055e21a030ae458054a1372030868e0e0540a3fca151d50b90)
    #assert (order_twist % r) == 0
    
    c2=ZZ(0x19e70d3618ca3a1ad2ad61580ef78d7f8e3ef002559614ed25f7ebe6c44b3be7ce44dd084e05d4a86d8ebb7fdf52ad09033887d09b4d2c4b99b0d9b077d21a3b8a42c4e1e6f5a9eea59baa2aeb4a1fee44392814554748b589c61531dbce7d1d80da227dd39c1d8f9eda95d6abf3ef8c39cbf0df977800f25b044f4e9c9327ecaaf050acb91f3a597504275948226574956a099ace2c3c541018fe4d10cc7bedd37f629e8e6348eb48c92998f3358bb74eb9b6630846b7a5e400154f633cf4216f0595cee44d02cea7dc01aae2140556e918748f067931d6f388069c063c11f0c9a51593693f88c98a12bd486d2fb4b77fca151d50b90)
    #assert c2*r == order_twist
    # Final exponentiation
    # hard part: (q^2-q+1)/r = W0 + W1*q
    W0 = ZZ(7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033)
    W1 = ZZ(86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986)

    # easy final exp costs 4 + 34 + inversion + 18 = 56 + inversion ~= 81 m
    W1_2naf=[0,1,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,-1,0,-1,0,0,1,0,0,-1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,1,0,0,0,0,-1,0,1,0,0,1,0,-1,0,-1,0,1,0,0,1,0,0,0,0,-1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,-1,0,0,1,0,0,0,1,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,0,-1,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,0,0,1]
    W0_2naf=[1,0,0,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,-1,0,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,-1,0,-1,0,0,-1,0,-1,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1,0,-1,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,-1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,-1,0,0,1,0,0,0,0,-1,0,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,-1,0,-1,0,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,-1,0,1,0,-1,0,-1,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,1,0,0,0,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,0,1,0,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,-1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,0,0,-1,0,-1,0,0,-1,0,-1,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,1,0,1,0,0,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,0,0,1,0,0,-1,0,1,0,0,-1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,-1,0,0,-1,0,0,0,-1,0,-1,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,-1,0,-1,0,1,0,0,0,0,0,1,0,-1,0,-1,0,0,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,0,0,0,-1,0,-1,0,-1,0,1,0,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,0,1,0,0,1,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,1,0,1]

    exponent = (q**2 - q + 1)//r
    e1 = exponent // q
    e0 = exponent % q
    e1 += 1
    e0 -= q
    print("e1 = {}\ne0 = {}\nW1 = {}\nW0 = {}".format(e1, e0, W1, W0))
    assert -W0 + W1*q == exponent
    log2_e = float(log(abs(e1), float(2.0)) + log(abs(e0), float(2.0)))
    log2_W = float(log(abs(W1), float(2.0)) + log(abs(W0), float(2.0)))
    print("log_2(e1) + log_2(e0) = {}".format(log2_e))
    print("log_2(W1) + log_2(W0) = {}".format(log2_W))
    M = Matrix(ZZ, 2, 2, [exponent, 0, -q, 1])
    R = M.LLL()
    print(R)
    for i in range(2):
        print("R[{}][0] + q*R[{}][1] = ({}) * ((q^2-q+1)//r) + ({})".format(i,i, (R[i][0] + R[i][1]*q) // exponent, (R[i][0] + R[i][1]*q) % exponent ))
    log2_R0 = float(log(abs(R[0][0]), float(2.0)) + log(abs(R[0][1]), float(2.0)))
    print("log_2(R00) + log_2(R01) = {}".format(log2_R0))
    log2_R1 = float(log(abs(R[1][0]), float(2.0)) + log(abs(R[1][1]), float(2.0)))
    print("log_2(R10) + log_2(R11) = {}".format(log2_R1))
    assert W0 == -e0
    assert W1 == e1

    test_curve(u0, a, b, q, r, tr, c, y, D, xi)


    #### this is possible to find a=-3 with another j-invariant, this is the curve below
    #An isogeneous curve to the Short Weierstrass curve SW6
    #such that a=-3
    #the curve has same order but distinct j-invariant
    #to find this curve, one computes the roots of hilbert_class_polynomial(-339)
    #the six roots are j-invariants of isogeneous curves, one is SW6 j-invariant
    a = -3
    b = 0x5d7ed1aa96167fb80c49c6724722ef3644eb763d432dbd9f89d6183d381098596b5eb1631e856b8fa0ea3954f036f79b658e281352a0d176a0869e5a8a7a55c0f1d10b1b81d1bd130bb9ffa162f5b8c2db5f027314b7758144f3cb19766658f64d4
    j_invariant = 0x31152b7fe0d9223b7c1d2ea70ba203a46b7120ac12255a04d32cc4d0fbc2c35bf4a01d29a189dbe41969b701f6c5f71e3c77242b402a461c170a7bcf4d55cd76f9968dca9b80049c543b0106ff4df53b3e6293db55f72a8fe560f76bb1d9ee89c171

    #find xi
    xi = 0
    Fq = GF(q)
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1)
    while not (z**6 - xi).is_irreducible():
        xi += 1
    
    test_curve(u0, a, b, q, r, tr, c, y, D, xi)
    
