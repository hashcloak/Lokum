from pairing import *

def miller_loop_opt_ate_cp6_782(Q, P, a, u):
    """
    Return frobenius(f_{u+1,Q}(P)) * f_{u*(u^2-u-1),Q}(P)

    Optimized optimal ate Miller loop for CP6_782
    Same formula as for BW6_BLS12 but with non-zero E.a4()
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate((Q[0], -Q[1]), P, a, -u)
    else:
        m_u, uQ = miller_function_ate(Q, P, a, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_j(uQ, (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**2 - u - 1
    m_uv, uvQ = miller_function_ate(uQ, P, a, v, m0=m_u)
    return m_uv * m_u1.frobenius()

def miller_loop_opt_ate_cp6_782_2naf(Q, P, a, u):
    """
    Return frobenius(f_{u+1,Q}(P))*f_{u*(u^2-u-1),Q}(P)

    Optimized optimal ate Miller loop for CP6_782
    Same formula as for BW6_BLS12 but with non-zero E.a4()
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate_2naf((Q[0], -Q[1]), P, a, -u)
    else:
        m_u, uQ = miller_function_ate_2naf(Q, P, a, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_j(uQ, (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**2 - u - 1
    m_uv, uvQ = miller_function_ate_2naf(uQ, P, a, v, m0=m_u)
    return m_uv * m_u1.frobenius()

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

def ate_pairing_cp6_782(Q, P, a, tr):
    T = tr-1
    m,S1 = miller_function_ate(Q, P, a, T)
    f = final_exp_cp6_782(m)
    return f

def ate_pairing_cp6_782_csb(Q, P, a, tr):
    T = tr-1
    m,S1 = miller_function_ate_csb(Q, P, a, T)
    f = final_exp_cp6_782(m)
    return f

def ate_pairing_cp6_782_2naf(Q, P, a, tr):
    T = tr-1
    m,S1 = miller_function_ate_2naf(Q, P, a, T)
    f = final_exp_cp6_782(m)
    return f

def tate_pairing_cp6_782(P, Q, a, r):
    m,S1 = miller_function_tate(P, Q, a, r)
    f = final_exp_cp6_782(m)
    return f

def tate_pairing_cp6_782_csb(P, Q, a, r):
    m,S1 = miller_function_tate_csb(P, Q, a, r)
    f = final_exp_cp6_782(m)
    return f

def tate_pairing_cp6_782_2naf(P, Q, a, r):
    m,S1 = miller_function_tate_2naf(P, Q, a, r)
    f = final_exp_cp6_782(m)
    return f

