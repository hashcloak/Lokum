from sage.all import Integer
"""
Pairing : Miller loop and final exponentiation
Formulas for add line and double line
Formulas for generic final exponentiation
and specific formulas for CP6_782 and BW6_761
preprint: https://eprint.iacr.org/2020/351
Optimized and secure pairing-friendly elliptic curves suitable for one layer proof composition
Youssef El Housni and Aurore Guillevic
"""

def double_line_j(S,P,a):
    """ Compute 2*S and l_{S,S}(P) in extended Jacobian coordinates

    Extended Jacobian coordinates (X,Y,Z,Z^2) correspond to the affine
    coordinates (x,y) = (X/Z^2,Y/Z^3)

    INPUT:
    - `S`: elliptic curve point in Jacobian coordinates (X,Y,Z,Z^2)
    - `P`: elliptic curve point in affine coordinates (x,y)
    - `a`: curve coefficient y^2=x^3+a*x+b

    RETURN:
    line tangent at S evaluated at P, point 2*S = (X',Y',Z',Z'^2)

    S and P are on the same elliptic curve (no twist is implemented).
    No vertical line is computed: it is assumed that the embedding
    degree is even and the points have the appropriate representation
    (the x-coordinates are in subfields).

    if S has coordinates in F_{p^k/d} and P has coordinates in F_p,
    it costs:
    a=0:  5*s_{k/d} + 5*m_{k/d} + 2*k/d*m
    a=-3: 4*s_{k/d} + 6*m_{k/d} + 2*k/d*m
    a:    6*s_{k/d} + 5*m_{k/d} + 2*k/d*m
    if S has coordinates in F_p and P has coordinates in F_{p^k/d},
    it costs:
    a=0:  5*s + 6*m + 2*k/d*m
    a=-3: 4*s + 7*m + 2*k/d*m
    a:    6*s + 6*m + 2*k/d*m
    """
    X,Y,Z,Z2 = S
    xP,yP = P
    t1 = Y**2                         # S
    t2 = 4*X*t1                       # M
    if a == 0:
        t3 = 3*X**2                   #   S  (a=0)
    elif a == -3:
        t3 = 3*(X-Z2)*(X+Z2)          #   M  (a=-3)
    else:
        t3 = 3*X**2 + a*Z2**2         #   2S (a any)
    X1= t3**2 - 2*t2                  # S
    Y1 = t3*(t2-X1)-8*t1**2           # M+S
    Z1 = Z*2*Y                        # M
    ld = Z1 * Z2                      # M
    ln = ld * yP - 2*t1 - t3*(Z2*xP-X)# M + 2*k/d*m_1
    # if S is in Fp and P in F_{p^{k/d}}, more efficient to compute:
    #ln = ld*yP-2*t1-(t3*Z2)*xP+t3*X  # 2M + 2*k/d*m_1
    return (ln, (X1,Y1,Z1,Z1**2))     # S

def add_line_j(S,Q,P):
    """ Compute S+Q and l_{S,Q}(P), S in extended Jacobian coordinates

    Jacobian coordinates (X,Y,Z,Z^2) correspond to the affine coordinates
    (x,y) = (X/Z^2,Y/Z^3)

    INPUT:
    - `S`: point in extended Jacobian coordinates (X, Y, Z, Z^2)
    - `Q`: point in affine coordinates (xQ, yQ)
    - `P`: point in affine coordinates (xP, yP)

    RETURN:
    line through S, Q evaluated at P, point S+Q=(X',Y',Z',Z'^2)

    If S,Q have coordinates in F_{p^k/d} and P has coordinates in F_p,
    it costs 10*m_{k/d} + 3*s_{k/d}
    """
    X,Y,Z,Z2 = S
    xQ,yQ = Q
    xP,yP = P
    t1 = xQ*Z2 - X                    #  M
    t2 = yQ*Z*Z2-Y                    # 2M
    t3 = t1**2                        #  S
    t4 = t1*t3                        #  M
    t5 = X*t3                         #  M
    XX = t2**2 - (t4+2*t5)            #  S
    YY = t2*(t5-XX)-Y*t4              # 2M
    ZZ = Z*t1                         #  M
    ld = ZZ
    ln = ld *(yP-yQ) - t2*(xP-xQ)     # 2M
    return (ln, (XX,YY,ZZ,ZZ**2))     #  S

def add_line_j_with_z(S,Q,P):
    """ computes S+Q and l_{S,Q}(P), Q,S in ext Jacobian coordinates

    INPUT:
    - `S`: point in extended Jacobian coordinates (X, Y, Z, Z^2)
    - `Q`: point in extended Jacobian coordinates (X', Y', Z', Z'^2)
    - `P`: point in affine coordinates (xP, yP)

    RETURN: line through S, Q evaluated at P, point S+Q=(X',Y',Z')
    affine coordinates satisfy (x,y) = (X/Z^2,Y/Z^3)
    If S,Q have coordinates in F_{p^k/d} and P has coordinates in F_p, it costs
    14*m_{k/d} + 3*s_{k/d} + 2*k/d*m:
    addition in 11*m_{k/d} + 3*s_{k/d}, evaluation in 3*m_{k/d} + 2*k/d*m.

    This is specifically for the additional terms for optimal ate pairing KSS18
    """
    Xs,Ys,Zs,Zs2 = S
    Xq,Yq,Zq,Zq2 = Q
    xP,yP = P
    # http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    # Zs2 == Z1**2
    # Zq2 == Z2**2
    U1 = Xs*Zq2     # M
    U2 = Xq*Zs2     # M
    Zq3 = Zq*Zq2    # M
    Zs3 = Zs*Zs2    # M
    S1 = Ys*Zq3     # M
    S2 = Yq*Zs3     # M
    # lambda = (S2-S1)/(Xs*Xq*(U2-U1))
    H = U2-U1       #     t1
    I = (2*H)**2    # S   4*t1^2
    J = H*I         # M   4*t1^3
    r = 2*(S2-S1)   #     2*t2
    V = U1*I        # M   Xs*Zq2*4*t1^2
    X3 = r**2-J-2*V # S   2*t2 - 4*t1^3 -2*Xs*Zq2*4*t1^2
    Y3 = r*(V-X3)-2*S1*J # 2M
    Z3 = ((Zq+Zs)**2-Zs2-Zq2)*H # S+M
    ld = Z3
    ln = ld*(yP*Zq3-Yq) - r*Zq*(xP*Zq2 - Xq) # 3 M + 2*k/d*m
    return (ln, (X3,Y3,Z3))

def double_line_j_csb(S,P,a):
    """ Compute 2*S and l_{S,S}(P) in Jacobian coordinates

    INPUT:
    - `S`: elliptic curve point in Jacobian coordinates (X,Y,Z)
    - `P`: elliptic curve point in affine coordinates (x,y)
    - `a`: curve coefficient y^2=x^3+a*x+b

    RETURN:
    line tangent at S evaluated at P, point 2*S = (X',Y',Z',Z'^2)

    Jacobian coordinates (X,Y,Z) correspond to the affine coordinates
    (x,y) = (X/Z^2,Y/Z^3)
    Chatterjee-Sarkar-Barua ICISC'04 formulas.
    10*m_{k/d} + 3*s_{k/d}
    assume S has coordinates in F_{p^k/d} and P has coordinates in F_p

    6m + 5s + k*m if a=0    # 5M + 5S + k*m if a=0
    7m + 4s + k*m if a=-3   # 6M + 4S + k*m if a=-3
    6m + 6s + k*m if a any  # 5M + 6S + k*m if a any
    """
    X1,Y1,Z1 = S
    x,y = P
                                # Tate, S in Fp           # ate, S in Fpk/2
    t1 = Y1**2                  # s                       # S
    t2 = 4*X1*t1                # m                       # M
    t3 = 8*t1**2                # s                       # S
    t4 = Z1**2                  # s                       # S
    if a == 0:
        t5 = 3*X1**2            # s                       # S
    elif a == -3:
        t5 = 3*(X1-t4)*(X1+t4)  # m                       # M
    else:
        t5 = 3*X1**2 + a*t4**2  # 2s                      # 2S
    X3 = t5**2 -2*t2            # s                       # S
    Y3 = t5*(t2-X3)-t3          # m                       # M
    Z3 = 2*Y1*Z1                # m                       # M
    l1 = Z3*t4*y                # m + k/2*m               # M + k/2*m
    #l0= -2*t1 - t5*(t4*x-X1)                             # k/2*m + M
    l0 = -2*t1 - t5*t4*x + t5*X1# 2m + k/2*m
    return (l0+l1,(X3,Y3,Z3))
                                # 6m + 5s + k*m if a=0    # 5M + 5S + k*m if a=0
                                # 7m + 4s + k*m if a=-3   # 6M + 4S + k*m if a=-3
                                # 6m + 6s + k*m if a any  # 5M + 6S + k*m if a any

def add_line_j_csb(S,Q,P):
    """ Compute S+Q and l_{S,Q}(P) in Jacobian coordinates

    INPUT:
    - `S`: elliptic curve point in Jacobian coordinates (X,Y,Z)
    - `Q`: elliptic curve point in affine coordinates (xQ,yQ)
    - `P`: elliptic curve point in affine coordinates (xP,yP)

    RETURN:
    line through S and Q evaluated at P, point S+Q = (X',Y',Z')

    Jacobian coordinates (X,Y,Z) correspond to the affine coordinates
    (x,y) = (X/Z^2,Y/Z^3)
    Chatterjee-Sarkar-Barua ICISC'04 formulas.
    10*m_{k/d} + 3*s_{k/d}
    if S,Q have coordinates in F_{p^k/2} and P has coordinates in F_p
    (ate pairing) it costs 10M + 3S + k/2*m
    if S,Q have coordinates in F_p and P has coordinates in F_{p^k/2}
    (Tate pairing) it costs 9m + k*m + 3s
    """
    X1,Y1,Z1 = S
    xQ,yQ = Q
    xP,yP = P
                                # Tate, S in Fp       # ate, S in F_{p^k/2}
    t1 = Z1**2                  # s                   # S
    t2 = Z1*t1                  # m                   # M
    t3 = xQ*t1                  # m                   # M
    t4 = yQ*t2                  # m                   # M
    t5 = t3-X1
    t6 = t4-Y1
    t7 = t5**2                  # s                   # S
    t8 = t5*t7                  # m                   # M
    t9 = X1*t7                  # m                   # M
    X3 = t6**2 - (t8+2*t9)      # s                   # S
    Y3 = t6*(t9-X3) - Y1*t8     # 2m                  # 2M
    Z3 = Z1*t5                  # m                   # M
    # Tate:
    #l1= Z3*yQ                  # k/2*m               #
    #l0= -Z3*yP - t6*(xQ - xP)  # m + k/2*m           #
    # ate:
    l1 = Z3*yP                  #                     # k/2*m
    l0 = -Z3*yQ - t6*(xP - xQ)  #                     # 2M
    return (l0+l1,(X3,Y3,Z3))   # 9m + k*m + 3s       # 10M + 3S + k/2*m

def double_line_cln_b0(S, P, a):
    """ Compute 2*S and l_{S,S}(P) in weight-(1,2) coordinates

    Assume the curve has coefficient b=0

    INPUT:
    - `S`: elliptic curve point in weight-(1,2) coordinates (X,Y,Z)
    - `P`: elliptic curve point in affine coordinates (x,y)
    - `a`: curve coefficient y^2=x^3+a*x

    RETURN:
    line tangent at S evaluated at P, point 2*S = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + a*x (b=0) Section 4
    Weight-(1,2) coordinates (X,Y,Z) with Y^2 = X^3*Z + a*X*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z^2)

    2*(k/d)*m + 2*m_{k/d} + 8*s_{k/d} + m_a
    assume S has coordinates in F_{p^k/d} and P has coordinates in F_p

    """
    (X1, Y1, Z1) = S
    (xP, yP) = P
    A = X1**2                        # S
    B = Y1**2                        # S
    C = Z1**2                        # S
    D = a*C                          # mult by a
    X3 = (A-D)**2                    # S
    E = 2*(A+D)**2-X3                # S
    F = (A-D+Y1)**2 - B - X3         # S
    Y3 = E*F                         # M
    Z3 = 4*B
    L10 = -2*Z1*(3*A+D)              # M
    L01 = 2*((Y1+Z1)**2 - B - C)     # S
    L00 = (X1 + A - D)**2 - X3 - A   # S
    return (L10*xP + L01*yP + L00, (X3,Y3,Z3)) # 2*k/d*m

def add_line_cln_b0(S, Q, P):
    """ Compute S+Q and l_{S,Q}(P) in weight-(1,2) coordinates

    Assume the curve has coefficient b=0

    INPUT:
    - `S`: elliptic curve point in weight-(1,2) coordinates (X,Y,Z)
    - `Q`: elliptic curve point in affine coordinates (xQ,yQ)
    - `P`: elliptic curve point in affine coordinates (xP,yP)

    RETURN:
    line through S and Q evaluated at P, point S+Q = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + a*x (b=0) Section 4
    Weight-(1,2) coordinates (X,Y,Z) with Y^2 = X^3*Z + a*X*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z^2)

    2*(k/d)*m + 10*m_{k/d} + 5*s_{k/d}
    assume S, Q have coordinates in F_{p^k/d} and P has coordinates in F_p

    addition in 8*m_{k/d} + 5*s_{k/d}
    line evaluation in 2*m_{k/d} + 2*k/d*m
    """
    (X1, Y1, Z1) = S
    (xQ, yQ) = Q
    (xP, yP) = P
    A = Z1**2                             # S
    E = xQ*Z1                             # M
    G = yQ*A                              # M
    H = (X1-E)
    I = 2*(Y1-G)
    II = I**2                             # S
    J = 2*Z1*H                            # M
    K = 4*J*H                             # M
    X3 = 2*II - (X1+E)*K                  # M
    Z3 = J**2                             # S
    Y3 = ((J+I)**2-Z3-II)*(X1*K-X3)-Y1*K**2 # 2S + 3M
    Z3 = 2*Z3
    L10 = -I
    L01 = J
    L00 = I*xQ - J*yQ                     # 2M
    return (L10*xP + L01*yP + L00, (X3,Y3,Z3)) # 2*k/d*m

def add_line_cln_b0_with_z(S, Q, P):
    """ Compute S+Q and l_{S,Q}(P) in weight-(1,2) coordinates

    Assume the curve has coefficient b=0

    INPUT:
    - `S`: elliptic curve point in weight-(1,2) coordinates (X,Y,Z)
    - `Q`: elliptic curve point in weight-(1,2) coordinates (X2,Y2,Z2)
    - `P`: elliptic curve point in affine coordinates (xP,yP)

    RETURN:
    line through S and Q evaluated at P, point S+Q = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + a*x (b=0) Section 4
    Weight-(1,2) coordinates (X,Y,Z) with Y^2 = X^3*Z + a*X*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z^2)

    2*(k/d)*m + 15*m_{k/d} + 7*s_{k/d}
    assume S, Q have coordinates in F_{p^k/d} and P has coordinates in F_p

    addition in 10*m_{k/d} + 7*s_{k/d}
    line evaluation in 5m_{k/d} + 2*k/d*m
    """
    (X1, Y1, Z1) = S
    (X2, Y2, Z2) = Q
    (xP, yP) = P
    A = Z1**2                             # S
    B = Z2**2                             # S
    C = (Z1+Z2)**2-A-B                    # S
    D = X1*Z2                             # M
    E = X2*Z1                             # M
    F = Y1*B                              # M
    G = Y2*A                              # M
    H = (D-E)
    I = 2*(F-G)
    II = I**2                             # S
    J = C*H                               # M
    K = 4*J*H                             # M
    X3 = 2*II - (D+E)*K                   # M
    Z3 = J**2                             # S
    Y3 = ((J+I)**2-Z3-II)*(D*K-X3)-F*K**2 # 2S + 3M
    Z3 = 2*Z3
    L10 = -I*B                            # M
    L01 = J*B                             # M
    L00 = I*X2*Z2 - J*Y2                  # 3M
    return (L10*xP + L01*yP + L00, (X3,Y3,Z3)) # 2*k/d*m

def double_line_cln_a0_cubic_twist(S, P, b):
    """ Compute 2*S and l_{S,S}(P) in projective coordinates

    Assume the curve has coefficient a=0, j=0, and odd embedding degree

    INPUT:
    - `S`: elliptic curve point in projective coordinates (X,Y,Z)
    - `P`: elliptic curve point in extended affine coordinates (x,y,x^2)
    - `b`: curve coefficient for E': y^2=x^3+b where S in E' (so either b or b_t)

    RETURN:
    line tangent at S evaluated at P, point 2*S = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + b (a=0) Section 6
    Projective coordinates (X,Y,Z) with Y^2*Z = X^3 + b*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z)

    3*(k/d)*m + 6*m_{k/d} + 7*s_{k/d} + m_b
    assume S has coordinates in F_{p^k/d} and P has coordinates in F_p
    """
    (X,Y,Z) = S
    (xP,yP,x2P) = P
    A = X**2                     # S
    B = Y**2                     # S
    C = Z**2                     # S
    D = b*C
    E = 3*D
    F = (X+Y)**2-A-B             # S
    G = (Y+Z)**2-B-C             # S
    H = 3*E
    X3 = F * (B-H)               # M
    Y3 = (B+H)**2 -3*(2*E)**2    # 2S
    Z3 = 4*B * G                 # M
    L10 = A * (B-H)              # M
    L20 = F * G                  # M
    L01 = -3*A * F               # M
    L00 = (B-D) * (B+H)          # M
    #assert L10 == X**2 * (Y**2 -9*b*Z**2)
    #assert L20 == 4*X*Y**2*Z
    #assert L01 == -6*X**3*Y
    #assert L00 == (Y**2-b*Z**2)*(Y**2+9*b*Z**2)
    return L10 * xP + L20 * x2P + L01 * yP + L00, (X3, Y3, Z3)

def add_line_cln_a0_cubic_twist(S, Q, P):
    """ Compute S+Q and l_{S,Q}(P) in projective coordinates

    Assume the curve has coefficient a=0, j=0, and odd embedding degree

    INPUT:
    - `S`: elliptic curve point in projective coordinates (X1,Y1,Z1)
    - `Q`: elliptic curve point in affine coordinates (X2,Y2)
    - `P`: elliptic curve point in extended affine coordinates (x,y,x^2)

    RETURN:
    line through S and Q evaluated at P, point S+Q = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + b (a=0) Section 6
    Projective coordinates (X,Y,Z) with Y^2*Z = X^3 + b_t*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z)

    3*(k/d)*m + 13*m_{k/d} + 5*s_{k/d}
    assume S, Q have coordinates in F_{p^k/d} and P has coordinates in F_p

    addition in 11*m_{k/d} + 2*s_{k/d}
    line evaluation in 2*m_{k/d} + 3*s_{k/d} + 3*k/d*m
    """
    (X1,Y1,Z1) = S
    (X2,Y2) = Q
    (xP,yP,x2P) = P
    D = Z1 * X2 - X1            # M
    E = Y1 - Z1 * Y2            # M
    F = D**2                    # S
    G = E**2                    # S
    H = -D * F                  # M
    I = F * X1                  # M
    J = H + Z1 * G - 2*I        # M
    K = Z1 * F * E              # 2M
    X3 = -D * J                 # M
    Y3 = E * (I - J) - (H * Y1) # 2M
    Z3 = Z1 * H                 # M
    L = X3**2                   # S
    M = Z3**2                   # S
    N = (X3 + Z3)**2 - L - M    # S
    L20 = 2*M
    L10 = N
    L00 = 2*(L + K * Y3)        # M
    L01 = -2*K * Z3             # M
    #assert L20 == 2*Z3**2
    #assert L10 == 2*X3*Z3
    #assert L00 == 2*X3**2 + 2*Z1*(X1-X2*Z1)**2*(Y1-Y2*Z1)*Y3
    #assert L01 == -2*Z1*(X1-X2*Z1)**2*(Y1-Y2*Z1)*Z3
    return L20 * x2P + L10 * xP + L01 * yP + L00, (X3, Y3, Z3)

def add_line_cln_a0_cubic_twist_with_z(S, Q, P):
    """ Compute S+Q and l_{S,Q}(P) in projective coordinates

    Assume the curve has coefficient a=0, j=0, and odd embedding degree

    INPUT:
    - `S`: elliptic curve point in projective coordinates (X1,Y1,Z1)
    - `Q`: elliptic curve point in projective coordinates (X2,Y2,Z2)
    - `P`: elliptic curve point in extended affine coordinates (x,y,x^2)

    RETURN:
    line through S and Q evaluated at P, point S+Q = (X',Y',Z')

    Costello Lange Naehrig PKC 2010
    eprint 2009/615 y^2 = x^3 + b (a=0) Section 6
    Projective coordinates (X,Y,Z) with Y^2*Z = X^3 + b_t*Z^3
    correspond to the affine coordinates (x,y) = (X/Z,Y/Z)

    3*(k/d)*m + 16*m_{k/d} + 5*s_{k/d}
    assume S, Q have coordinates in F_{p^k/d} and P has coordinates in F_p

    addition in 14*m_{k/d} + 2*s_{k/d}
    line evaluation in 2*m_{k/d} + 3*s_{k/d} + 3*k/d*m
    """
    (X1,Y1,Z1) = S
    (X2,Y2,Z2) = Q
    (xP,yP,x2P) = P
    A = X1 * Z2                 # M
    B = Y1 * Z2                 # M
    C = Z1 * Z2                 # M
    D = Z1 * X2 - A             # M
    E = B - Z1 * Y2             # M
    F = D**2                    # S
    G = E**2                    # S
    H = -D * F                  # M
    I = F * A                   # M
    J = H + C * G - 2*I         # M
    K = C * F * E               # 2M
    X3 = -D * J                 # M
    Y3 = E * (I - J) - (H * B)  # 2M
    Z3 = C * H                  # M
    L = X3**2                   # S
    M = Z3**2                   # S
    N = (X3 + Z3)**2 - L - M    # S
    L20 = 2*M
    L10 = N
    L00 = 2*(L + K * Y3)        # M
    L01 = -2*K * Z3             # M
    #assert L20 == 2*Z3**2
    #assert L10 == 2*X3*Z3
    #assert L00 == 2*X3**2 + 2*Z1*Z2*(X1*Z2-X2*Z1)**2*(Y1*Z2-Y2*Z1)*Y3
    #assert L01 == -2*Z1*Z2*(X1*Z2-X2*Z1)**2*(Y1*Z2-Y2*Z1)*Z3
    return L20 * x2P + L10 * xP + L01 * yP + L00, (X3, Y3, Z3)

def bits_2naf(x):
    """
    This functions returns the binary non-adjacent form of x. It is
    uniquely defined. The i-th item in the returned list is the
    multiplier attached to 2^i in the 2-NAF.

    sage: bits_2naf(3456780)
    [0, 0, -1, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 1, 0, -1, 0, 1]
    https://gitlab.inria.fr/smasson/cocks-pinch-variant/blob/master/enumerate_sparse_T.py
    paper eprint 2019/431, DOI 10.1007/s10623-020-00727-w
    """
    L = []
    xx = Integer(x)
    assert x >= 0
    while xx > 0 :
        rr = xx % 4
        if rr == 3 :
            rr = -1
        else:
            rr = rr%2
        L.append(rr)
        xx -= rr
        xx,rr = xx.quo_rem(2)
        assert rr == 0
    assert x == sum([r*2**i for i,r in enumerate(L)])
    return L

def miller_function_ate(Q,P,a,T,m0=1):
    """
    computes the Miller function f_{T,Q}(P)

    INPUT:
    - `Q`: r-torsion point on E(Fpk) or E(Fqd) in affine coordinates
    - `P`: r-torsion point on E(Fp) (mapped to E(Fpk) or E(Fqd) in affine coordinates?)
    - `T`: scalar, T=(t-1) for ate pairing for example
    - `a`: curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    - `m0`: optional parameter, for multi-exponentiation optimization,
            this is not needed for simple ate pairing, this is for CP6_782 and BW6_761.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    m = m0
    with_m0 = m0 != 1
    S = (Q[0],Q[1],1,1) # extended Jacobian coordinates
    QQ = (Q[0],Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ = (Q[0],-Q[1])
        S = (Q[0],-Q[1],1,1)
    loop = Integer(T).digits(2)
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_j(S,PP,a)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_j(S,QQ,PP)
            if with_m0:
                m = m*m0
            m = m*ln
    return m, S

def miller_function_tate(P,Q,a,r):
    """
    computes the Miller function f_{r,P}(Q)
    P, Q are r-torsion points in affine coordinates,
    r is an Integer
    a is the curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    """
    return miller_function_ate(P,Q,a,r)

def miller_function_ate_csb(Q,P,a,T,m0=1):
    """
    computes the Miller function f_{T,Q}(P)
    Q,P are r-torsion points in affine coordinates,
    T is a scalar, T=(t-1) for ate pairing for example
    a is the curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    m0 is an optional parameter, for multi-exponentiation optimization,
    this is not needed for simple ate pairing, this is for CP6_782 and BW6_761.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    m = m0
    with_m0 = m0 != 1
    S = (Q[0],Q[1],1)
    QQ = (Q[0],Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ = (Q[0],-Q[1])
        S = (Q[0],-Q[1],1)
    loop = Integer(T).digits(2)
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_j_csb(S,PP,a)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_j_csb(S,QQ,PP)
            if with_m0:
                m = m*m0
            m = m*ln
    return m, S

def miller_function_tate_csb(P,Q,a,r):
    return miller_function_ate_csb(P,Q,a,r)

def miller_function_ate_2naf(Q,P,a,T,m0=1):
    """
    computes the Miller function f_{T,Q}(P)
    Q,P are r-torsion points in affine coordinates,
    T is a scalar, T=(t-1) for ate pairing for example
    a is the curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    m0 is an optional parameter, for multi-exponentiation optimization,
    this is not needed for simple ate pairing, this is for CP6_782 and BW6_761.
    The loop iterates over T in 2-NAF representation.
    If T < 0, to avoid inversion, computes f_{-T, -Q}(P), and the
    byproduct is [-T](-Q)=[T]Q.
    Assumes the embedding degree k is even, do not compute vertical lines.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    S = (Q[0],Q[1],1,1)
    QQ = (Q[0],Q[1])
    negQQ = (Q[0],-Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ, negQQ = negQQ, QQ
        S = (Q[0],-Q[1],1,1)
    loop = bits_2naf(T)
    has_minus_1 = len([i for i in loop if i < 0]) > 0
    m = m0
    with_m0 = m0 != 1
    if with_m0 and has_minus_1:
        m0_inv = 1/m0 # this costs one inversion in Fq^k
    else:
        m0_inv = 1
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_j(S,PP,a)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_j(S,QQ,PP)
            if with_m0:
                m = m*m0
            m = m*ln
        elif bi == -1:
            ln, S = add_line_j(S,negQQ,PP)
            if with_m0:
                m = m*m0_inv
            m = m*ln
    return m, S

def miller_function_tate_2naf(P,Q,a,r):
    return miller_function_ate_2naf(P,Q,a,r)

def miller_function_ate_cln_b0(Q, P, a, T, m0=1):
    """ Miller loop with Costello-Lauter-Naehrig formulas """
    m = m0
    with_m0 = m0 != 1
    S = (Q[0],Q[1],1)
    QQ = (Q[0],Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ = (Q[0],-Q[1])
        S = (Q[0],-Q[1],1)
    loop = Integer(T).digits(2)
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_cln_b0(S,PP,a)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_cln_b0(S,QQ,PP)
            if with_m0:
                m = m*m0
            m = m*ln
    return m, S

def miller_function_tate_cln_b0(P,Q,a,r):
    """
    computes the Miller function f_{r,P}(Q)
    P, Q are r-torsion points in affine coordinates,
    r is an Integer
    a is the curve coefficient in y^2=x^3+a*x (short Weierstrass, b=0)
    """
    return miller_function_ate_cln_b0(P,Q,a,r)

def miller_function_ate_2naf_cln_b0(Q,P,a,T,m0=1):
    """
    computes the Miller function f_{T,Q}(P)
    Q,P are r-torsion points in affine coordinates,
    T is a scalar, T=(t-1) for ate pairing for example
    a is the curve coefficient in y^2=x^3+a*x (short Weierstrass, b=0)
    m0 is an optional parameter, for multi-exponentiation optimization,
    this is not needed for simple ate pairing.
    The loop iterates over T in 2-NAF representation.
    If T < 0, to avoid inversion, computes f_{-T, -Q}(P), and the
    byproduct is [-T](-Q)=[T]Q.
    Assumes the embedding degree k is even, do not compute vertical lines.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    S = (Q[0],Q[1],1)
    QQ = (Q[0],Q[1])
    negQQ = (Q[0],-Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ, negQQ = negQQ, QQ
        S = (Q[0],-Q[1],1)
    loop = bits_2naf(T)
    has_minus_1 = len([i for i in loop if i < 0]) > 0
    m = m0
    with_m0 = m0 != 1
    if with_m0 and has_minus_1:
        m0_inv = 1/m0 # this costs one inversion in Fq^k
    else:
        m0_inv = 1
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_cln_b0(S,PP,a)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_cln_b0(S,QQ,PP)
            if with_m0:
                m = m*m0
            m = m*ln
        elif bi == -1:
            ln, S = add_line_cln_b0(S,negQQ,PP)
            if with_m0:
                m = m*m0_inv
            m = m*ln
    return m, S

def miller_function_tate_2naf_cln_b0(P,Q,a,r):
    return miller_function_ate_2naf_cln_b0(P,Q,a,r)

### for curves of odd embedding degree multiple of 3 with cubic twist

def miller_function_ate_cln_a0_cubic_twist(Q, P, b_t, T, m0=1):
    """ Miller loop with Costello-Lauter-Naehrig formulas """
    m = m0
    with_m0 = m0 != 1
    S = (Q[0], Q[1], 1)
    QQ = (Q[0], Q[1])
    PP = (P[0], P[1], P[0]**2)
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ = (Q[0], -Q[1])
        S = (Q[0], -Q[1], 1)
        # consider the inverse of the vertical line at Q evaluated at P
        #m = m * (P[0] - Q[0])
        inv_vertical_Q_P = PP[2] + PP[0]*QQ[0] + QQ[0]**2
        m = m * inv_vertical_Q_P
    loop = Integer(T).digits(2)
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_cln_a0_cubic_twist(S, PP, b_t)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_cln_a0_cubic_twist(S, QQ, PP)
            if negative_T:
                m = m * inv_vertical_Q_P
            if with_m0:
                m = m*m0
            m = m*ln
    return m, S

def miller_function_tate_cln_a0_cubic_twist(P, Q, b, r, m0=1):
    """
    computes the Miller function f_{r,P}(Q)
    P, Q are r-torsion points in affine coordinates,
    r is an Integer
    b is the curve coefficient in y^2=x^3+b (short Weierstrass, a=0)
    """
    if (r % 2) == 0:
        R = r // 2
        b0 = 0
    else:
        R = r-1
        b0 = 1
    f, S = miller_function_ate_cln_a0_cubic_twist(P, Q, b, R, m0)
    # final step: we will get possibly 0
    # if b0 == 1: S+P = 0, the line through S and P is a vertical at P
    QQ = (Q[0], Q[1])
    PP = (P[0], P[1], P[0]**2)
    if b0 == 1:
        f = f * (S[0]-Q[0]*S[2])
        ln, S = add_line_cln_a0_cubic_twist(S,QQ,PP)
    else:
        f = f**2 * (S[0]-Q[0]*S[2])
        ln, S = double_line_cln_a0_cubic_twist(S, PP, b_t)
    return f, S

def miller_function_ate_2naf_cln_a0_cubic_twist(Q, P, b_t, T, m0=1):
    """
    computes the Miller function f_{T,Q}(P)
    Q,P are r-torsion points in affine coordinates,
    T is a scalar, T=(t-1) for ate pairing for example
    b_t is the curve coefficient in y^2=x^3+b_t (short Weierstrass, a=0)
    m0 is an optional parameter, for multi-exponentiation optimization,
    this is not needed for simple ate pairing.
    The loop iterates over T in 2-NAF representation.
    If T < 0, to avoid inversion, computes f_{-T, -Q}(P), and the
    byproduct is [-T](-Q)=[T]Q.
    Assumes the embedding degree k is odd, compute vertical lines.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is computed as TODO, hence
    f_{-|T|,Q} = f_{|T|,-Q} * TODO
    """
    S = (Q[0], Q[1], 1)
    QQ = (Q[0], Q[1])
    negQQ = (Q[0], -Q[1])
    PP = (P[0], P[1], P[0]**2)
    negative_T = (T < 0)
    m = m0
    if negative_T:
        T = -T
        QQ, negQQ = negQQ, QQ
        S = (Q[0], -Q[1], 1)
    loop = bits_2naf(T)
    has_minus_1 = len([i for i in loop if i < 0]) > 0
    if negative_T or has_minus_1 :
        inv_vertical_Q_P = PP[2] + PP[0]*QQ[0] + QQ[0]**2
    if negative_T:
        # consider the vertical line at Q evaluated at P
        #m = m / (P[0] - Q[0])
        m = m * inv_vertical_Q_P
    with_m0 = m0 != 1
    if with_m0 and has_minus_1:
        m0_inv = 1/m0 # this costs one inversion in Fq^k
    else:
        m0_inv = 1
    for i in range(len(loop)-2, -1, -1):
        bi = loop[i]
        ln, S = double_line_cln_a0_cubic_twist(S, PP, b_t)
        m = m**2 * ln
        if bi == 1:
            ln, S = add_line_cln_a0_cubic_twist(S, QQ, PP)
            if negative_T:
                m = m * inv_vertical_Q_P
            if with_m0:
                m = m*m0
            m = m*ln
        elif bi == -1:
            ln, S = add_line_cln_a0_cubic_twist(S, negQQ, PP)
            if not negative_T:
                m = m * inv_vertical_Q_P
            if with_m0:
                m = m*m0_inv
            m = m*ln
    return m, S

def miller_function_tate_2naf_cln_a0_cubic_twist(P, Q, b_t, r):
    return miller_function_ate_2naf_cln_a0_cubic_twist(P, Q, b_t, r)

###

def psi_sextic_m_twist(Q, alpha):
    """
    Returns the 6-th D-twist of a point Q in affine coordinates
    """
    if len(Q) == 2:
        return (Q[0]/alpha**2, Q[1]/alpha**3)
    if len(Q) == 3:
        return (Q[0]/alpha**2, Q[1]/alpha**3, Q[2])
    if len(Q) == 4:
        return (Q[0]/alpha**2, Q[1]/alpha**3, Q[2], Q[3])

def psi_sextic_d_twist(Q, alpha):
    """
    Returns the 6-th M-twist of a point Q in affine coordinates
    """
    if len(Q) == 2:
        return (Q[0]*alpha**2, Q[1]*alpha**3)
    if len(Q) == 3:
        return (Q[0]*alpha**2, Q[1]*alpha**3, Q[2])
    if len(Q) == 4:
        return (Q[0]*alpha**2, Q[1]*alpha**3, Q[2], Q[3])

#### AKLGL paper Eurocrypt 2011 ####
def double_line_h_a0_twist6_aklgl(S,P,b_t,D_twist=False):
    """
    Computes 2*S and l_{S,S}(P) in Homogeneous coordinates (X,Y,Z) and (x,y) = (X/Z,Y/Z)

    b' is the curve parameter for E', the sextic twist of E.
    E': y^2 = x^3 + b' = x^3 + b/xi and E' is a D-twist of E
    xi is the non-residue s.t. Fp6 = Fq[w]/(w^6-xi)

    For BLS12 curves, xi in Fp2
    For BW6_761 curve, xi in Fp.
    Faster Explicit Formulas for Computing Pairings over Ordinary Curves
    Catherine H. Gebotys, Koray Karabina, Patrick Longa, Julio Lopez, Diego F. Aranha
    http://www.iacr.org/archive/eurocrypt2011/66320047/66320047.pdf
    https://eprint.iacr.org/2010/526 Section 4 [Aranha Karabina Longa Gebotys Lopez Eurocrypt'11]
    Fp2 = Fp[i]/(i^2-beta), beta = -1
    Fp6 = Fp2[v]/(v^3-xi), xi = (1+i)
    Fp12 = Fp6[w]/(w^2 - v)
    E:  y^2 = x^3 + 2
    xi = (1+i)
    E': y^2 = x^3 + 2/(1+i) = x^3 + (1-i) -> b' = 1-i and E' is a D-twist
    """
    X1,Y1,Z1 = S
    xP,yP = P
    # Homogeneous coordinates: Z*Y^2 = X^3 + b'*Z^3 <=> (Y/Z)^2 = (X/Z)^3 + b'
    # note that Y^2 - b'*Z^2 = X^3/Z
    # D-twist: y^2 = x^3 + b/xi (multiply by xi=w^6) => w^6*y^2 = w^6*x^3 + b
    #                                                  <=> (w^3*y)^2 = (w^2*x)^3+b
    # psi: (x',y') -> (x'*w^2,y'*w^3) <-> (x',y'*w,1/w^2)
    # so that Y^2*Z is in Fq (no w^i term)
    #X3 = X1*Y1/2 * (Y1**2-9*b'*Z1**2)
    #Y3 = ((Y1**2+9*b'*Z1**2)/2)**2 - 27*b'**2*Z1**4
    #Z3 = 2*Y1**3*Z1
    #ln = ((-2*Y1*Z1*yP)*v*w + (3*X1**2*xP)*v**2 + xi*(3*b'*Z1**2-Y1**2))
    # with v=w^2, xi = w^6, this is
    #ln = xi*(3*b'*Z1**2-Y1**2) + (-2*Y1*Z1*yP)*w^3 + (3*X1**2*xP)*w**4

    A = X1*Y1/2            # M
    B = Y1**2              # S
    C = Z1**2              # S
    D = 3*C                #     3*Z1^2
    E = b_t*D              #     3*b_t*Z1^2
    F = 3*E                #     9*b_t*Z1^2
    X3 = A*(B-F)           # M   A*(B-9*b_t*C) = X1*Y1/2*(Y1^2-9*b_t*Z1^2)
    G = (B+F)/2            #     (Y1^2 + 9*b_t*Z1^2)/2
    Y3 = G**2-3*E**2       # 2S  (Y1^2 + 9*b_t*Z1^2)^2/4 -3*(3*b_t*Z1^2)^2
    H = (Y1+Z1)**2 - (B+C) # S   2*Y1*Z1
    Z3 = B*H               # M   2*Y1^3*Z1
    I = E-B                #     3*b_t*Z1^2-Y1^2
    J = X1**2              # S
    l3 = I
    l0 = H*(-yP)       # k/d*m
    l1 = 3*J*xP        # k/d*m
    if D_twist:
        ln = [l0,l1,0,l3,0,0]
        # AKLGL: multiply line by w^3 i.e. rotate by 3 to the right the vector -> now we have l3*w^6 = l3*xi at position 0
        #ln = (l3*xi,0,0,l0,l1,0)
        # and (l3*xi,0,l1) + (0,l0,0) as (Fq3)^2
        #     (l3*xi,0,xP*ll1) + (0,-yP*ll0,0)
    else:
        ln = [l3,0,l1,l0,0,0]
    return ln,(X3,Y3,Z3)   # 3M + 6S + 2*k/d*m

def double_line_h_a0_twist6_aklgl_no_div2(S,P,b_t,D_twist=False):
    """
    Computes 2*S and l_{S,S}(P) in Homogeneous coordinates (X,Y,Z) and (x,y) = (X/Z,Y/Z)

    b' is the curve parameter for E', the sextic twist of E.
    E': y^2 = x^3 + b' = x^3 + b/xi and E' is a D-twist of E or
    E': y^2 = x^3 + b' = x^3 + b*xi and E' is a M-twist of E

    Faster Explicit Formulas for Computing Pairings over Ordinary Curves
    Catherine H. Gebotys, Koray Karabina, Patrick Longa, Julio Lopez, Diego F. Aranha
    http://www.iacr.org/archive/eurocrypt2011/66320047/66320047.pdf
    https://eprint.iacr.org/2010/526 Section 4 [Aranha Karabina Longa Gebotys Lopez]
    https://eprint.iacr.org/2010/526 Section 4 [Aranha Karabina Longa Gebotys Lopez Eurocrypt'11]
    """
    X1,Y1,Z1 = S
    xP,yP = P
    # Homogeneous coordinates: Z*Y^2 = X^3 + b'*Z^3 <=> (Y/Z)^2 = (X/Z)^3 + b'

    A = X1*Y1              # M
    B = Y1**2              # S
    C = Z1**2              # S
    D = 3*C                #     3*Z1^2
    E = b_t*D              #     3*b_t*Z1^2
    F = 3*E                #     9*b_t*Z1^2
    X3 = 2*A*(B-F)         # M   2*A*(B-9*b_t*C) = 2*X1*Y1*(Y1^2-9*b_t*Z1^2)
    G = (B+F)              #     (Y1^2 + 9*b_t*Z1^2)
    Y3 = G**2-3*(2*E)**2   # 2S  (Y1^2 + 9*b_t*Z1^2)^2 -3*4*(3*b_t*Z1^2)^2
    H = (Y1+Z1)**2 - (B+C) # S   2*Y1*Z1
    Z3 = 4*B*H             # M   2*Y1^3*Z1
    I = E-B                #     3*b_t*Z1^2-Y1^2
    J = X1**2              # S
    l3 = I
    l0 = H*(-yP)       # k/d*m
    l1 = 3*J*xP        # k/d*m
    if D_twist:
        ln = [l0,l1,0,l3,0,0]
        # AKLGL: multiply line by w^3 i.e. rotate by 3 to the right the vector -> now we have l3*w^6 = l3*xi at position 0
        #ln = (l3*xi,0,0,l0,l1,0)
        # and (l3*xi,0,l1) + (0,l0,0) as (Fq3)^2
        #     (l3*xi,0,xP*ll1) + (0,-yP*ll0,0)
    else:
        ln = [l3,0,l1,l0,0,0]
    return ln,(X3,Y3,Z3)   # 3M + 6S + 2*k/d*m

def double_line_aklgl_test(S,P,b_t,w,D_twist=False):
    X1,Y1,Z1 = S
    xP,yP = P
    X3 = (X1*Y1)/2 * (Y1**2-9*b_t*Z1**2)
    Y3 = ((Y1**2+9*b_t*Z1**2)/2)**2 - 27*b_t**2*Z1**4
    Z3 = 2*Y1**3*Z1
    if D_twist:
        l = (-2*Y1*Z1*yP) + (3*X1**2*xP)*w + (3*b_t*Z1**2-Y1**2)*w**3
        lxy = ((3*b_t*Z1**2-Y1**2)*w**3, (3*X1**2)*w, (-2*Y1*Z1))
    else:
        l = (-2*Y1*Z1*yP)*w**3 + (3*X1**2*xP)*w**2 + (3*b_t*Z1**2-Y1**2)
        lxy = ((3*b_t*Z1**2-Y1**2), (3*X1**2)*w**2, (-2*Y1*Z1)*w**3)
    return l, (X3,Y3,Z3), lxy

def add_line_h_a0_twist6_aklgl(S,Q,P,D_twist=False):
    """ computes S+Q and l_{S,Q}(P), Q,S in Homogeneous coordinates, P affine

    INPUT:
    - `S`: point in Homogeneous coordinates (X, Y, Z)
    - `Q`: point in affine coordinates (xQ, yQ)
    - `P`: point in affine coordinates (xP, yP)

    RETURN: line through S, Q evaluated at P, point S+Q=(X',Y',Z')
    affine coordinates satisfy (x,y) = (X/Z,Y/Z)
    If S,Q have coordinates in F_{p^k/d} and P has coordinates in F_p, it costs
    11*m_{k/d} + 2*s_{k/d} + 2*k/d*m_1
    Algorithm 12 in eprint 2010/526
    """
    (X1,Y1,Z1) = S
    (X2,Y2) = Q
    (xP,yP) = P
    t1 = Z1*X2 ; t2 = Z1*Y2     # 2M
    t1 = X1-t1 ; t2 = Y1-t2
    # D = X1 - Z1*X2 -> t1
    # E = Y1 - Y2*Z1 -> t2
    t3 = t1**2                  # S
    # F = D^2 -> t1^2 -> t3
    X3 = t3*X1 ; t4 = t2**2     # M+S
    # I = X1*F -> X3
    # G = E^2 -> t2^2 -> t4
    t3 = t1*t3 ; t4 = t4*Z1     # 2M
    # H = D*F -> t1*t3 -> t3
    # J = H + Z1*G - 2*I -> t3 + t4 -2*X3
    t4 = t4+t3
    t4 = t4-X3
    t4 = t4-X3 # J -> t4
    X3 = X3-t4
    T1 = t2*X3 ; T2 = t3*Y1     # 2M
    T2 = T1 - T2
    Y3 = T2; X3=t1*t4; Z3=t3*Z1 # 2M
    lx = -t2*xP                 # k/d*m  lx=-(Y1-Z1*Y2)*xP
    T1 = t2*X2                  # M
    T2 = t1*Y2                  # M
    l0 = T1 - T2                #        l0=(Y1-Z1*Y2)*X2-(X1-Z1*X2)*Y2
    ly = t1*yP                  # k/d*m  ly=(X1-Z1*X2)*yP
    if D_twist:
        ln = [ly,lx,0,l0,0,0] # *w^3 -> (l0*xi,0,0,ly,lx,0)
    else:
        ln = [l0,0,lx,ly,0,0]
    return ln,(X3,Y3,Z3) # 11M + 2S + 2*k/d*m

def add_line_h_a0_twist6_aklgl_test(S,Q,P,w,D_twist=False):
    """
    Fp2 = Fp[i]/(i^2-beta), beta = -1
    Fp6 = Fp2[v]/(v^3-xi), xi = (1+i)
    Fp12 = Fp6[w]/(w^2 - v)
    """
    (X1,Y1,Z1) = S
    (X2,Y2) = Q
    (xP,yP) = P
    T = Y1-Y2*Z1
    L = X1 - X2*Z1
    X3 = L*(L**3+Z1*T**2-2*X1*L**2)
    Y3 = T*(3*X1*L**2 - L**3 - Z1*T**2) - Y1*L**3
    Z3 = Z1*L**3
    if D_twist:
        l = L*(yP) - (T*xP)*w + (T*X2-L*Y2)*w**3
        lxy = ((T*X2-L*Y2)*w**3, -T*w, L)
    else:
        l = (T*X2-L*Y2) - (T*xP)*w**2 + L*(yP)*w**3
        lxy = ((T*X2-L*Y2), -T*w**2, L*w**3)
    return l, (X3,Y3,Z3), lxy

def add_line_h_a0_twist6_aklgl_with_z(S,Q,P,D_twist=False):
    """ computes S+Q and l_{S,Q}(P), Q,S in Homogeneous coordinates, P affine

    INPUT:
    - `S`: point in Homogeneous coordinates (X, Y, Z)
    - `Q`: point in Homogeneous coordinates (xQ, yQ, zQ)
    - `P`: point in affine coordinates (xP, yP)

    RETURN: line through S, Q evaluated at P, point S+Q=(X',Y',Z')
    affine coordinates satisfy (x,y) = (X/Z,Y/Z)
    If S,Q have coordinates in F_{p^k/d} and P has coordinates in F_p, it costs
    16*m_{k/d} + 2*s_{k/d} + 2*k/d*m_1
    This is specifically for the additional terms for optimal ate pairing KSS18.
    Algorithm 12 in eprint 2010/526 + Z2 != 1
    """
    (X1,Y1,Z1) = S
    (X2,Y2,Z2) = Q
    (xP,yP) = P
    yP_ = -yP
    t1 = Z1*X2 ; t2 = Z1*Y2       # 2M
    X1Z2 = X1*Z2                  # M
    Y1Z2 = Y1*Z2                  # M
    t1 = X1Z2-t1 ; t2 = Y1Z2-t2
    # D = X1*Z2 - Z1*X2 -> t1 (lambda)
    # E = Y1*Z2 - Y2*Z1 -> t2 (theta)
    t3 = t1**2                    # S
    # F = D^2 -> t1^2 -> t3
    X3 = t3*X1Z2 ; t4 = t2**2     # M+S
    # I = X1*F -> X3
    # G = E^2 -> t2^2 -> t4
    Z1Z2 = Z1*Z2                  # M
    t3 = t1*t3 ; t4 = t4*Z1Z2     # 2M
    # H = D*F -> t1*t3 -> t3
    # J = H + Z1*G - 2*I -> t3 + t4 -2*X3
    t4 = t4+t3
    t4 = t4-X3
    t4 = t4-X3 # J -> t4
    X3 = X3-t4
    T1 = t2*X3 ; T2 = t3*Y1Z2     # 2M
    T2 = T1 - T2
    Y3 = T2; X3=t1*t4; Z3=t3*Z1Z2 # 2M
    lx = -t2*Z2*xP              # M + k/d*m  lx=-Z2*(Y1*Z2-Z1*Y2)*xP
    T1 = t2*X2                  # M
    T2 = t1*Y2                  # M
    l0 = T1 - T2                #        l0=(Y1*Z2-Z1*Y2)*X2-(X1*Z2-Z1*X2)*Y2
    ly = t1*Z2*yP               # M + k/d*m  ly=Z2*(X1*Z2-Z1*X2)*yP
    if D_twist:
        ln = [ly,lx,0,l0,0,0] # *w^3 -> (l0*xi,0,0,ly,lx,0)
    else:
        ln = [l0,0,lx,ly,0,0]
    return ln,(X3,Y3,Z3) # 16M + 2S + 2*k/d*m

def sparse_mult_M6_twist(l0,l2,l3, f, xi, Fq6):
    # source: PhD thesis A. Guillevic 2013 p. 91 Sect. 3.2.2
    # https://tel.archives-ouvertes.fr/tel-00921940
    #print("f= {}\nf.polynomial()={}\nf.polynomial().list()={}\n".format(f, f.polynomial(), f.polynomial().list()))
    # problem of type in Python sometimes here
    try:
        coeffs = f.polynomial().list()
    except AttributeError as err:
        coeffs = f.list()
    if len(coeffs) < 6:
        coeffs += [0]*(6-len(coeffs))
    (f0,f1,f2,f3,f4,f5) = coeffs
    l0f0 = l0*f0
    l2f2 = l2*f2
    l3f5 = l3*f5
    h2 = xi*l3f5 + (l0+l2)*(f0+f2) - l0f0 - l2f2
    l2f4 = l2*f4
    l3f3 = l3*f3
    h0 = l0f0 + xi*(l2f4 + l3f3)
    l2f1 = l2*f1
    h3 = l2f1 + (l0+l3)*(f0+f3) - l0f0 - l3f3
    l0f4 = l0*f4
    l3f1 = l3*f1
    h4 = l2f2 + l0f4 + l3f1
    l0f5 = l0*f5
    h5 = l0f5 + (l2+l3)*(f2+f3) - l2f2 - l3f3
    h1 = (l0+l2+l3)*(f1 + xi*(f4+f5)) - xi*(l0f4 + l0f5 + l2f4 + l3f5) - l2f1 - l3f1
    return Fq6([h0,h1,h2,h3,h4,h5])

def sparse_mult_D6_twist(l0,l1,l3, f, xi, Fq6):
    # source: PhD thesis A. Guillevic 2013 p. 91 Sect. 3.2.2
    # https://tel.archives-ouvertes.fr/tel-00921940
    #print("f= {}\nf.polynomial()={}\nf.polynomial().list()={}\n".format(f, f.polynomial(), f.polynomial().list()))
    try:
        coeffs = f.polynomial().list()
    except AttributeError as err:
        coeffs = f.list()
    if len(coeffs) < 6:
        coeffs += [0]*(6-len(coeffs))
    (f0,f1,f2,f3,f4,f5) = coeffs
    l1f5 = l1*f5
    l0f0 = l0*f0
    l3f3 = l3*f3
    h0 = l0f0 + xi*(l1f5 + l3f3)
    l1f1 = l1*f1
    l3f4 = l3*f4
    h1 = (l0+l1)*(f0+f1) - l0f0 - l1f1 + xi*l3f4
    l0f2 = l0*f2
    l3f5 = l3*f5
    h2 = l0f2 + l1f1 + xi*l3f5
    l1f2 = l1*f2
    h3 = l1f2 + (l0+l3)*(f0+f3) - l0f0 - l3f3
    l0f4 = l0*f4
    h4 = l0f4 + (l1+l3)*(f1+f3) - l1f1 - l3f3
    h5 = (l0+l1+l3)*(f2+f4+f5) - (l0f4 + l1f2 + l1f5 + l3f4 + l0f2 + l3f5)
    return Fq6([h0,h1,h2,h3,h4,h5])

def sparse_sparse_mult_D6_twist(l0, l1, l3, h0, h1, h3, xi, Fq6):
    """
    sparse-sparse multiplication in Fq6

    INPUT:
    l0, l1, l3 are the coefficients of l = l0 + l1*a + l3*a^3 mod a^6-xi
    h0, h1, h3 are the coefficients of h = h0 + h1*a + h3*a^3 mod a^6-xi

    RETURN:
    l*h = l0*h0 + l3*h3*xi + (l0*h1 + l1*h0)*a + (l1*h1)*a^2
          + (l0*h3 + l3*h0)*a^3 + (l1*h3 + l3*h1)*a^4

    Total cost: 6*m_{k/6} + 1*mult_by_xi
    """
    c0 = l0*h0
    c6 = l3*h3
    r0 = c0 + c6*xi
    r2 = l1*h1
    r1 = (l0+l1)*(h0+h1) - r2 - c0
    r3 = (l0+l3)*(h0+h3) - c0 - c6
    r4 = (l1+l3)*(h1+h3) - r2 - c6
    r5 = 0
    return Fq6([r0, r1, r2, r3, r4, r5])

def sparse_sparse_mult_M6_twist(l0, l2, l3, h0, h2, h3, xi, Fq6):
    """
    sparse-sparse multiplication in Fq6

    INPUT:
    l0, l2, l3 are the coefficients of l = l0 + l2*a^2 + l3*a^3 mod a^6-xi
    h0, h2, h3 are the coefficients of h = h0 + h2*a^2 + h3*a^3 mod a^6-xi

    RETURN:
    l*h = l0*h0 + l3*h3*xi + (l0*h2 + l2*h0)*a^2 + (l0*h3 + l3*h0)*a^3
          + (l2*h2)*a^4 + (l2*h3 + l3*h2)*a^5

    Total cost: 6*m_{k/6} + 1*mult_by_xi
    """
    c0 = l0*h0
    c6 = l3*h3
    r0 = c0 + c6*xi
    r1 = 0
    r4 = l2*h2
    r2 = (l0+l2)*(h0+h2) - r4 - c0
    r3 = (l0+l3)*(h0+h3) - c0 - c6
    r5 = (l2+l3)*(h2+h3) - r4 - c6
    return Fq6([r0, r1, r2, r3, r4, r5])

def miller_function_ate_aklgl(Q,P,b_t,T,Fq6,D_twist=False,m0=1,xi=None):
    """
    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    m = m0
    with_m0 = m0 != 1
    S = (Q[0],Q[1],1)
    QQ = (Q[0],Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ = (Q[0],-Q[1])
        S = (Q[0],-Q[1],1)
    loop = Integer(T).digits(2)
    # very first step: no "m = m**2" needed
    bi = loop[len(loop)-2]
    ln, S = double_line_h_a0_twist6_aklgl(S,PP,b_t,D_twist=D_twist)
    if with_m0:
        m = m**2
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
    else:
        m = Fq6(ln)
    if bi == 1:
        if with_m0:
            m = m*m0
        ln, S = add_line_h_a0_twist6_aklgl(S,QQ,PP,D_twist=D_twist)
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)

    for i in range(len(loop)-3, -1, -1):
        bi = loop[i]
        ln, S = double_line_h_a0_twist6_aklgl(S,PP,b_t,D_twist=D_twist)
        m = m**2
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
        if bi == 1:
            if with_m0:
                m = m*m0
            ln, S = add_line_h_a0_twist6_aklgl(S,QQ,PP,D_twist=D_twist)
            if D_twist:
                l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
                m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
            else:
                l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
                m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
    return m, S

def miller_function_ate_2naf_aklgl(Q,P,b_t,T,Fq6,D_twist=False,m0=1,xi=None):
    """
    computes the Miller function f_{T,Q}(P) with AKLGL Dbl/Add formulas
    Q,P are r-torsion points in affine coordinates,
    T is a scalar, T=(t-1) for ate pairing for example
    b_t is the 6-twist curve coefficient in y^2=x^3+b_t (short Weierstrass)
    m0 is an optional parameter, for multi-exponentiation optimization,
    this is not needed for simple ate pairing, this is for CP6_782 and BW6_761.
    The loop iterates over T in 2-NAF representation.

    If T < 0, then f_{|T|, -Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|T|,Q} = f_{|T|,-Q}.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    S = (Q[0],Q[1],1)
    QQ = (Q[0],Q[1])
    negQQ = (Q[0],-Q[1])
    PP = (P[0],P[1])
    negative_T = (T < 0)
    if negative_T:
        T = -T
        QQ, negQQ = negQQ, QQ
        S = (Q[0],-Q[1],1)
    loop = bits_2naf(T)
    has_minus_1 = len([i for i in loop if i < 0]) > 0
    m = m0
    with_m0 = m0 != 1
    if with_m0 and has_minus_1:
        m0_inv = 1/m0 # this costs one inversion in Fq^k
    else:
        m0_inv = 1
    # very first step: no "m = m**2" needed
    bi = loop[len(loop)-2]
    ln, S = double_line_h_a0_twist6_aklgl(S,PP,b_t,D_twist=D_twist)
    if with_m0:
        m = m**2
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
    else:
        m = Fq6(ln)
    if bi == 1:
        if with_m0:
            m = m*m0
        ln, S = add_line_h_a0_twist6_aklgl(S,QQ,PP,D_twist=D_twist)
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
    elif bi == -1:
        if with_m0:
            m = m*m0_inv
        ln, S = add_line_h_a0_twist6_aklgl(S,negQQ,PP,D_twist=D_twist)
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)

    for i in range(len(loop)-3, -1, -1):
        bi = loop[i]
        ln, S = double_line_h_a0_twist6_aklgl(S,PP,b_t,D_twist=D_twist)
        m = m**2
        if D_twist:
            l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
            m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
        else:
            l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
            m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
        if bi == 1:
            if with_m0:
                m = m*m0
            ln, S = add_line_h_a0_twist6_aklgl(S,QQ,PP,D_twist=D_twist)
            if D_twist:
                l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
                m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
            else:
                l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
                m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
        elif bi == -1:
            if with_m0:
                m = m*m0_inv
            ln, S = add_line_h_a0_twist6_aklgl(S,negQQ,PP,D_twist=D_twist)
            if D_twist:
                l0 = ln[0] ; l1 = ln[1] ; l3 = ln[3]
                m = sparse_mult_D6_twist(l0,l1,l3, m, xi, Fq6)
            else:
                l0 = ln[0] ; l2 = ln[2] ; l3 = ln[3]
                m = sparse_mult_M6_twist(l0,l2,l3, m, xi, Fq6)
    return m, S

#### Exponentiations ####

# https://www.ams.org/journals/mcom/2013-82-281/S0025-5718-2012-02625-1/S0025-5718-2012-02625-1.pdf
# Koray Karabina, Squaring in Cyclotomic Subgroups
# Mathematics of Computation Volume 82, Number 281, January 2013, Pages 555–579S 0025-5718(2012)02625-1
# Article electronically published on June 27, 2012


def multi_exp(f0, e0, f1, e1):
    # returns f_0^e0 * f_1^e1
    # assumes e0 > 0, e1 > 1
    if e0 < e1:
        ee = e0
        e0 = e1
        e1 = ee
        ff = f0
        f0 = f1
        f1 = ff
    e0 = Integer(e0).digits(2)
    e1 = Integer(e1).digits(2)
    f = f0
    S = 0
    M = 0
    for i in range(len(e0)-2, len(e1)-2, -1): # for i= # e0-1 to # e1 by -1 do
        f = f**2;         S += 1
        bi = e0[i]
        if bi == 1:
            f = f * f0;  M += 1
    f = f * f1;          M += 1
    f0f1 = f0*f1;        M += 1
    for i in range(len(e1)-2, -1, -1): # i= # e1-1 to 1 by -1 do
        f = f**2;         S += 1
        bi = e0[i]
        ci = e1[i]
        if bi == 1 and ci == 0:
            f *= f0;     M += 1
        elif bi == 0 and ci == 1:
            f *= f1;     M += 1
        elif bi == 1 and ci == 1:
            f *= f0f1;   M += 1
    print("multi-exp: {}S + {}M".format(S,M))
    return f

def exp_2naf(f, f_inv, e_2naf):
    # returns f^e_2naf
    # assumes e_2naf > 0
    g = f
    for i in range(len(e_2naf)-2, -1, -1):
        g = g**2
        bi = e_2naf[i]
        if bi == 1:
            g *= f
        elif bi == -1:
            g *= f_inv
    return g

def multi_exp_2naf(f0, f0_inv, f1, f1_inv, e0_2naf, e1_2naf, verbose=False):
    # returns f_0^e0_2naf * f_1^e1_2naf
    # assumes e0 > 0, e1 > 1
    f = f0
    S = 0
    M = 0
    for i in range(len(e0_2naf)-2, len(e1_2naf)-2, -1): # i= # e0_2naf-1 to # e1_2naf by -1 do
        f = f**2;               S += 1
        bi = e0_2naf[i]
        if bi == 1:
            f *= f0;           M += 1
        elif bi == -1:
            f *= f0_inv;       M += 1
    # precomputations
    f = f * f1;                M += 1
    f0f1 = f0*f1;              M += 1
    f0if1 = f0_inv*f1;         M += 1
    f0f1i = f0*f1_inv;         M += 1
    f0if1i = f0_inv*f1_inv;    M += 1
    for i in range(len(e1_2naf)-2, -1, -1): # i= # e1_2naf-1 to 1 by -1 do
        f = f**2;               S += 1
        bi = e0_2naf[i]
        ci = e1_2naf[i]
        if bi == 1 and ci == 0:
            f *= f0;           M += 1
        elif bi == 0 and ci == 1:
            f *= f1;           M += 1
        elif bi == 1 and ci == 1:
            f *= f0f1;         M += 1
        elif bi == -1 and ci == 0:
            f *= f0_inv;       M += 1
        elif bi == 0 and ci == -1:
            f *= f1_inv;       M += 1
        elif bi == 1 and ci == -1:
            f *= f0f1i;        M += 1
        elif bi == -1 and ci == 1:
            f *= f0if1;        M += 1
        elif bi == -1 and ci == -1:
            f *= f0if1i;       M += 1
    if verbose:
        print("multi-exp: {}S + {}M".format(S,M))
    return f, S,M

def final_exp_easy_k6(m):
    """cost: f3 + i6 + 2*m6 + f
    """
    mp3 = m.frobenius(3) # assuming the extension is in one layer
    im = 1/m
    f = mp3*im # m^(q^3-1)
    f = f * f.frobenius() # m^((q^3-1)*(q+1))
    return f

def final_exp_easy_k9(m):
    # embedding degree k=9, easy part is (p^9-1)/Phi_9(p) = (p^3-1)
    """cost: f3 + i9 + m9
    """
    mp3 = m.frobenius(3) # assuming the extension is in one layer
    im = 1/m
    f = mp3*im # m^(q^3-1)
    return f

def final_exp_easy_k12(m):
    # embedding degree k=12, easy part is (p^12-1)/Phi_12(p) = (p^6-1)*(p^2+1)
    mp6 = m.frobenius(6)
    im = 1/m
    f = mp6*im # m^(q^6-1)
    f = f * f.frobenius(2) # m^((q^6-1)*(q^2+1))
    return f

def final_exp_easy_k15(m):
    # embedding degree k=15, easy part is (p^15-1)/Phi_15(p) = (p^5-1)*(p^2+p+1)
    # cost: f5 + i15 + 3*m15
    mp5 = m.frobenius(5)
    im = 1/m
    f = mp5*im # m^(q^5-1)
    f = f * f.frobenius() * f.frobenius(2) # m^((q^5-1)*(q^2+q+1))
    return f

def final_exp_easy_k16(m):
    # embedding degree k=16, easy part is (p^16-1)/Phi_16(p) = (p^8-1)
    mp8 = m.frobenius(8)
    im = 1/m
    f = mp8*im # m^(q^8-1)
    return f

def final_exp_easy_k18(m):
    # embedding degree k=18, easy part is (p^18-1)/Phi_18(p) = (p^9-1)*(p^3+1)
    mp9 = m.frobenius(9)
    im = 1/m
    f = mp9*im # m^(q^9-1)
    f = f * f.frobenius(3) # m^((q^9-1)*(q^3+1))
    return f

def final_exp_easy_k24(m):
    # embedding degree k=24, easy part is (p^24-1)/Phi_24(p) = (p^12-1)*(p^4+1)
    mp12 = m.frobenius(12)
    im = 1/m
    f = mp12*im # m^(q^12-1)
    f = f * f.frobenius(4) # m^((q^12-1)*(q^4+1))
    return f

def final_exp_easy_k27(m):
    # embedding degree k=27, easy part is (p^27-1)/Phi_27(p) = (p^9-1)
    # cost: f9 + i27 + m27
    mp9 = m.frobenius(9) # assuming the extension is in one layer
    im = 1/m
    f = mp9*im # m^(q^9-1)
    return f

# On final exponentation, see
# - Michael Scott, Naomi Benger, Manuel Charlemagne, Luis J. Dominguez Perez, and Ezekiel J. Kachisa.
#   On the final exponentiation for calculating pairings on ordinary elliptic curves.
#   In Pairing-Based Cryptography -- Pairing 2009, LNCS 5671, pp. 78--88, Springer, https://doi.org/10.1007/978-3-642-03298-1\_6, 2009.
# - Castaneda et al work
# - Efficient Final Exponentiation via Cyclotomic Structure for Pairings over Families of Elliptic Curves
#   Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
#   https://eprint.iacr.org/2020/875

def final_exp_hard_bls9(m, u):
    """
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    page 13
    exponent = (u-1)^2*(q^2+u*q+u^2)*(q^3+u^3+1) + 3
             = (u-1)^2*(q^2+u*(q+u))*(q^3+u^3+1) + 3
    cost 2*exp(u-1) + 5*exp(u) + 6M + S + 3f
    """
    m0 = m**((u-1)**2)                           # 2 exp(u-1)
    m1 = m0**u * m0.frobenius()                  # exp(u) + M + f2
    m1 = m1**u * m0.frobenius(2)                 # exp(u) + M + f
    m1 = m1.frobenius(3) * m1**(u**3) * m1       # 3*exp(u) + 2M + f
    return m1 * m**2 * m                         # 2 M + S

def final_exp_hard_bls15(m, u):
    """
    exponent*3 = (u-1)^2*(u^2+u+1)*(q*(q^6+q^3+q+1) + (u-1)*(q^6+u*q^5+u^2*q^4+(u^3+1)*q^3+u*(u^3+1)*q^2+(u^2*(u^3+1)+1)*q+u^3*(u^3+1)+u+1)) + 3
    cost 3*exp(u-1) + 8*exp(u) + 17 M + S + 10 f

    alternatively, note that (u-1)^2*(u^2+u+1) = (u-1)*(u^3-1)
    cost exp(u-1) + 3*exp(u) + M + inv_cyclo instead of 2*exp(u-1) + 2*exp(u) + 2*M
    but inv_cyclo costs one M.
    """
    m0 = m**((u-1)**2)               # 2 exp(u-1)
    m1 = m0**u                       # exp(u)
    m1 = (m1 * m0)**u                # exp(u) + M
    m1 = m1 * m0                     # M
    m1q = m1.frobenius()             # f
    m2 = m1q * m1q.frobenius() * m1q.frobenius(3) * m1q.frobenius(6) # 3M + 3f
    m3 = m1**(u-1)                   # exp(u-1)
    mu = m3**u                       # exp(u)
    mb = m3.frobenius(6) * mu.frobenius(5) # M + 2f
    mu = mu**u                       # exp(u)
    mb = mb * mu.frobenius(4)        # M + f
    mu = mu**u * m3                  # exp(u) + M
    mb = mb * mu.frobenius(3)        # M + f
    mu = mu**u                       # exp(u)
    mb = mb * mu.frobenius(2)        # M + f
    mu = mu**u * m3                  # exp(u) + M
    mb = mb * mu.frobenius()         # M + f
    mu = mu**u                       # exp(u)
    mb = mb * mu * m3                # 2 M
    return m**2 * m * m2 * mb        # 3 M + S

def final_exp_hard_bls12_ghammam_fouotsa(m,u):
    # formulas from GF18
    # Loubna Ghammam and Emmanuel Fouotsa,
    # Improving the computation of the optimal ate pairing for a high security level.
    # J. Appl. Math. Comput. 59, 21-36 (2019). https://doi.org/10.1007/s12190-018-1167-y
    # http://eprint.iacr.org/2016/130
    # exponent = 3*(p^4-p^2+1)/r
    # l3 = (u-1)**2
    # l2 = l3*u
    # l1 = l2*u-l3
    # l0 = l1*u+3
    # l0+px*(l1+px*(l2+px*l3)) == 3*exponent
    # cost 2*exp(u-1) + 3*exp(u) + 6M + S + 3f + f6
    m3 = m**(u-1)
    m3 = m3**(u-1)
    m2 = m3**u
    res = m3.frobenius() * m2
    m1 = m2**u * m3.frobenius(6) # m3^(-1) = m3^(p^6)
    res = res.frobenius() * m1
    m0 = m1**u
    m0 = m0 * m**2 * m
    res = res.frobenius() * m0
    return res

def final_exp_hard_bls12(m, u):
    """
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    page 14

    exponent = (u-1)^2/3 * (q + x) * (q^2 + u^2 - 1) + 1
    3*exponent = (u-1)^2 * (q + x) * (q^2 + u^2 - 1) + 3
    cost 2*exp(u-1) + 3*exp(u) + 5 M + S + 2f + f6
    """
    m1 = m**(u-1)
    m1 = m1**(u-1)
    m2 = m1**u
    m1 = m1.frobenius() * m2
    m1 = m1.frobenius(2) * (m1**u)**u * m1.frobenius(6)
    return m1 * m**2 * m

def final_exp_hard_2naf_bls12(m,u):
    # TODO with eprint 2020/875 formula
    if u < 0: # (u-1)^2 == (-u+1)^2
        e_u1 = bits_2naf(-u+1)
    else:
        e_u1 = bits_2naf(u-1)
    m_inv = m.frobenius(6)
    m3 = exp_2naf(m, m_inv, e_u1)
    m_inv = m3.frobenius(6)
    m3 = exp_2naf(m3, m_inv, e_u1)
    if u > 0:
        e_u = bits_2naf(u)
    else:
        e_u = bits_2naf(-u)
    m3_inv = m3.frobenius(6)
    m2 = exp_2naf(m3, m3_inv, e_u)
    m2_inv = m2.frobenius(6)
    if u < 0:
        m2, m2_inv = m2_inv, m2
    res = m3.frobenius() * m2

    if u < 0:
        m1_inv = exp_2naf(m2, m2_inv, e_u)*m3
        m1 = (m1_inv).frobenius(6) # m1^(-1) = m1^(p^6)
    else:
        m1 = exp_2naf(m2, m2_inv, e_u)*m3_inv # m3^(-1) = m3^(p^6)
        m1_inv = m1.frobenius(6)
    res = res.frobenius() * m1

    m0 = exp_2naf(m1, m1_inv, e_u)
    if u < 0:
        m0 = m0.frobenius(6)
    m0 *= m**2 * m
    res = res.frobenius() * m0
    return res

def final_exp_hard_bls24(m, u):
    """
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    page 15

    (p^8-p^4+1)/r = (u-1)^2*(u+p)*(u^2+p^2)*(u^4+p^4-1) + 3
    cost 2 exp(u-1) + 7*exp(u) + 6 M + S + f + f2 + f4 + f12
    """
    m1 = m**(u-1)                                         # exp(u-1)
    m1 = m1**(u-1)                                        # exp(u-1)
    m1 = m1.frobenius() * m1**u                           # exp(u) + M + f
    m1 = m1.frobenius(2) * (m1**u)**u                     # 2 exp(u) + M + f2
    m1 = m1**(u**4) * m1.frobenius(4) * m1.frobenius(12)  # 4 exp(u) + 2 M + f4 + f12
    return m1 * m**2 * m                                  # 2 M + S

def final_exp_hard_bls24_ghammam_fouotsa(m,u):
    # formulas from GF18
    # Loubna Ghammam and Emmanuel Fouotsa,
    # Improving the computation of the optimal ate pairing for a high security level.
    # J. Appl. Math. Comput. 59, 21-36 (2019). https://doi.org/10.1007/s12190-018-1167-y
    # http://eprint.iacr.org/2016/130
    #l7 = x**2-2*x+1
    #l6 = l7*x
    #l5 = l6*x
    #l4 = l5*x
    #l3 = l4*x-l7
    #l2 = l3*x
    #l1 = l2*x
    #l0 = l1*x + 3
    #l0+px*(l1+px*(l2+px*(l3+px*(l4+px*(l5+px*(l6+px*l7)))))) == 3*exponent
    # u**2-2*u+1 = (u-1)^2
    m7 = m**(u-1)
    m7 = m7**(u-1)
    m6 = m7**u
    res = m7.frobenius() * m6
    m5 = m6**u
    res = res.frobenius() * m5
    m4 = m5**u
    res = res.frobenius() * m4
    m3 = m4**u * m7.frobenius(12) # m7^(-1) = m7^(p^12)
    res = res.frobenius() * m3
    m2 = m3**u
    res = res.frobenius() * m2
    m1 = m2**u
    res = res.frobenius() * m1
    m0 = m1**u
    m0 *= m**2 * m
    res = res.frobenius() * m0
    return res

def final_exp_hard_2naf_bls24(m,u):
    #l7 = x**2-2*x+1
    #l6 = l7*x
    #l5 = l6*x
    #l4 = l5*x
    #l3 = l4*x-l7
    #l2 = l3*x
    #l1 = l2*x
    #l0 = l1*x + 3
    #l0+px*(l1+px*(l2+px*(l3+px*(l4+px*(l5+px*(l6+px*l7)))))) == 3*exponent
    # u**2-2*u+1 = (u-1)^2
    if u < 0: # (u-1)^2 == (-u+1)^2
        e_u1 = bits_2naf(-u+1)
    else:
        e_u1 = bits_2naf(u-1)
    m_inv = m.frobenius(12)
    m7 = exp_2naf(m, m_inv, e_u1)
    m_inv = m7.frobenius(12)
    m7 = exp_2naf(m7, m_inv, e_u1)
    if u > 0:
        e_u = bits_2naf(u)
    else:
        e_u = bits_2naf(-u)
    m7_inv = m7.frobenius(12)
    m6 = exp_2naf(m7, m7_inv, e_u)
    m6_inv = m6.frobenius(12)
    if u < 0:
        m6, m6_inv = m6_inv, m6
    res = m7.frobenius() * m6
    m5 = exp_2naf(m6, m6_inv, e_u)
    m5_inv = m5.frobenius(12)
    if u < 0:
        m5, m5_inv = m5_inv, m5
    res = res.frobenius() * m5
    m4 = exp_2naf(m5, m5_inv, e_u)
    m4_inv = m4.frobenius(12)
    if u < 0:
        m4, m4_inv = m4_inv, m4
    res = res.frobenius() * m4
    if u < 0:
        m3_inv = exp_2naf(m4, m4_inv, e_u)*m7
        m3 = (m3_inv).frobenius(12) # m7^(-1) = m7^(p^12)
    else:
        m3 = exp_2naf(m4, m4_inv, e_u)*m7_inv # m7^(-1) = m7^(p^12)
        m3_inv = m3.frobenius(12)
    res = res.frobenius() * m3
    m2 = exp_2naf(m3, m3_inv, e_u)
    m2_inv = m2.frobenius(12)
    if u < 0:
        m2, m2_inv = m2_inv, m2
    res = res.frobenius() * m2
    m1 = exp_2naf(m2, m2_inv, e_u)
    m1_inv = m1.frobenius(12)
    if u < 0:
        m1, m1_inv = m1_inv, m1
    res = res.frobenius() * m1
    m0 = exp_2naf(m1, m1_inv, e_u)
    if u < 0:
        m0 = m0.frobenius(12)
    m0 *= m**2 * m
    res = res.frobenius() * m0
    return res

def final_exp_hard_bls27_zhang_lin(m, u):
    """
    Zhang Lin INDOCRYPT'2012 Section 4.3
    Analysis of Optimum Pairing Products at High Security Levels
    (p^18+p^9+1)/r = 3 + (u-1)^2 * (p^9 + u^9 + 1)
    * (p^8 +u*p^7 + u^2*p^6 + u^3*p^5 +u^4*p^4 + u^5*p^3 + u^6*p^2 + u^7*p + u^8)
    cost 17 exp(u) + 12 M + S + f9 + f8 + f7 + f6 + f5 + f4 + f3 + f2 + f
    """
    m1 = m.frobenius(8)                        # f8
    m_ui = m
    for i in range(1, 8):                      # 7*exp(u) + 7 M + f7 + f6 + f5 + f4 + f3 + f2 + f
        m_ui = m_ui**u
        m1 = m1 * m_ui.frobenius(8-i)
    m_ui = m_ui**u # m**(u**8)                 # exp(u)
    m1 = m1 * m_ui                             # M
    m2 = m1**((u-1)**2)                        # 2 exp(u-1)
    m3 = m2**(u**9) * m2.frobenius(9) * m2     # 9*exp(u) + 2 M + f9
    m4 = (m**2 * m) * m3                       # 2 M + S
    return m4

def final_exp_hard_bls27(m, u):
    """
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya

    (p^18+p^9+1)/r = 3+(u-1)^2*(p^2+p*u+u^2)*(p^6+p^3*u^3+u^6)*(p^9+u^9+1)
                   = 3+(u-1)^2*(p^2+(p+u)*u)*(p^6+(p^3+u^3)*u^3)*(p^9+u^9+1)
    cost 2 exp(u-1) + 17*exp(u) + 8 M + S + f + f2 + f3 + f6 + f9
    """
    m1 = (m**(u-1))**(u-1)                  # 2 exp(u-1)
    m2 = m1**u * m1.frobenius()             # exp(u) + M + f
    m2 = m2**u * m1.frobenius(2)            # exp(u) + M + f2
    m1 = m2**(u**3) * m2.frobenius(3)       # 3*exp(u) + M + f3
    m1 = m1**(u**3) * m2.frobenius(6)       # 3*exp(u) + M + f6
    m2 = m1 * m1**(u**9) * m1.frobenius(9)  # 9*exp(u) + 2 M + f9
    return m2 * m**2 * m                    # 2 M + S

def final_exp_hard_bls48(m, u):
    """
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    page 17

    (p^16-p^8+1)/r = (u-1)^2*(u+p)*(u^2+p^2)*(u^4+p^4)(u^8+p^8-1) + 3
    cost 2 exp(u-1) + 15*exp(u) + 7 M + S + f + f2 + f4 + f8 + f24
    """
    m1 = m**(u-1)                                         # exp(u-1)
    m1 = m1**(u-1)                                        # exp(u-1)
    m1 = m1.frobenius() * m1**u                           # exp(u) + M + f
    m1 = m1.frobenius(2) * (m1**u)**u                     # 2 exp(u) + M + f2
    m1 = m1**(u**4) * m1.frobenius(4)                     # 4 exp(u) + M + f4
    m1 = m1**(u**8) * m1.frobenius(8) * m1.frobenius(24)  # 8 exp(u) + 2 M + f8 + f24
    return m1 * m**2 * m                                  # 2 M + S

def final_exp_bls12(m,u):
    f = final_exp_easy_k12(m)
    g = final_exp_hard_bls12(f,u)
    return g

def final_exp_2naf_bls12(m,u):
    f = final_exp_easy_k12(m)
    g = final_exp_hard_2naf_bls12(f,u)
    return g

def final_exp_bls24(m,u):
    f = final_exp_easy_k24(m)
    g = final_exp_hard_bls24(f,u)
    return g

def final_exp_2naf_bls24(m,u):
    f = final_exp_easy_k24(m)
    g = final_exp_hard_2naf_bls24(f,u)
    return g

def ate_pairing_bls12_aklgl(Q,P,b_t,u0,Fq6,map_Fp12_Fp12_A,D_twist=False):
    m,S = miller_function_ate_aklgl(Q,P,b_t,u0,Fq6,D_twist=D_twist,m0=1)
    # convert m from tower field to absolute field
    m = map_Fp12_Fp12_A(m)
    f = final_exp_bls12(m,u0)
    return f

def ate_pairing_bls12_2naf_aklgl(Q,P,b_t,u0,Fq6,map_Fp12_Fp12_A,D_twist=False):
    m,S = miller_function_ate_2naf_aklgl(Q,P,b_t,u0,Fq6,D_twist=D_twist,m0=1)
    # convert m from tower field to absolute field
    m = map_Fp12_Fp12_A(m)
    f = final_exp_2naf_bls12(m,u0)
    return f

def ate_pairing_bls24_aklgl(Q,P,b_t,u0,Fq6,map_Fp24_Fp24_A,D_twist=False):
    m,S = miller_function_ate_aklgl(Q,P,b_t,u0,Fq6,D_twist=D_twist,m0=1)
    # convert m from tower field to absolute field
    m = map_Fp24_Fp24_A(m)
    f = final_exp_bls24(m,u0)
    return f

def ate_pairing_bls24_2naf_aklgl(Q,P,b_t,u0,Fq6,map_Fp24_Fp24_A,D_twist=False):
    m,S = miller_function_ate_2naf_aklgl(Q,P,b_t,u0,Fq6,D_twist=D_twist,m0=1)
    # convert m from tower field to absolute field
    m = map_Fp24_Fp24_A(m)
    f = final_exp_2naf_bls24(m,u0)
    return f

####### membership testing, co-factor multiplication
# Formulas from Fuentes-Castaneda, Knapp, Rodriguez-Henriquez
# Faster hashing to G2.
# In: Miri, Vaudenay (eds.) SAC 2011. LNCS, vol. 7118, pp. 412-430.
# Springer, Heidelberg (Aug 2012). https://doi.org/10.1007/978-3-642-28496-0_25
# Example for BLS12-381: https://eprint.iacr.org/2019/814
# Sean Bowe, Faster Subgroup Checks for BLS12-381

def bw6_phi(P, omega):
    """ Return Phi(P(x,y)) = (omega*x, y) where omega^2+omega+1=0 mod p
    """
    Fp = P[0].parent()
    return P.curve()([Fp(omega)*P[0], P[1]])

