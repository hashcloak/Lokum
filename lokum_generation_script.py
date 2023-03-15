import cocks_pinch as cp
import complex_multiplication as cm
from utils import print_curve
from sage.all import divisors, Integer, log, e, exp, euler_phi

def test_finite_field_nfs(q, r, k, cost, sieving_dim=2, samples=100000, special=False, conj=False, sarkarsingh=False, jouxlercier=False, qx=None, u=None, max_coeff=2, deg_f=None, deg_phi_base=None, B0_alpha=800, B1_alpha=2000, compute_alpha=True):
    """
    run NFS for the field GF(q^k) with r a prime divisor of the cyclotomic subgroup
    INPUT:
    - `q`: prime integer, characteristic of the target field
    - `r`: prime integer, order of subgroup of GF(q^k)
    - `k`: extension degree
    - `cost`: expected cost of DL
    - `sieving_dim`: dimension of sieving
    - `samples`: number of samples in the Monte Carlo simulation
    - `special`: Special-NFS
    - `conj`: Conjugation-NFS
    - `sarkarsingh`: Sarkar-Singh-NFS
    - `jouxlercier`: Joux-Lercier-NFS
    - `qx`: polynomial for the special form of q if special is chosen
    - `u`: seed such that q = qx(u) for the special SNFS method
    - `max_coeff`: for Conj or Sarkar-Singh or Joux-Lercier methods
    - `deg_f`: for JL and GJL (jouxlercier) deg_f > k, or Sarkar-Singh: deg_f = deg_phi
    - `deg_phi_base`: for Sarkar-Singh: deg_phi_top*deg_phi_base = k, gcd(aux0, aux1) = phi_top, f=Resultant(aux1, phi_base), g = Resultant(aux0, phi_base) 
    """
    inv_zeta = float(1.0/zeta(float(sieving_dim)))

    poly_init = Polyselect(p=q, k=k)
    # 2. NFS
    if special:
        print("special-NFS")
        res_poly = poly_init.NFS_Special(deg_g=k, poly_p=qx, u=u)
    elif conj:
        print("conj-NFS")
        res_poly = poly_init.NFS_Conjugation(deg_g=k, sieving_dim=sieving_dim, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, verbose=2, number_results=1, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
        # f, g, aux0, aux1, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score
    elif sarkarsingh:
        print("Sarkar-Singh-NFS")
        res_poly = poly_init.NFS_SarkarSingh_JL(deg_phi_top=k//deg_phi_base, deg_aux1=deg_f//deg_phi_base, deg_phi_base=deg_phi_base, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, sieving_dim=sieving_dim, max_test_poly_aux1=100, B1_alpha=B1_alpha)
    elif jouxlercier:
        print("(Generalized)Joux-Lercier-NFS")
        assert deg_f > k
        res_poly = poly_init.GeneralizedJouxLercier(deg_f=deg_f, max_coeff=max_coeff, deg_phi=k, monic=False, sieving_dim=sieving_dim, compute_alpha=compute_alpha, B1_alpha=B1_alpha, max_test_polys=100)

    if res_poly is None or len(res_poly) == 0:
        raise ValueError("Error in Polynomial selection")
    if conj:
        print("res_poly[0] = {}".format(res_poly[0]))
        f, g, aux0, aux1, max_fi, max_gi, aut_g = res_poly[0][:7]
        aut = aut_g
    elif jouxlercier:
        f, g, max_fi, max_gi = res_poly[0][:4]
        aut = 1
    else:
        f, g, max_fi, max_gi, aut_fg = res_poly[:5]
        aut = aut_fg

    print("f = {}".format(f))
    print("g = {}".format(g))
    print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, log(float(max_gi) ,2.0)))
    print("aut = {}".format(aut))

    # computing alpha, this takes at least a few seconds

    if sieving_dim == 2:
        alpha_f = float(alpha2d(f, B1_alpha))
        alpha_g = float(alpha2d(g, B1_alpha))
    elif sieving_dim == 3:
        alpha_f = float(alpha3d(f, B1_alpha))
        alpha_g = float(alpha3d(g, B1_alpha))
    else:
        alpha_f = 0.5
        alpha_g = 0.5
    sum_alpha = alpha_f+alpha_g
    print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f, alpha_g, sum_alpha))
    print("    ({:.8f}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta), float(alpha_f),float(alpha_g),float(sum_alpha)))

    #initialisation of data
    Fp = poly_init.get_Fp()
    Fpz = poly_init.get_Fpz()
    simul = Simulation_NFS(q,r,Fp,Fpz,sieving_dim,f,g,Rx,cost,aut,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

    simul.print_params()
    simul.simulation(samples=samples) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
    simul.print_results()
    print("#::::::::::::::")
    # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
    simul.adjust_cost(samples=samples)
    print("############")

def test_finite_field_tnfs(q, r, k, cost, samples=100000, special=False, conj=False, sarkarsingh=False, jouxlercier=False, all_deg_h=None, qx=None, u=None, max_coeff=2, deg_f=None, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False):
    """
    run TNFS for the field GF(q^k) with r a prime divisor of the cyclotomic subgroup
    INPUT:
    - `q`: prime integer, characteristic of the target field
    - `r`: prime integer, order of subgroup of GF(q^k)
    - `k`: extension degree
    - `cost`: expected cost of DL
    - `samples`: number of samples in the Monte Carlo simulation
    - `special`: Special-TNFS, requires q of special form given with qx, u, such that q = qx(u)
    - `conj`: Conjugation-TNFS
    - `sarkarsingh`: Sarkar-Singh-TNFS, requires k composite
    - `jouxlercier`: Joux-Lercier-TNFS
    - `all_deg_h`: a list of targeted degrees of subfield in TNFS. If k is prime, it is deg_h=k by default.
    - `qx`: polynomial for the special form of q if special is chosen
    - `u`: seed such that q = qx(u) for the special STNFS method
    - `max_coeff`: for Conj or Sarkar-Singh or Joux-Lercier methods
    - `deg_f`: for JL and GJL (jouxlercier)
    """
    poly_init = Polyselect(p=q, k=k)
    # 1. TNFS
    if all_deg_h is None:# test all possible degrees of subfields, so that it divides k
        # if k is prime, the only possible degree is k. deg = 1 corresponds to NFS (below)
        all_deg_h = [i for i in range(k,1,-1) if (k % i) == 0]
    for deg_h in all_deg_h:
        deg_g = k // deg_h
        poly_init.compute_h(deg_h=deg_h)
        list_h = poly_init.get_h()[deg_h]
        print("polynomials h")
        for item in list_h:
            inv_zeta_Kh, w, hc = item
            hc_str = pretty_print_coeffs_from_coeffs(hc)
            h_str = pretty_print_poly_from_coeffs(hc)
            print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
        print("")
        for item in list_h:
            inv_zeta_Kh, w, hc = item
            hc_str = pretty_print_coeffs_from_coeffs(hc)
            h_str = pretty_print_poly_from_coeffs(hc)
            print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
            h = ZZy(hc)

            # actually with D=3, one can have aux = x^2+x+1, and with D=1, aux=x^2+1, but that's all.
            if special:
                print("special-TNFS")
                res_poly = poly_init.TNFS_Special(deg_g=deg_g, h=h, poly_p=px, u=u, with_y=(gcd(deg_g, deg_h)==1))
            elif conj:
                print("conj-TNFS")
                res_poly = poly_init.TNFS_Conjugation(deg_g=deg_g, h=h, with_y=(gcd(deg_g, deg_h) > 1), max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, alpha_test_principal=alpha_test_principal,verbose=2, number_results=1, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
            elif sarkarsingh:
                print("Sarkar-Singh-TNFS")
                res_poly = poly_init.TNFS_SarkarSingh_JL(deg_aux1=deg_aux1, deg_g=deg_g, h=h, with_y=(gcd(deg_g,deg_h) > 1), max_coeff=max_coeff, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
            elif jouxlercier:
                print("Joux-Lercier-TNFS")
                res_poly = poly_init.TNFS_GJL(deg_phi=k//deg_h, deg_f=deg_f, h=h, with_y=False, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, alpha_test_principal=alpha_test_principal, B0_alpha=B0_alpha, B1_alpha=B1_alpha, number_results=1)

            if res_poly == None:
                raise ValueError("Error in Polynomial selection")
            f, g, max_fij, max_gij, aut_fg = res_poly[:5]
            Kh = NumberField(h, names=('ah',)); (ah,) = Kh._first_ngens(1)
            print("h = {} # {}".format(h, h.list()))
            print("inv_zeta_Kh, w = {:.6f},{:2d}".format(inv_zeta_Kh, w))
            assert (ZZ(h.resultant(f.resultant(g))) % q**k) == 0
            if f.parent() is not Rxy:
                f = f(x)
                g = g(x)

            if (gcd(deg_g,deg_h) == 1):
                aut_h = automorphism_factor(hc)
            else:
                aut_h = 1
            aut = aut_h*aut_fg
            print("f = {}".format(f))
            print("g = {}".format(g))
            print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fij, float(log(max_gij,2))))
            print("aut = {}".format(aut))

            # computing alpha, this takes at least a few seconds

            if f.degree() <= 12:# and gcd(deg_g,deg_h) == 1 and gcd(f.degree(),deg_h) == 1:
                alpha_f = float(alpha_TNFS_2d(f,h,1000))
                alpha_g = float(alpha_TNFS_2d(g,h,1000))
            else:
                alpha_f = 0.5
                alpha_g = 0.5
            sum_alpha = alpha_f+alpha_g
            print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f,alpha_g,sum_alpha))
            print("    ({:.8f},{:2d}, {:40s}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta_Kh),int(w),str(h),float(alpha_f),float(alpha_g),float(sum_alpha)))

            #initialisation of data
            Fp = poly_init.get_Fp()
            Fpz = poly_init.get_Fpz()
            simul = Simulation_TNFS(q,r,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

            simul.print_params()
            simul.simulation(samples=samples) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
            simul.print_results()
            print("#::::::::::::::")
            # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
            simul.adjust_cost(samples=samples)
            print("############")

if __name__ == '__main__':
    r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    potential_ks = divisors(r-1)
    for potential_k in potential_ks[:-1]:
        if r % potential_k == 1:
            r, potential_k, D = cp.gen_params_from_r(r, potential_k)
            q,t,r,k,D = cp.run(r, potential_k, D)
            q, r = Integer(q), Integer(r)
            print("q: {}, size of q: {}, r: {}, size of r: {}, k: {}".format(q, q.nbits(), r, r.nbits(), k))
            c, kappa = Integer(32), Integer(-8)
            l = q.nbits() * k * log(2)/log(e)
            bits_of_sec = log((2**(kappa) * exp(pow(c / k, 1/3) * pow(l, 1/3) * pow(log(l), 2 / 3))).n())/log(2).n()
            print("Estimated bits of security (conservative): ", bits_of_sec)
            print("Ideal trace size vs actual trace", (log(r, 2)/euler_phi(k)).n(), t)
            print_curve(q,t, r, k, D)
            #print(cm.make_curve(q,t,r,k,D))
            print("-------------------------------------------------------------------------------------------------")

