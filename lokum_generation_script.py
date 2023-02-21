import cocks_pinch as cp
import complex_multiplication as cm
from utils import print_curve
from sage.all import divisors, Integer, log, e, exp

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
        print("Estimated bits of security: ", bits_of_sec)
        print_curve(q,t, r, k, D)
        print(cm.make_curve(q,t,r,k,D))
        print("-------------------------------------------------------------------------------------------------")

