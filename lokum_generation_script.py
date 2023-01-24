import cocks_pinch as cp
from utils import print_curve
from sage.all import divisors

r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
potential_ks = divisors(r-1)
for potential_k in potential_ks[:-1]:
    if r % potential_k == 1:
        r, potential_k, D = cp.gen_params_from_r(r, potential_k)
        q,t,r,k,D = cp.run(r, potential_k, D)
        print_curve(q,t, r, k, D)
