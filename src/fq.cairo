from starkware.cairo.common.uint256 import SHIFT
from starkware.cairo.common.cairo_builtins import BitwiseBuiltin
from starkware.cairo.commoon.registers import get_fp_and_pc
from src.curve import P0, P1, P2, P3, P4, N_LIMBS, N_LIMBS_UNREDUCED, DEGREE, BASE
from src.bigint import BigInt5, UnreducedBigInt5, UnreducedBigInt9, bigint5_mul

const BASE_MIN_1 = BASE - 1;

func fq_zero() -> BigInt5 {
    let res = BigInt5(0, 0, 0, 0, 0);
    return res;
}

func fq_eq_zero(x: BigInt5*) -> felt {
    if (x.d0 != 0) {
        return 0;
    }
    if (x.d1 != 0) {
        return 0;
    }
    if (x.d2 != 0) {
        return 0;
    }
    if (x.d3 != 0) {
        return 0;
    }
    if (x.d4 != 0) {
        return 0;
    }
    return 1;
}

func fq_eq_one(x: BigInt5*) -> felt {
    if (x.d0 != 1) {
        return 0;
    }
    if (x.d1 != 0) {
        return 0;
    }
    if (x.d2 != 0) {
        return 0;
    }
    if (x.d3 != 0) {
        return 0;
    }
    if (x.d4 != 0) {
        return 0;
    }
    return 1;
}