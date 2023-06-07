%builtins output range_check

from starkware.cairo.common.uint256 import SHIFT
from src.curve import P0, P1, P2, P3, P4
from src.bigint import BigInt5, UnreducedBigInt9, bigint5_mul
from starkware.cairo.common.registers import get_fp_and_pc
from starkware.cairo.common.serialize import serialize_word

func main{output_ptr: felt*, range_check_ptr}() {
    alloc_locals;
    let (__fp__, _) = get_fp_and_pc();
    //local zero: BigInt3 = BigInt3(0, 0, 0);
    //local Xb: BigInt3 = BigInt3(BASE - 2, BASE - 2, 12);
    //local Yb: BigInt3 = BigInt3(9, 10, 11);
    local zero: BigInt5 = BigInt5(0, 0, 0, 0, 0);
    let p_lokum: BigInt5 = BigInt5(P0, P1, P2, P3, P4);
    let zero_mul_x: UnreducedBigInt9 = bigint5_mul(zero, p_lokum);
    let zero_mul_p: UnreducedBigInt9 = bigint5_mul(p_lokum, p_lokum);

    tempvar x0 = zero_mul_p.d0;
    tempvar x1 = zero_mul_p.d1;
    tempvar x2 = zero_mul_p.d2;
    tempvar x3 = zero_mul_p.d3;
    tempvar x4 = zero_mul_p.d4;
    tempvar x5 = zero_mul_p.d5;
    tempvar x6 = zero_mul_p.d6;
    tempvar x7 = zero_mul_p.d7;
    tempvar x8 = zero_mul_p.d8;

    serialize_word(x0);
    serialize_word(x1);
    serialize_word(x2);
    serialize_word(x3);
    serialize_word(x4);
    serialize_word(x5);
    serialize_word(x6);
    serialize_word(x7);
    serialize_word(x8);
    //let xxx = fq_bigint3.sub(&Xb, &Yb);
    //let res = fq_bigint3.mul(&Xb, &Yb);
    // let res = fq_bigint3.mul(&Xb, &zero);
    // let res = fq_bigint3.mulo(&Xb, &Yb);
    //let res_inv = fq_bigint3.inv(&Yb);
    //let res_inv = fq_bigint3.inv(&Xb);

    //%{ value = 456 + 456*2**86 + 15*2**(86*2) %}
    //let value = nd();
    //let (__fp__, _) = get_fp_and_pc();
    //tempvar y = fp + 1;
    return ();
}
