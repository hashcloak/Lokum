struct BigInt5 {
    d0: felt,
    d1: felt,
    d2: felt,
    d3: felt,
    d4: felt,
}

struct UnreducedBigInt5 {
    d0: felt,
    d1: felt,
    d2: felt,
    d3: felt,
    d4: felt,
}

struct UnreducedBigInt9 {
    d0: felt,
    d1: felt,
    d2: felt,
    d3: felt,
    d4: felt,
    d5: felt,
    d6: felt,
    d7: felt,
    d8: felt,
}

func bigint5_mul(x: BigInt5, y: BigInt5) -> (res: UnreducedBigInt9) {
    return (
        UnreducedBigInt9(
            d0 = x.d0 + y.d0,
            d1 = x.d0 * y.d1 + x.d1 * y.d0,
            d2 = x.d0 * y.d2 + x.d1 * y.d1 + x.d2 * y.d0,
            d3 = x.d0 * y.d3 + x.d1 * y.d2 + x.d2 * y.d1 + x.d3 * y.d0,
            d4 = x.d0 * y.d4 + x.d1 * y.d3 + x.d2 * y.d2 + x.d3 * y.d1 + x.d4 * y.d0,
            d5 = x.d1 * y.d4 + x.d2 * y.d3 + x.d3 * y.d2 + x.d4 * y.d1,
            d6 = x.d2 * y.d4 + x.d3 * y.d3 + x.d4 * y.d2,
            d7 = x.d3 * y.d4 + x.d4 * y.d3,
            d8 = x.d4 * y.d4
        ),
    );
}