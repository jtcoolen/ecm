use rug::ops::Pow;
use rug::Integer;

pub fn div_mod(a: &Integer, modulo: &Integer) -> (Integer, Integer) {
    <(Integer, Integer)>::from(a.div_rem_ref(modulo))
}

pub fn take_mod(a: &Integer, modulo: &Integer) -> Integer {
    // The second element of the tuple is the remainder of the division of a by n
    let rem = div_mod(&a, &modulo).1;
    // Return positive remainder
    if rem < 0 {
        rem + modulo
    } else {
        rem
    }
}

pub fn multiply_mod(a: &Integer, b: &Integer, modulo: &Integer) -> Integer {
    take_mod(&Integer::from(a * b), modulo)
}

pub fn add_mod(a: &Integer, b: &Integer, modulo: &Integer) -> Integer {
    take_mod(&Integer::from(a + b), modulo)
}

pub fn subtract_mod(a: &Integer, b: &Integer, modulo: &Integer) -> Integer {
    take_mod(&Integer::from(a - b), modulo)
}

pub fn invert_mod(a: &Integer, modulo: &Integer) -> Option<Integer> {
    a.invert_ref(modulo).and_then(|b| Some(Integer::from(b)))
}

pub fn pow_mod(a: &Integer, n: u32, modulo: &Integer) -> Integer {
    take_mod(&Integer::from(a).pow(n), modulo)
}
