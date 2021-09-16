use crate::arith::modular_arithmetic::*;
use rug::{rand::RandState, Assign, Integer};

// TODO: optimize
pub fn eratosthenes(primes: &mut Vec<bool>, limit: usize) {
    primes[0] = false;
    primes[1] = false;
    primes[2] = true;

    let slimit = f64::sqrt(limit as f64) as usize;
    for i in (3..slimit).step_by(2) {
        if primes[i / 2] {
            for j in ((i * i)..limit).step_by(2 * i) {
                primes[j / 2] = false;
            }
        }
    }
}

pub fn bits_amount(n: &Integer) -> u32 {
    let mut amount = 0;
    let mut m = n.clone();
    while m != 0 {
        m >>= 1;
        amount += 1;
    }
    amount
}

pub fn bits_amount_u64(n: u64) -> u32 {
    let mut amount = 0;
    let mut m = n;
    while m != 0 {
        m >>= 1;
        amount += 1;
    }
    amount
}

/// Get bits of an integer from most significant to least significant
pub fn bits(n: &Integer) -> Vec<char> {
    format!("{:b}", n).chars().collect()
}

pub fn randint(rand: &mut RandState, min: &Integer, max: &Integer) -> Integer {
    let min_bits = bits_amount(min);
    let max_bits = bits_amount(max);
    Integer::from(Integer::random_bits(max_bits - min_bits + 1, rand)) + min_bits
}

/// Returns (e, bool) where e is the largest integer such that
/// |y| >= |x^e| and bool=true if y=x^e.
pub fn integer_log(y: u64, x: u64) -> Option<(u32, bool)> {
    match (x, y) {
        (x, y) if (x == 1 || y == 0) => None,
        (x, y) if (0..2).contains(&x) => {
            let e = bits_amount_u64(y) - 1;
            Some((e, x.pow(e) == y))
        }
        _ => {
            let mut r = 0;
            let mut e = 0;
            let mut yy = y;
            while yy >= x {
                let mut d = x;
                let mut m = 1;
                while yy >= d {
                    let div = div_mod(&Integer::from(yy), &Integer::from(d));
                    yy = div.0.to_u64().unwrap();
                    r = r | div.1.to_u64().unwrap();
                    e += m;
                    if yy > d {
                        d *= d;
                        m *= 2;
                    }
                }
            }
            Some((e, r == 0 && yy == 1))
        }
    }
}

pub fn fast_pow(a: &Integer, n: &Integer) -> Integer {
    let mut a_copy = Integer::new();
    a_copy.assign(a);
    let mut n_copy = Integer::new();
    n_copy.assign(n);
    if a_copy == 0 {
        a_copy
    } else {
        let mut acc = Integer::from(1);
        while n_copy > 1 {
            if n_copy.is_even() {
                a_copy.square_mut();
                n_copy >>= 1;
            } else {
                acc = acc * &a_copy;
                a_copy.square_mut();
                n_copy -= 1;
                n_copy >>= 1;
            }
        }
        acc = acc * &a_copy;
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_log_tests() {
        assert_eq!(integer_log(125, 5), Some((3, true)));
        assert_eq!(integer_log(17, 9), Some((1, false)));
    }
}
