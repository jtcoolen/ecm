pub mod arith;
use crate::arith::misc::*;
use crate::arith::modular_arithmetic::*;
use crate::arith::montgomery_point::MontgomeryPoint;
use log::{debug, info};
use rug::{rand::RandState, Assign, Integer};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;

/// Lenstra's Elliptic Curve Method for Factorization (ECM).
/// Returns a nontrivial factor of n upon success.

/// Notes:
/// The bound b1 (resp. b2) for stage 1 (resp. stage 2) must be even,
/// and is usually taken s.t. b2 ~ 100*b1.

/// The boolean found_factor is shared by all threads and set to false initially.

/// Implements Algorithm 7.4.4 (Inversionless ECM) from the book
/// Prime Numbers from R. Crandall and C. B. Pomerance.
pub fn inversionless_ecm(
    n: &Integer,
    max_curves: &Option<Integer>,
    primes: &Vec<bool>,
    b1: u64,
    b2: u64,
    sigma: &Option<Integer>,
    thread_no: usize,
    found_factor: &AtomicBool,
) -> Option<Integer> {
    debug!("max_curves={:?}", max_curves);
    debug!("B1={}", b1);
    debug!("B2={}", b2);

    let mut rand = RandState::new();
    rand.seed(&Integer::from(thread_no));

    let d: usize = (b2 as f64).sqrt() as usize;

    let mut points = Vec::new();
    let mut beta = Vec::new();
    for _ in 0..(d + 1) {
        points.push(MontgomeryPoint::default());
        beta.push(Integer::default());
    }
    let mut curve = Integer::from(0);

    let mut infinite = false;
    let mut limit = Integer::from(0);
    match max_curves {
        Some(l) => limit = Integer::from(l),
        None => infinite = true,
    };

    while curve < limit || infinite {
        curve += 1;
        // We check before trying a new curve if an other
        // thread has already found a factor, in which case
        // we return
        if found_factor.load(Ordering::Relaxed) {
            return None;
        }
        info!("Curve {}", curve);

        // Choose a random curve using Suyama's parametrization
        // If sigma is provided, only one iteration is needed since
        // there's only one curve to try out
        let sigma = match sigma {
            Some(s) => {
                // only one iteration
                infinite = false;
                limit = Integer::from(0);
                Integer::from(s)
            }
            None => randint(&mut rand, &Integer::from(6), &Integer::from(n - 1)),
        };
        debug!("Sigma={}", sigma);
        let v = multiply_mod(&Integer::from(4), &sigma, n);
        let u = subtract_mod(&Integer::from(&sigma).square(), &Integer::from(5), n);

        let diff = subtract_mod(&v, &u, n);
        let u_cubed = pow_mod(&u, 3, n);

        match invert_mod(&(4 * Integer::from(&u_cubed * &v)), n) {
            None => {
                let a: Integer = 4 * Integer::from(&u_cubed * &v);
                info!("Sigma={}", sigma);
                debug!("found factor\n\n\n");
                found_factor.swap(true, Ordering::Relaxed);
                // if a is not invertible mod n then by Bezout the GCD of a and n is > 1
                return Some(a.gcd(n));
            }
            Some(inv) => {
                // c determines the curve in Montgomery form y^2 = x^3 + cx^2 + x
                let c = take_mod(
                    &Integer::from((pow_mod(&diff, 3, n) * Integer::from(3 * &u + &v) * &inv) - 2),
                    n,
                );
                // Initial point in Montgomery form [X:Z]=[u^3 mod n : v^3 mod n]
                let mut q = MontgomeryPoint::new2(
                    Integer::from(&u_cubed),
                    pow_mod(&v, 3, &n),
                    Integer::from(&c),
                    Integer::from(n),
                );

                // Stage 1
                info!("Stage 1");
                let mut k = Integer::from(1);
                for p_i in 2..(b1 + 1) {
                    if primes[p_i as usize] {
                        // will fail if b1 is bigger than a usize
                        match integer_log(b1, p_i) {
                            // find largest integer a s.t. p_i^a is <= to our first bound b1
                            Some(a) => {
                                // Compute Q = [p_i^a] Q using Montgomery's ladder algo
                                // TODO: Maybe implement some sort of FFT?
                                k *= fast_pow(&Integer::from(p_i), &Integer::from(a.0));
                            }
                            None => return None,
                        }
                    }
                }
                q = q.montgomery_ladder(&k);
                let mut g = Integer::from(&q.z).gcd(n);

                if 1 < g && g < *n {
                    info!("Sigma={}", sigma);
                    debug!("found factor {}\n\n\n", g);
                    found_factor.swap(true, Ordering::Relaxed);
                    return Some(g);
                }

                // Stage 2
                info!("Stage 2");
                points[1] = q.double();
                points[2] = points[1].double();
                beta[1] = multiply_mod(&points[1].x, &points[1].z, n);
                beta[2] = multiply_mod(&points[2].x, &points[2].z, n);

                // Compute points[idx] = 2*idx.q
                for idx in 3..(d + 1) {
                    points[idx] = points[idx - 1].addh(&points[1], &points[idx - 2]);
                    // Keep the products X*Z
                    beta[idx] = multiply_mod(&points[idx].x, &points[idx].z, n);
                }

                g.assign(1);
                let b = b1 - 1;
                let mut t = q.montgomery_ladder(&Integer::from(b - 2 * (d as u64)));
                let mut s = q.montgomery_ladder(&Integer::from(b));

                for r in (b..b2).step_by(2 * d) {
                    let alpha = take_mod(&Integer::from(&s.x * &s.z), n);
                    let min = r + 2;
                    let max = r + 2 * (d as u64) + 1;
                    for i in min..max {
                        if primes[i as usize] {
                            let delta: usize = ((i as usize) - (r as usize)) / 2; // Distance to next prime
                            let f = Integer::from(
                                Integer::from(&s.x - &points[d].x)
                                    * Integer::from(&s.z + &points[d].z),
                            ) - &alpha
                                + &beta[delta];
                            g = multiply_mod(&g, &f, n);
                        }
                    }

                    let tmp = s.clone();
                    s = s.addh(&points[d], &t);
                    t = tmp;
                }
                g = g.gcd(&n);

                if 1 < g && g < *n {
                    info!("Sigma={}", sigma);
                    debug!("found factor {}\n\n\n", g);
                    found_factor.swap(true, Ordering::Relaxed);
                    return Some(g);
                }
            }
        }
    }
    None
}

pub fn ecm_singlethreaded(
    n: &Integer,
    max_curves: &Option<Integer>,
    b1: u64,
    b2: u64,
    sigma: &Option<Integer>,
) -> Option<Integer> {
    let d: usize = (b2 as f64).sqrt() as usize;

    let limit: usize = b2 as usize + 2 * d + 1; // not correct, assumes 64-bit architecture
    info!("Computing up to the {}th prime", limit);
    let mut primes = vec![true; limit];
    eratosthenes(&mut primes, limit);
    info!("Done");

    inversionless_ecm(
        &n,
        &max_curves,
        &primes,
        b1,
        b2,
        &sigma,
        0,
        &AtomicBool::new(false),
    )
}

pub fn ecm_multithreaded(
    n: &Integer,
    max_curves: &Option<Integer>,
    b1: u64,
    b2: u64,
    sigma: &Option<Integer>,
    nthreads: usize,
) -> Option<Integer> {
    let d: usize = (b2 as f64).sqrt() as usize;

    let limit: usize = b2 as usize + 2 * d + 1; // not correct, assumes 64-bit architecture
    info!("Computing up to the {}th prime", limit);
    let mut primes = vec![true; limit];
    eratosthenes(&mut primes, limit);
    info!("Done");

    let found_factor = Arc::new(AtomicBool::new(false));

    let mut children = vec![];
    let n = Arc::new(Integer::from(n));
    for i in 0..nthreads {
        let n = n.clone();
        let curves = max_curves.clone();
        let primes = primes.clone();
        let sigma = sigma.clone();
        let found_factor = Arc::clone(&found_factor);
        // Spin up another thread
        children.push(thread::spawn(move || -> Option<Integer> {
            inversionless_ecm(&n, &curves, &primes, b1, b2, &sigma, i, &found_factor)
        }))
    }
    let mut found = None;
    for child in children {
        // Wait for the thread to finish.
        // Once a thread found a factor, it notifies
        // all the other threads to terminate by means
        // of the found_factor atomic boolean.
        match child.join() {
            Ok(None) | Err(_) => (),
            Ok(f) => found = f,
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Find a factor of Fermat numbers F_n = 2^(2^n) + 1
    /// for n = 5 up to 8.
    #[test]
    fn check_ecm_fermat() {
        let b1 = 10000;
        let b2 = 100 * b1;
        for i in 5..8 {
            let fermat = Integer::from(Integer::u_pow_u(2, 2u32.pow(i))) + 1;
            match ecm_singlethreaded(&fermat, &None, b1, b2, &None) {
                Some(factor) => {
                    print!("got {}\n", factor);
                    assert_eq!(div_mod(&fermat, &factor).1, Integer::from(0))
                }
                None => (),
            }
        }
    }
}
