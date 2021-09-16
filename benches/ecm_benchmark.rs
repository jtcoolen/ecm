use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ecm::arith::misc::*;
use ecm::inversionless_ecm;
use log::info;
use rug::Integer;
use std::sync::atomic::AtomicBool;

pub fn ecm_f6_benchmark(c: &mut Criterion) {
    let limit = 10000000;
    info!("Computing up to the {}th prime", limit);
    let mut primes = vec![true; limit];
    eratosthenes(&mut primes, limit);
    info!("Done.");
    let fermat = Integer::from(Integer::u_pow_u(2, 2u32.pow(6))) + 1;
    let b1 = 10000;
    let b2 = 100 * b1;
    c.bench_function("ecm F_6", |b| {
        b.iter(|| {
            inversionless_ecm(
                black_box(&fermat),
                black_box(&None),
                black_box(&primes),
                black_box(b1),
                black_box(b2),
                black_box(&None),
                black_box(0),
                black_box(&AtomicBool::new(false)),
            )
        })
    });
}

criterion_group!(benches, ecm_f6_benchmark);
criterion_main!(benches);
