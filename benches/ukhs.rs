#[macro_use]
extern crate criterion;

use criterion::{Bencher, Criterion, Fun};
use rand::distributions::{Distribution, Uniform};
use ukhs::UKHS;

fn ukhs_bench(c: &mut Criterion) {
    let range = Uniform::from(0..4);
    let mut rng = rand::thread_rng();
    let seq = (0..10000)
        .map(|_| match range.sample(&mut rng) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        })
        .collect::<String>();

    let ukhs_it = Fun::new("ukhs_iterator", |b: &mut Bencher, i: &String| {
        let ukhs = UKHS::new(7, 20).unwrap();
        b.iter(|| {
            //  iter.for_each(drop);
            let iter = ukhs.iter_sequence(i.as_bytes());
            let _res: Vec<(String, String)> = iter.collect();
        })
    });

    let ukhs_hash_it = Fun::new("ukhs_hash_iterator", |b: &mut Bencher, i: &String| {
        let ukhs = UKHS::new(7, 20).unwrap();
        b.iter(|| {
            let iter = ukhs.hash_iter_sequence(i.as_bytes()).unwrap();
            //  iter.for_each(drop);
            let _res: Vec<(u64, u64)> = iter.collect();
        })
    });

    let functions = vec![ukhs_it, ukhs_hash_it];
    c.bench_functions("ukhs", functions, seq);
}

criterion_group!(benches, ukhs_bench);
criterion_main!(benches);
