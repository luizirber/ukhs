#[macro_use]
extern crate criterion;

use criterion::{Bencher, Criterion, Fun};
use rand::distributions::{Distribution, Uniform};
use ukhs::{UKHSHashIterator, UKHSIterator};

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
        b.iter(|| {
            let iter = UKHSIterator::new(i.as_bytes(), 7, 20).unwrap();
            //  iter.for_each(drop);
            let _res: Vec<String> = iter.collect();
        })
    });

    let ukhs_hash_it = Fun::new("ukhs_hash_iterator", |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = UKHSHashIterator::new(i.as_bytes(), 7, 20).unwrap();
            //  iter.for_each(drop);
            let _res: Vec<u64> = iter.collect();
        })
    });

    let functions = vec![ukhs_it, ukhs_hash_it];
    c.bench_functions("ukhs", functions, seq);
}

criterion_group!(benches, ukhs_bench);
criterion_main!(benches);
