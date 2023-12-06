mod exercise1;

fn main() {
    use std::time::Instant;
    let now = Instant::now();
    exercise1::exercise1();
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}

