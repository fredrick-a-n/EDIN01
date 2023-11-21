mod exercise4;

fn main() {
    use std::time::Instant;
    let now = Instant::now();
    exercise4::exercise4();
    println!("Took {}ms", now.elapsed().as_millis());
}
