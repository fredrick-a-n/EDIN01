use std::io::*;

fn main() {
    exercise5();
}

fn exercise5() {
    let mut z2 = vec![0, 0, 0, 1];
    let mut z5 = vec![0, 0, 0, 1];
    // create file sequence.txt
    for i in 0..10003 {
        z2.push(
            match z2.clone()[i..] {
                // Special case of 1000->0000->0001
                [0, 0, 0, 0] => 1,
                [1, 0, 0, 0] => 0,
                // General case x^4 + x^3 + 1
                [a,b,_,_] => (-(a+b) as i32).rem_euclid(2),
                _ => unreachable!()
            }
        );
        z5.push(
            match z5.clone()[i..] {
                // Special case of 2000->0000->0001
                [0, 0, 0, 0] => 1,
                [2, 0, 0, 0] => 0,
                // General case 2x^4 + 2x^3 + x^2 + 1
                [a,b,c, _] => (-(2*a+2*b+c) as i32).rem_euclid(5),
                _ => unreachable!()
            }
        );
    }
    // print to file
    let mut file = std::fs::File::create("output.txt").unwrap();
    file.write_all(
        z2.iter().zip(z5.iter()).skip(4)
        .map(|(a,b)| {
            return (5*a + b).to_string();
        })
        .reduce(|a,b| 
            a + "\n" + &b
        ).unwrap().as_bytes()
    ).unwrap();
}