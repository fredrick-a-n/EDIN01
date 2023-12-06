

pub fn exercise1() {
    let num = include_str!("../task06.txt").chars().map(|c| c.to_digit(10).unwrap() as u8).collect::<Vec<u8>>();

    // Primitive polynomials
    let prim: Vec<Vec<u8>> = vec![
        vec![1,0,1,1,0,0,1,1,0,1,0,1,1], 
        vec![1,0,1,0,1,1,0,0,1,1,0,1,0,1,0],
        vec![1,1,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0]];

    // Generate de Bruijn sequences for the primitive polynomials
    let seq: Vec<Vec<u8>> = prim.iter()
        .map(|p| lfsr(p, vec![0; p.len()], 2_usize.pow(p.len() as u32))).collect();

    // Find the position of the states with maximum correlation to the given number
    let pos = seq.iter().map(|s| max_p(s, &num)).collect::<Vec<usize>>();

    // Get the specific states of the sequence with the maximum correlation
    let states = seq.iter().zip(pos.iter()).map(|(s, i)| s[*i..*i+num.len()].to_vec()).collect::<Vec<Vec<u8>>>();

    // Confirm that the three sequences generate the given number
    if check_seqs(num.clone(), &states){
        println!("Found!");
        println!("State 1: {:?}", seq[0][pos[0]..pos[0]+prim[0].len()].to_vec());
        println!("State 2: {:?}", seq[1][pos[1]..pos[1]+prim[1].len()].to_vec());
        println!("State 3: {:?}", seq[2][pos[2]..pos[2]+prim[2].len()].to_vec());
    } else {
        println!("Not found!");
    }
}

// Checks if the sequences generate the given number with majority vote
fn check_seqs(num: Vec<u8>, seq: &Vec<Vec<u8>>) -> bool{
    num == seq[0].iter()
        .zip(seq[1].iter())
        .zip(seq[2].iter())
        .map(|((x,y),z)| (*x+*y+*z) / 2)
        .collect::<Vec<u8>>()
}

// Finds the position of the maximum correlation between the given sequence and the given number
fn max_p(seq: &Vec<u8>, num: &Vec<u8>) -> usize {
    let mut dists: Vec<f32> = Vec::new();
    for i in 0..(seq.len()-num.len()) {
        let j = i + num.len();
        let state = seq[i..j].to_vec();
        let dist = distance(state, num.clone());
        dists.push(dist);
    }
    let (pos, max) = dists.iter().enumerate().max_by(|(_, x),(_, y)| x.partial_cmp(y).unwrap()).unwrap();
    println!("Found max: {}", max);
    pos
}


// Generates a lfsr sequence of length len, starting with the state init and using prim as the primitive polynomial
fn lfsr(prim: &Vec<u8>, init: Vec<u8>, len: usize) -> Vec<u8> {
    let mut seq: Vec<u8> = init.clone();
    for _ in 0..len {
        let last = seq.as_slice()[seq.len()-prim.len()..].to_vec();
        if last[1..].to_vec() == vec![0_u8; prim.len()-1] {
            // special case of 0 state
            seq.push(if last[0]==1 {0} else {1});
        } else {
            // general case
            seq.push(and(last,prim.clone()).iter().sum::<u8>() % 2);
        }
    }
    seq
}

// Bitwise and of two vectors
fn and(a: Vec<u8>, b: Vec<u8>) -> Vec<u8> {
    let mut out: Vec<u8> = vec![0; a.len()];
    for i in 0..a.len() {
        out[i] = a[i] * b[i];
    }
    out
}

// Calculates the distance between two vectors
fn distance(a: Vec<u8>, b: Vec<u8>) -> f32 {
    let n = a.len() as f32;
    1.0 - hamming(a,b) as f32 / n as f32
}

// Calculates the hamming distance between two vectors
fn hamming(a: Vec<u8>, b: Vec<u8>) -> u8 {
    let mut out: u8 = 0;
    for i in 0..a.len() {
        out += (a[i] ^ b[i]) as u8;
    }
    out
}