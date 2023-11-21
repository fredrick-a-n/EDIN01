use std::{fs, vec};
use std::path::Path;
use std::f64::consts::{E, FRAC_1_SQRT_2};
use std::cmp;
use gcd::Gcd;
use std::env;
use num_bigint::BigUint;

pub fn exercise4() {
    env::set_var("RUST_BACKTRACE", "1");
    let _contents = fs::read_to_string(Path::new("./task06.txt")).expect("Should have been able to read the file");
    // let num: u128 = contents.chars()
    //     .filter(|a| a.to_digit(10).is_some())
    //     .collect::<String>().parse().unwrap();
    let num: u128 = 31741649;
    let factor = quadratic_sieve(num);
    let factor_2 = num / factor;

    if factor == 1 {
        print!("No factors found. Try increasing size of B.")
    } else {
        print!("The factors are: {} and {}", factor, factor_2);
    }
}  


fn quadratic_sieve(num: u128) -> u128{
    // let b = gen_b(num as f64); b can be far too large in this case, the program will run out of memory
    let b = cmp::min(1000, gen_b(num as f64)); // can be increased if needed
    //let b_smooth_vec:Vec<u128> = gen_base(num, b);
    let b_smooth_vec:Vec<u128> = gen_primes(b);
    println!("gen_b: {}", gen_b(num as f64));

    // Need better way to genereate b-smooth numbers

    let mut r_vec: Vec<(u128, u128)> = Vec::new();
    let mut k = 0;
    let comp = cmp::max(gen_b(num as f64) as usize, b as usize) + 10;
    while r_vec.len() < comp{
        let mut j = 0;
        while j < k && r_vec.len() < comp {
            let r = gen_r(num, j, k);
            let r2 = (r*r)%num;
            if r2 > 0 {
                if check_smoothness(r2, b_smooth_vec.clone()) {
                    r_vec.push((r,r2));
                }   
            }
            j += 1;
        }
        k += 1;
    }


    //let mut r_vec: Vec<(u128, u128)> = find_smooth(&b_smooth_vec, num, gen_b(num as f64));
    
    let (squares, matrix) = create_matrix(b_smooth_vec, r_vec.clone());


    for q in squares.clone() {
        let x = r_vec[q].1;
        let factor = x.gcd(num);
        if factor != 1 {
            println!("Found square factor.");
            return factor;
        }
    }
    println!("Found no square factors. Continuing...");


    let (sol_rows, marks, g_m) = gauss(&matrix);

    println!("Gauss done. Finding solution vector...");

    let mut solution_vec = solve_row(sol_rows.clone(), &g_m, 0, marks.clone());

    let mut tries = 0;
    let mut factor = solve(solution_vec.clone(), r_vec.clone(), num);

    for (i, _val) in sol_rows.clone().iter().enumerate() {
        if factor == 1 || factor == num {
            tries += 1;
            println!("Didn't work. Trying different solution vector...");
            solution_vec = solve_row(sol_rows.clone(), &g_m, i as u128, marks.clone());
            factor = solve(solution_vec, r_vec.clone(), num);
        } else {
            println!("Found factors after {} tries!", tries);
            return factor;
        }
    }

    println!("found no factors after {} tries.", tries);
    return 1;
}

fn check_smoothness(num: u128, b_smooth_vec: Vec<u128>) -> bool {
    let mut n = num;
    let mut i = 0;
    while i < b_smooth_vec.len() {
        if n >= b_smooth_vec[i] && n % b_smooth_vec[i] == 0 {
            n = n / b_smooth_vec[i];
        } else {
            i += 1;
        }
        if n == 1 {
            break;
        }
        if n == 0 {
            break;
        }
    }
    return n == 1 ;
}

// from https://youtu.be/x5LTBsmGfFc
fn gen_b(n: f64) -> u128 {
    let k = n.log2();
    return E.powf(FRAC_1_SQRT_2*(k*k.log2()).sqrt()) as u128;
}


fn gen_base(n: u128, b: u128) -> Vec<u128> {
    let mut base: Vec<u128> = Vec::new();
    let primes = gen_primes(b);
    for p in primes.clone() {
        if quad_residue(n as i128, p as i128) == 1 {
            base.push(p);
        }
    }
    // for p in primes: # such that N is a quadratic residue mod p
    //   if quad_residue(N,p) == 1:
    //     factor_base.append(p)
    return base;
}

// sieving for primes
fn gen_primes(num: u128) -> Vec<u128> {
    let mut primes: Vec<u128> = Vec::new();
    let mut sieve = vec![true; num as usize];
    let mut i = 2;
    while i * i <= num {
        if sieve[i as usize] {
            let mut j = i * i;
            while j < num {
                sieve[j as usize] = false;
                j += i;
            }
        }
        i += 1;
    }
    for i in 2..num {
        if sieve[i as usize] {
            primes.push(i);
        }
    }
    return primes;
}

// from project description
fn gen_r(n: u128, j: u128, k: u128) -> u128 {
    let root = ((n*k) as f64).sqrt().floor() as u128;
    return root + j;
}


fn create_matrix(b_smooth_vec: Vec<u128>, r_vec: Vec<(u128, u128)>) -> (Vec<usize>, Vec<Vec<u8>>) {
    let mut matrix: Vec<Vec<u8>> = vec![vec![0; b_smooth_vec.len()]; r_vec.len()];
    //let mut matrix: Vec<Vec<u128>> = Vec::new();
    let mut squares: Vec<usize> = Vec::new();
    for i in 0..((r_vec.len()-1) as usize) {
        let mut square_count = 0;
        for j in 0..((b_smooth_vec.len() - 1) as usize) {
            let mut e = 0;
            let mut n = r_vec[i].1;
            while n % b_smooth_vec[j] == 0 {
                e += 1;
                n = n / b_smooth_vec[j];
            }
            matrix[i][j] = (e % 2) as u8;
            square_count += e%2;
        }
        if square_count == 0 {
            squares.push(i)
        }
    }

    return (squares, matrix);
}


fn gauss(m: &Vec<Vec<u8>>) -> (Vec<(Vec<u8>, usize)>, Vec<bool>, Vec<Vec<u8>>) {
    let mut marks: Vec<bool> = vec![false; m[0].len()];
    let mut g_m = m.clone();
    
    for (i, _) in g_m.clone().iter().enumerate() {
        for (j, _) in g_m[i].clone().iter().enumerate() {
            if g_m[i][j] == 1 {
                marks[j] = true;
                for (k, _) in g_m.clone().iter().enumerate() {
                    if k != i && g_m[k][j] == 1 {
                        for (l, _) in g_m[k].clone().iter().enumerate() {
                            g_m[k][l] = (g_m[k][l] + g_m[i][l]) % 2;
                        }
                    }
                }
                break;
            }
        }
    }


    // println!("marks: {:?}", marks);
    // println!("g_m: {:?}", g_m);
    // println!("m: {:?}", m);

    g_m = transpose(g_m);

    // for (i, val) in g_m.iter().enumerate(){
    //     println!("{}: {:?}", marks[i], val);
    // }

    let mut sol_rows: Vec<(Vec<u8>, usize)> = Vec::new();
    for (i, mark) in marks.iter().enumerate() {
        if !mark {

            sol_rows.push((g_m[i].clone(), i));
        }
    }
    println!("sol_rows: {:?}", sol_rows);

    return (sol_rows, marks, g_m);
}


fn transpose(matrix: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut new_matrix: Vec<Vec<u8>> = Vec::new();
    for i in 0..matrix[0].len() {
        let mut new_row: Vec<u8> = Vec::new();
        for row in matrix.clone() {
            new_row.push(row[i]);
        }
        new_matrix.push(new_row);
    }
    return new_matrix;
}

fn solve_row(sol_rows: Vec<(Vec<u8>, usize)>, m: &Vec<Vec<u8>>, k: u128, marks: Vec<bool>) -> Vec<usize> {
    let mut solution_vec: Vec<usize> = Vec::new();
    let mut indices: Vec<usize> = Vec::new();
    let free_row = sol_rows[k as usize].0.clone();
    for (i, num) in free_row.iter().enumerate() {
        if *num == 1 {
            indices.push(i);
        }
    }
    for (r, row) in m.iter().enumerate() {
        for i in indices.clone() {
            if row[i] == 1 && marks[r] {
                solution_vec.push(r);
                break;
            }
        }
    }
    solution_vec.push(sol_rows[k as usize].1);
    return solution_vec;
}


// def solve_row(sol_rows,M,marks,K=0):
//     solution_vec, indices = [],[]
//     free_row = sol_rows[K][0] # may be multiple K
//     for i in range(len(free_row)):
//         if free_row[i] == 1: 
//             indices.append(i)
//     for r in range(len(M)): #rows with 1 in the same column will be dependent
//         for i in indices:
//             if M[r][i] == 1 and marks[r]:
//                 solution_vec.append(r)
//                 break
            
//     solution_vec.append(sol_rows[K][1])       
//     return(solution_vec)


fn solve(solution_vec: Vec<usize>, r_vec: Vec<(u128,u128)>, n: u128) -> u128 {
    println!("solution_vec: {:?}", solution_vec);
    let mut solution_nums: Vec<u128> = Vec::new();
    let mut x_nums: Vec<u128> = Vec::new();
    for i in solution_vec.clone() {
        solution_nums.push(r_vec[i].1);
        x_nums.push(r_vec[i].0);
    }

    let mut a_square = BigUint::from(1u128);
    for n in solution_nums.clone() {
        a_square *= n;
    }

    let mut b = BigUint::from(1u128);
    for n in x_nums.clone() {
        b *= n;
    }

    let a = a_square.sqrt();

    // println!("a: {}", a);
    // println!("b: {}", b);


    let factor = bgcd(b-a,BigUint::from(n));

    let factor = biguint_to_u64(factor);

    // println!("factor: {}", factor);
    return factor as u128;
}

fn bgcd(a: BigUint, b: BigUint) -> BigUint {
    let mut a = a.clone();
    let mut b = b.clone();
    while b != BigUint::from(0 as u8) {
        let t = b.clone();
        b = a % b;
        a = t;
    }
    return a;
}


// def solve(solution_vec,smooth_nums,xlist,N):
    
//     solution_nums = [smooth_nums[i] for i in solution_vec]
//     x_nums = [xlist[i] for i in solution_vec]
    
//     Asquare = 1
//     for n in solution_nums:
//         Asquare *= n
        
//     b = 1
//     for n in x_nums:
//         b *= n

//     a = isqrt(Asquare)
    
//     factor = gcd(b-a,N)
//     return factor






// fn find_smooth(b_smooth_vec: &Vec<u128>, num: u128, interval: u128) -> Vec<(u128,u128)> {
//     let mut sieve_seq: Vec<i128> = Vec::new();
//     let root = (num as f64).sqrt().floor() as u128;
//     for i in root..(root+interval) {
//         sieve_seq.push(((i*i) as i128) - (num as i128));
//     }
//     let mut sieve_list = sieve_seq.clone();

//     if b_smooth_vec[0] == 2 {
//         let mut i = 0;
//         while sieve_list[i] % 2 != 0 {
//             i += 1;
//         }
//         for j in (i..sieve_list.len()).step_by(2) {
//             while sieve_list[j] % 2 == 0 {
//                 sieve_list[j] /= 2;
//             }
//         }
        
//     }

//     for p in b_smooth_vec[1..].iter() {
//         println!("p: {}  n: {}", p, num as i128);
//         let residues =  ModSqrt::mod_sqrt(num as i128, *p as i128).unwrap();
        
//         //stonelli(BigInt::from(num), BigInt::from(*p));
//         for i in (i128::abs((residues.0-root as i128) % (*p as i128))..(sieve_list.len() as i128)).step_by(*p as usize) {
//             while sieve_list[i as usize] % (*p as i128) == 0 {
//                 sieve_list[i as usize] /= *p as i128;
//             }
//         }
//         for i in (i128::abs((residues.1-root as i128) % (*p as i128))..(sieve_list.len() as i128)).step_by(*p as usize) {
//             while sieve_list[i as usize] % (*p as i128) == 0 {
//                 sieve_list[i as usize] /= *p as i128;
//             }
//         }
        
//     }

//     let mut r_vec: Vec<(u128, u128)> = Vec::new();
//     for (i, num) in sieve_list.iter().enumerate() {
//         if r_vec.len() >= b_smooth_vec.len() + 10 {
//             break;
//         }
//         if *num == 1 {
//             r_vec.push((root+i as u128, sieve_seq[i] as u128));
//         }
//     }

//     return r_vec;
    
// }



// def find_smooth(factor_base,N,I):
// # tries to find B-smooth numbers in sieve_seq, using sieving

//     def sieve_prep(N,sieve_int):
//     # generates a sequence from Y(x) = x^2 - N, starting at x = root 
//         sieve_seq = [x**2 - N for x in range(root,root+sieve_int)]
//         #sieve_seq_neg = [x**2 - N for x in range(root,root-sieve_int,-1)]
//         return sieve_seq

//     sieve_seq = sieve_prep(N,I)
//     sieve_list = sieve_seq.copy() # keep a copy of sieve_seq for later
//     if factor_base[0] == 2:
//         i = 0
//         while sieve_list[i] % 2 != 0:
//             i += 1
//         for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
//             while sieve_list[j] % 2 == 0: #account for powers of 2
//                 sieve_list[j] //= 2
//     #print("")
//     for p in factor_base[1:]: #not including 2
//         residues = STonelli(N,p) #finds x such that x^2 = n (mod p). There are two start solutions
        
//         for r in residues:
//             for i in range((r-root) % p, len(sieve_list), p): # Now every pth term will also be divisible
//                 while sieve_list[i] % p == 0: #account for prime powers
//                     sieve_list[i] //= p
//     xlist = [] #original x terms
//     smooth_nums = []
//     indices = [] # index of discovery
    
//     for i in range(len(sieve_list)):
//         if len(smooth_nums) >= len(factor_base)+T: #probability of no solutions is 2^-T
//             break
//         if sieve_list[i] == 1: # found B-smooth number
//             smooth_nums.append(sieve_seq[i])
//             xlist.append(i+root)
//             indices.append(i)

//     return(smooth_nums,xlist,indices)


fn quad_residue(n: i128, p: i128) -> i128 {
    let l = 1;
    let q = (p-1)/2;
    let mut x = q.pow(l as u32);
    if x == 0 {
        return 1;
    }
    let mut a = n % p;
    let mut z = 1;
    while x != 0 {
        if x % 2 == 0 {
            a = (a*a) % p;
            x /= 2;
        } else {
            x -= 1;
            z = (z*a) % p;
        }
    }
    return z;
}



fn biguint_to_u64(num: BigUint) -> u64 {
    if num.bits() > 64 {
        panic!("factor too large");
    }
    // let mut arr: [u8; 16] = [0; 16];
    // for (i, byte) in factor.to_bytes_be().iter().enumerate() {
    //     arr[15-i] = *byte;
    //     println!("byte: {}", byte);
    // }
    // let factor = u128::from_be_bytes(arr);
   return num.to_u64_digits()[0].clone();
}


// // Arguments n, p as described in Wikipedia (WP)
// fn ts(n: i64, p: i64) -> (i64, i64, bool) {
    
//     if ls(n, p) != 1 { return (0, 0, false) }
    
//     // WP step 1, factor out powers two.
//     // variables Q, S named as at WP.
//     let mut q = p - 1;
//     let mut s = 0;
//     while q & 1 == 0 {
//         s += 1;
//         q >>= 1
//     }
    
//     // WP step 1, direct solution
//     if s == 1 {
//         let r1 = power_mod(n, (p+1)/4, p);
//         return (r1, p - r1, true)
//     }
    
//     // WP step 2, select z, assign c
//     let mut z = 2;
//     while ls(z, p) != p-1 { z += 1 }
//     let mut c = power_mod(z, q, p);
    
//     // WP step 3, assign R, t, M
//     let mut r = power_mod(n, (q+1)/2, p);
//     let mut t = power_mod(n, q, p);
//     let mut m = s;
    
//     // WP step 4, loop
//     loop {
//     // WP step 4.1, termination condition
//         if t == 1 { return (r, p - r, true) }
        
//         // WP step 4.2, find lowest i...
//         let mut i = 0;
//         let mut z = t;
//         while z != 1 && i < e {
//             b = b * b % p;
//             e -= 1
//         }
//         r = r * b % p;
//         c = b * b % p; // more convenient to compute c before t
//         t = t * c % p;
//         m = i;
//     }
// }

// // Legendre symbol, returns 1, 0, or -1 mod p
// fn ls(a: i64, p: i64) -> i64 {
//     power_mod(a, (p-1)/2, p)
// }

// fn power_mod(base: i64, exp: i64, modulus: i64) -> i64 {
//     if exp < 0 { unimplemented!() } let mut b = base % modulus; let mut result = 1; let mut e = exp; while e > 0 {
//     if e & 1 != 0 { result = (result * b) % modulus }
//         b = (b * b) % modulus;
//         e >>= 1
//     }
//     modulo(result, modulus)
// }

// fn modulo(n: i64, m: i64) -> i64 {
//     if n >= 0 { n % m } else { - (-n) % m }
// }