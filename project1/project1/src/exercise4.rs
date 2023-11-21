use num::{Integer, One, Zero, ToPrimitive};
use num_bigint::BigUint;
use std::{cmp, sync::atomic};

pub fn exercise4() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_millis()
        .init();

    let num: BigUint = BigUint::parse_bytes(b"323", 10).unwrap();

    let factors = factor(num.clone()).unwrap();
    log::info!("Two factors of {num} are {factors:?}", factors = factors);
}  

fn factor(n: BigUint) -> Option<(BigUint, BigUint)> {
    // Increase b if needed
    let b = cmp::min(10000, cmp::max(2*n.to_f64().unwrap().sqrt().floor() as usize, 50));
    log::info!("Let's find some factors of N={n}");
    log::info!("Got bound B = {b}");
    // Find all prime numbers up to B
    let p: Vec<BigUint> = find_primes(b)
        .into_iter()
    // Filter all primes that aren't quadratic residue of n, makes the program a lot faster in general but makes it break for some small numbers.
        .filter(|p| legendre(&n, &BigUint::from(*p)) == 1)
        .map(BigUint::from)
        .collect();
    log::debug!("Filtered prime list is {p:?}");
    log::info!("Found {p_len} primes", p_len = p.len());

    // Find enough B-smooth numbers
    let x = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));

    let thread_count = 8;
    let prog: atomic::AtomicUsize = Default::default();
    let max = p.len() + 1;
    log::info!("Waiting for threads to finish, please wait");
    std::thread::scope(|t| {
        for i in 0..thread_count {
            let mut j = 1 + i as u32;
            let prog = &prog;
            let n = n.clone();
            let s = n.sqrt();
            let p = &p;
            let x = x.clone();
            t.spawn(move || {
                while prog.load(atomic::Ordering::SeqCst) < max  {
                    //let r: BigUint = s.clone() * BigUint::from((j as f64).sqrt().floor() as u32) + BigUint::from(j);
                    //let r2 = r.pow(2) % &n;
                    let r = &s + j;
                    let r2: BigUint = r.pow(2) - n.clone();
                    if let Some(r2_vec) = check_smoothness(r2.clone(), p) {
                        x.lock().unwrap().push((r, r2_vec));
                        prog.fetch_add(1, atomic::Ordering::SeqCst);
                    }

                    j += thread_count;
                }
            });
        }
    });
    let x = x.lock().unwrap().clone();
    log::debug!("Found B-smooth numbers: {:?}", x);

    let flat: Vec<i8> = x
        .iter()
        .flat_map(|x| x.1.clone())
        .map(|x| (x % 2) as i8)
        .collect();
    let matrix = ndarray::Array2::from_shape_vec((x.len(), p.len()), flat).unwrap();
    log::debug!("x: ");
    for y in x.iter() {
        log::debug!("{:?}", y);
    }

    // Each Coeff is a potential factor
    for coeffs in gauss_reduction(matrix) {
        log::debug!("Trying coeffs: {:?}", coeffs);
        assert!(coeffs.len() == x.len());
        let a = x
            .iter()
            .zip(coeffs.iter())
            .filter(|(_, coeff)| **coeff == 1)
            .map(|(a, _)| a.0.clone())
            .reduce(|x, y| x * y)
            .unwrap();
        let b_coeff_sq = x
            .iter()
            .zip(coeffs.iter())
            .map(|(a, coeff)| a.1.iter().map(|x| *coeff as u32 * (*x) as u32).collect())
            .reduce(|x: Vec<u32>, y: Vec<u32>| x.into_iter().zip(y).map(|(a, b)| a + b).collect())
            .unwrap();
        log::debug!("b_coeff_sq: {:?}", b_coeff_sq);
        let b_coeff = b_coeff_sq.into_iter().map(|x| x / 2); // THIS is why it works, we are SURE all coeffs are at this point even,
                                                             // thus, we can divide by 2, thus taking the sqrt
        let b = b_coeff
            .zip(&p)
            .map(|(c, p)| p.pow(c))
            .reduce(|x, y| x * y)
            .unwrap();
        log::debug!("a = {a}, b = {b}", a = a, b = b);
        log::debug!("a^2 - b^2 mod n = {}", (a.pow(2) - b.pow(2)) % &n);
        let (a, b) = if a > b { (a, b) } else { (b, a) };
        let diff = &a - &b;
        let sum = &a + &b;

        let factor1 = diff.gcd(&n);
        let factor2 = sum.gcd(&n);
        log::debug!("a={:?} mod n, b={:?} mod n", a % &n, b % &n);
        let f1 = factor1.clone().max(factor2.clone());
        let f2 = &n / f1.clone();

        log::info!("Might have found something: {factor1}, {factor2} => n = {f1}*{f2} = {n}");
        if f2 < n && f1 < n {
            return Some(if f1 < f2 { (f2, f1) } else { (f1, f2) });
        }
    }
    return None
}

// Checks if B-smooth, returns vector of exponants if so
fn check_smoothness(mut n: BigUint, p: &Vec<BigUint>) -> Option<Vec<u8>> {
    let mut r = vec![0_u8; p.len()];

    for (i, p) in p.iter().enumerate() {
        let mut j = 0;
        loop {
            let (d, r) = n.div_rem(p);
            if r == BigUint::zero() {
                n = d;
            } else {
                break;
            }
            j += 1;
        }
        r[i] = j;
    }
    if n == BigUint::one() {
        return Some(r)
    } else {
        return None
    }
}

// https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
fn find_primes(b: usize) -> Vec<u32> {
    log::info!("Sieving up to {b}");
    let mut p = Vec::new();
    let mut sieve = vec![true; b];
    for i in 2usize..b {
        if !sieve[i] {
            continue;
        }

        log::trace!("Found prime {i}");
        p.push(i as u32);
        for j in (i * i..b).step_by(i) {
            sieve[j] = false;
        }
    }
    log::debug!("Found primes {p:?}");
    return p
}

// Check if quadratic residue
fn legendre(n: &BigUint, p: &BigUint) -> u32 {
    if *p == 2u8.into() {
        return 1;
    }
    let (f, r) = (p - 1u8).div_rem(&BigUint::from(2u8));
    assert!(r == Zero::zero());
    return n.modpow(&f, p).to_u32().unwrap()
}

fn gauss_reduction(mut m: ndarray::Array2<i8>) -> Vec<Vec<i8>> {
    let mut vars = ndarray::Array2::eye(m.nrows());

    fn swap(m: &mut ndarray::ArrayViewMut2<i8>, index1: usize, index2: usize) {
        if index1 == index2 {
            return;
        }
        let (index1, index2) = if index1 > index2 {
            (index2, index1)
        } else {
            (index1, index2)
        };

        let mut it = m.axis_iter_mut(ndarray::Axis(0));
        ndarray::Zip::from(it.nth(index1).unwrap())
            .and(it.nth(index2 - (index1 + 1)).unwrap())
            .for_each(std::mem::swap);
    }

    fn gauss_reduction_rec(
        mut m: ndarray::ArrayViewMut2<i8>,
        mut vars: ndarray::ArrayViewMut2<i8>,
    ) {
        if m.is_empty() {
            return;
        }

        log::debug!("Gaussian reduction on matrix of size {:?}", m.shape());
        let mut first_one_pos = None;

        // Find first one if any
        for i in 0..m.nrows() {
            if m[[i, 0]] == 1 {
                first_one_pos = Some(i);
                break;
            }
        }
        if let Some(first_one_pos) = first_one_pos {
            swap(&mut m, 0, first_one_pos);
            swap(&mut vars, 0, first_one_pos);
            let mut it = m.axis_iter_mut(ndarray::Axis(0));
            let mut it_var = vars.axis_iter_mut(ndarray::Axis(0));
            let r0 = it.next().unwrap();
            let r0_var = it_var.next().unwrap();
            for (mut r, mut r_var) in it.zip(it_var) {
                if r[0] == 1 {
                    r -= &r0;
                    r.mapv_inplace(|x| i8::abs(x % 2));
                    r_var -= &r0_var;
                    r_var.mapv_inplace(|x| i8::abs(x % 2));
                }
            }

            let (_, subm) = m.split_at(ndarray::Axis(0), 1);
            let (_, subm_var) = vars.split_at(ndarray::Axis(0), 1);
            let (_, subm) = subm.split_at(ndarray::Axis(1), 1);
            //let (_, subm_var) = subm_var.split_at(ndarray::Axis(1), 1);
            gauss_reduction_rec(subm, subm_var)
        } else {
            log::trace!("Pivot not found");
            let (_, subm) = m.split_at(ndarray::Axis(1), 1);
            //let (_, subm_var) = vars.split_at(ndarray::Axis(1), 1);
            gauss_reduction_rec(subm, vars)
        }
    }
    gauss_reduction_rec(m.view_mut(), vars.view_mut());
    log::info!("Before gauss reduction \n {:?}", m);
    log::info!("Gaussian reduction done \n {:?}", vars);

    let mut v = vec![];
    for (r, r_var) in m
        .axis_iter(ndarray::Axis(0))
        .zip(vars.axis_iter(ndarray::Axis(0)))
    {
        if r == ndarray::Array1::zeros([r.len()]) {
            log::trace!("Got coeff candidate: {:?}", r_var.to_vec());
            let v_var: Vec<i8> = r_var.to_vec();
            v.push(v_var);
        }
    }

    if v.is_empty() {
        log::error!("No coeff found ?! {m}");
    }
    log::debug!("Gaussian reduction done 2 \n {:?}", v);
    return v
}
