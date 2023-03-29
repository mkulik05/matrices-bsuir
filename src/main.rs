
use num_bigint::BigInt;
use num_bigint::{ToBigInt, RandBigInt};
use num_traits::{Zero, One};
use rand::Rng;
use std::collections::btree_map::Values;
use std::{time::{Duration, Instant}, primitive};

struct Matrix {
    els: Vec<Box<BigInt>>,
    size: (usize, usize)
}

impl Matrix {
    fn new(n: usize, VALUES_RANGE: &(BigInt, BigInt)) -> Self {
        let mut els = Vec::new();
        let mut rng = rand::thread_rng();
        for _i in 0..(n*n) { 
            let v = rng.gen_bigint_range(&VALUES_RANGE.0, &VALUES_RANGE.1);
            // println!("{}", v);
            els.push(Box::new(v))
        }
        Self{els, size: (n, n)}
    }
    fn zero_filled(n: usize) -> Self {
        let mut els = Vec::new();
        for _i in 0..(n*n) { 
            els.push(Box::new(Zero::zero()))
        }
        Self{els, size: (n, n)}
    }
    fn get_el(&self, x: usize, y: usize) -> Box<BigInt> {
        self.els[y * self.size.0 + x ].clone()
    }
    fn set_el(&mut self, x: usize, y: usize, value: Box<BigInt>) {
        self.els[y * self.size.0 + x ] = value;
    }
    fn multiply(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for k in 0..n {
                for j in 0..m.size.1 {
                    res.set_el(i, j, Box::new(*res.get_el(i, j) + (*self.get_el(i, k)) * &(*m.get_el(k, j))))
                }
            }
        }
        res
    }
    fn add(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n {
                res.set_el(i, j, Box::new(*m.get_el(i, j) + &*self.get_el(i, j)))
            }
        }
        res
    }

    fn sub(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n {
                res.set_el(i, j, Box::new(*self.get_el(i, j) - &*m.get_el(i, j)))
            }
        }
        res
    }

    fn strassen_recursion(&self, m: &Matrix, base_case: i32) -> Matrix {
        let n = self.size.0;
        if n <= (base_case as usize) {
            return self.multiply(m)
        } else {
            let new_n = n / 2;
            let mut a11 = Matrix::zero_filled(new_n);
            let mut a12 = Matrix::zero_filled(new_n);
            let mut a21 = Matrix::zero_filled(new_n);
            let mut a22 = Matrix::zero_filled(new_n);

            let mut b11 = Matrix::zero_filled(new_n);
            let mut b12 = Matrix::zero_filled(new_n);
            let mut b21 = Matrix::zero_filled(new_n);
            let mut b22 = Matrix::zero_filled(new_n);

            let mut a_res = Matrix::zero_filled(new_n);
            let mut b_res = Matrix::zero_filled(new_n);

            for i in 0..new_n {
                for j in 0..new_n { 
                    a11.set_el(i, j, self.get_el(i, j));
                    a12.set_el(i, j, self.get_el(i, j + new_n));
                    a21.set_el(i, j, self.get_el(i + new_n, j));
                    a22.set_el(i, j, self.get_el(i + new_n, j + new_n));

                    b11.set_el(i, j, m.get_el(i, j));
                    b12.set_el(i, j, m.get_el(i, j + new_n));
                    b21.set_el(i, j, m.get_el(i + new_n, j));
                    b22.set_el(i, j, m.get_el(i + new_n, j + new_n));
                }
            }

            a_res = a11.add(&a22);
            b_res = b11.add(&b22);
            let p1 = a_res.strassen_recursion(&b_res, base_case);

            a_res = a21.add(&a22);
            let p2 = a_res.strassen_recursion(&b11, base_case);

            b_res = b12.sub(&b22);
            let p3 = a11.strassen_recursion(&b_res, base_case);

            b_res = b21.sub(&b11);
            let p4 = a22.strassen_recursion(&b_res, base_case);

            a_res = a11.add(&a12);
            let p5 = a_res.strassen_recursion(&b22, base_case);

            a_res = a21.sub(&a11);
            b_res = b11.add(&b12);
            let p6 = a_res.strassen_recursion(&b_res, base_case);

            a_res = a12.sub(&a22);
            b_res = b21.add(&b22);
            let p7 = a_res.strassen_recursion(&b_res, base_case);

            let c12 = p3.add(&p5);
            let c21 = p2.add(&p4);

            a_res = p1.add(&p4);
            b_res = a_res.add(&p7);
            let c11 = b_res.sub(&p5);

            a_res = p1.add(&p3);
            b_res = a_res.add(&p6);
            let c22 = b_res.sub(&p2);

            let mut c = Matrix::zero_filled(n);
            for i in 0..new_n {
                for j in 0..new_n {
                    c.set_el(i, j, c11.get_el(i, j));
                    c.set_el(i, j + new_n, c12.get_el(i, j));
                    c.set_el(i + new_n, j, c21.get_el(i, j));
                    c.set_el(i + new_n, j + new_n, c22.get_el(i, j));
                }
            }

            return c
        }
    }

    fn multiply_strassen (&self, m: &Matrix, base_case: i32) -> Matrix {
        let n = self.size.0;
        let power_base: i32 = 2;
        let new_n = power_base.pow((n as f32).log2().ceil() as u32);
        let mut a_prep = Matrix::zero_filled(new_n as usize);
        let mut b_prep = Matrix::zero_filled(new_n as usize);
        for i in 0..n {
            for j in 0..n { 
                a_prep.set_el(i, j, self.get_el(i, j));
                b_prep.set_el(i, j, m.get_el(i, j));
            }
        }
        let c_prep = a_prep.strassen_recursion(&b_prep, base_case);
        let mut c = Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n { 
                c.set_el(i, j, c_prep.get_el(i, j));
            }
        }
        
        c

    }

    fn is_equal(&self, other: &Matrix) -> bool {
        if self.size != other.size {
            return false;
        }
        let n = self.size;
        for i in 0..(n.0 * n.1) {
            if self.els[i] != other.els[i] {
                return false
            }
        }

        true
    }
}

// impl std::fmt::Display for Matrix { 
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         let mut output = String::new();
//         for i in 0..self.size.0 {
//             for j in 0..self.size.1 {
//                 let value = self.get_el(i, j);
//                 output.push_str(&format!("{:5}", value));
//             }
//             output.push('\n');
//         }
//         write!(f, "{}", output)
//     }
// }

fn bench_mul(m1: &Matrix, m2: &Matrix, min_time: f64, min_run_times: i32) -> f64 {
    let start = Instant::now();
    for _i in 0..min_run_times {
        let _res = m1.multiply(&m2);
    }
    let mut duration = start.elapsed();
    let mut n = 0;
    if duration.as_secs_f64() < min_time {
        let start = Instant::now();
        for i in 0..((min_time / duration.as_secs_f64() + 1.0) as i32) {
            let _res = m1.multiply(&m2);
            n += 1;
        }
        duration = start.elapsed();
    } else {
        n = min_run_times;
    }
    duration.as_secs_f64() / (n as f64)
}

fn bench_strassen(m1: &Matrix, m2: &Matrix, min_time: f64, min_run_times: i32, base_case: i32) -> f64 {
    let start = Instant::now();
    for _i in 0..min_run_times {
        let _res = m1.multiply_strassen(&m2, base_case);
    }
    let mut duration = start.elapsed();
    let mut n = 0;
    if duration.as_secs_f64() < min_time {
        let start = Instant::now();
        for i in 0..((min_time / duration.as_secs_f64() + 1.0) as i32) {
            let _res = m1.multiply_strassen(&m2, base_case);
            n += 1;
        }
        duration = start.elapsed();
    } else {
        n = min_run_times;
    }
    duration.as_secs_f64() / (n as f64)
}

fn main() {
    let values_range: (BigInt, BigInt) = ((-10000000000000000000000000000 as i128).to_bigint().unwrap(), (100000000000000000000000000000 as i128).to_bigint().unwrap());
    let base_case = 1;
    let min_run_times1 = 2;
    let min_run_times2 = 2;
    let min_time = 0.5;
    let mut n = 1;
    while n < 3500 {
        let m1 = Matrix::new(n, &values_range);
        let m2 = Matrix::new(n, &values_range);

        let t1 = bench_mul(&m1, &m2, min_time, min_run_times1);
        let formatted = format!("loops;{n};{}", t1).replace(".", ",");
        println!("{}", formatted);


        let t2 = bench_strassen(&m1, &m2, min_time, min_run_times2, base_case);
        let formatted = format!("stassen;{n};{}", t2).replace(".", ",");
        println!("{}", formatted);

        if n > 400 {
            n += 10;
        } else if n > 200 {
            n += 8;
        } else if n > 100 {
            n += 5;
        } else {
            n += 1;
        }
    }
}