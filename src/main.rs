
use rand::Rng;
use std::{time::{Duration, Instant}, primitive};

const VALUES_RANGE: (i32, i32) = (-10, 10);
const BASE_CASE: usize = 64;

struct Matrix {
    els: Vec<i32>,
    size: (usize, usize)
}

impl Matrix {
    fn new(n: usize) -> Self {
        let mut els = Vec::new();
        let mut rng = rand::thread_rng();
        for _i in 0..(n*n) { 
            els.push(rng.gen_range(VALUES_RANGE.0..=VALUES_RANGE.1))
        }
        Self{els, size: (n, n)}
    }
    fn zero_filled(n: usize) -> Self {
        let mut els = Vec::new();
        for _i in 0..(n*n) { 
            els.push(0)
        }
        Self{els, size: (n, n)}
    }
    fn get_el(&self, x: usize, y: usize) -> i32 {
        self.els[y * self.size.0 + x ]
    }
    fn set_el(&mut self, x: usize, y: usize, value: i32) {
        self.els[y * self.size.0 + x ] = value;
    }
    fn multiply(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for k in 0..n {
                for j in 0..m.size.1 {
                    res.set_el(i, j, res.get_el(i, j) + self.get_el(i, k) * m.get_el(k, j))
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
                res.set_el(i, j, self.get_el(i, j) + m.get_el(i, j))
            }
        }
        res
    }

    fn sub(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n {
                res.set_el(i, j, self.get_el(i, j) - m.get_el(i, j))
            }
        }
        res
    }

    fn strassen_recursion(&self, m: &Matrix) -> Matrix {
        let n = self.size.0;
        if n <= BASE_CASE {
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
            let p1 = a_res.strassen_recursion(&b_res);

            a_res = a21.add(&a22);
            let p2 = a_res.strassen_recursion(&b11);

            b_res = b12.sub(&b22);
            let p3 = a11.strassen_recursion(&b_res);

            b_res = b21.sub(&b11);
            let p4 = a22.strassen_recursion(&b_res);

            a_res = a11.add(&a12);
            let p5 = a_res.strassen_recursion(&b22);

            a_res = a21.sub(&a11);
            b_res = b11.add(&b12);
            let p6 = a_res.strassen_recursion(&b_res);

            a_res = a12.sub(&a22);
            b_res = b21.add(&b22);
            let p7 = a_res.strassen_recursion(&b_res);

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

    fn multiply_strassen (&self, m: &Matrix) -> Matrix {
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
        let c_prep = a_prep.strassen_recursion(&b_prep);
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

impl std::fmt::Display for Matrix { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut output = String::new();
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                let value = self.get_el(i, j);
                output.push_str(&format!("{:5}", value));
            }
            output.push('\n');
        }
        write!(f, "{}", output)
    }
}

fn main() {

    let run_times1 = 1;
    let run_times2 = 1;
    let mut n = 1024;
    while n < 1500 {
        let m1 = Matrix::new(n);
        let m2 = Matrix::new(n);

        let start = Instant::now();
        for _i in 0..run_times1 {
            let _res = m1.multiply(&m2);
        }
        let duration = start.elapsed();
        let formatted = format!("loops;{n};{}", duration.as_secs_f64() / (run_times1 as f64)).replace(".", ",");
        println!("{}", formatted);


        let start = Instant::now();
        for _i in 0..run_times2 {
            let _res = m1.multiply_strassen(&m2);
        }
        let duration = start.elapsed();
        let formatted = format!("strassen;{n};{}", duration.as_secs_f64() / (run_times2 as f64)).replace(".", ",");
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