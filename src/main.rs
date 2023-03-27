
use rand::Rng;
use std::{time::{Duration, Instant}, primitive};

const VALUES_RANGE: (i32, i32) = (1, 10);
const BASE_CASE: usize = 2;

struct Matrix {
    els: Vec<i32>,
    size: (usize, usize)
}

struct OperationCount {
    muls: u32,
    sums: u32,
    sums0: u32,
    muls0: u32
}

impl OperationCount {
    fn add(&self, other: &OperationCount) -> Self {
        OperationCount{muls: self.muls + other.muls, sums: self.sums + other.sums, sums0: self.sums0 + other.sums0, muls0: self.muls0 + other.muls0}
    }
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
    fn multiply(&self, m: &Matrix) -> (Matrix, OperationCount) {
        let mut sums = 0;
        let mut muls = 0;
        let mut sums0 = 0;
        let mut muls0 = 0;
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for k in 0..n {
                for j in 0..m.size.1 {
                    // if res.get_el(i, j) == 0 || self.get_el(i, k) * m.get_el(k, j) == 0 {
                    //     // sums += 1;
                    // } else {
                    //     sums += 1;
                    // }
                    sums += 1;
                    if self.get_el(i, k) == 0 || m.get_el(k, j) == 0 {
                        muls0 += 1;
                    } else {
                        muls += 1;
                    }
                    res.set_el(i, j, res.get_el(i, j) + self.get_el(i, k) * m.get_el(k, j))
                }
            }
        }
        sums -= (n*n) as u32;
        (res, OperationCount{sums, sums0, muls, muls0})
    }

    fn add(&self, m: &Matrix) -> (Matrix, OperationCount) {
        let mut sums = 0;
        let muls = 0;
        let mut sums0 = 0;
        let muls0 = 0;
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n {
                if self.get_el(i, j)== 0 || m.get_el(i, j) == 0 {
                    sums0 += 1;
                } else {
                    sums += 1;
                }
                res.set_el(i, j, self.get_el(i, j) + m.get_el(i, j))
            }
        }
        (res, OperationCount{sums, sums0, muls0, muls})
    }

    fn sub(&self, m: &Matrix) -> (Matrix, OperationCount) {
        let mut sums = 0;
        let muls = 0;
        let mut sums0 = 0;
        let muls0 = 0;
        let n = self.size.0;
        let mut res =  Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n {
                if self.get_el(i, j) == 0 || m.get_el(i, j) == 0 {
                    sums0 += 1;
                } else {
                    sums += 1;
                }
                res.set_el(i, j, self.get_el(i, j) - m.get_el(i, j))
            }
        }
        (res, OperationCount{sums, sums0, muls0, muls})
    }

    fn strassen_recursion(&self, m: &Matrix) -> (Matrix, OperationCount) {
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
            let mut ops = OperationCount{sums: 0, sums0:0, muls: 0, muls0:0};
            let mut ops1;
            (a_res, ops1) = a11.add(&a22);
            ops = ops.add(&ops1);
            (b_res, ops1) = b11.add(&b22);
            ops = ops.add(&ops1);
            let (p1, mut ops1) = a_res.strassen_recursion(&b_res);
            ops = ops.add(&ops1);

            (a_res, ops1) = a21.add(&a22);
            ops = ops.add(&ops1);
            let (p2, mut ops1) = a_res.strassen_recursion(&b11);
            ops = ops.add(&ops1);

            (b_res, ops1) = b12.sub(&b22);
            ops = ops.add(&ops1);
            let (p3, mut ops1) = a11.strassen_recursion(&b_res);
            ops = ops.add(&ops1);

            (b_res, ops1) = b21.sub(&b11);
            ops = ops.add(&ops1);
            let (p4, mut ops1) = a22.strassen_recursion(&b_res);
            ops = ops.add(&ops1);

            (a_res, ops1) = a11.add(&a12);
            ops = ops.add(&ops1);
            let (p5, mut ops1) = a_res.strassen_recursion(&b22);
            ops = ops.add(&ops1);

            (a_res, ops1) = a21.sub(&a11);
            ops = ops.add(&ops1);
            (b_res, ops1) = b11.add(&b12);
            ops = ops.add(&ops1);
            let (p6, mut ops1) = a_res.strassen_recursion(&b_res);
            ops = ops.add(&ops1);

            (a_res, ops1) = a12.sub(&a22);
            ops = ops.add(&ops1);
            (b_res, ops1) = b21.add(&b22);
            ops = ops.add(&ops1);
            let (p7, mut ops1) = a_res.strassen_recursion(&b_res);
            ops = ops.add(&ops1);

            let (c12, mut ops1) = p3.add(&p5);
            ops = ops.add(&ops1);
            let (c21, mut ops1) = p2.add(&p4);
            ops = ops.add(&ops1);

            (a_res, ops1) = p1.add(&p4);
            ops = ops.add(&ops1);
            (b_res, ops1) = a_res.add(&p7);
            ops = ops.add(&ops1);
            let (c11, mut ops1) = b_res.sub(&p5);
            ops = ops.add(&ops1);

            (a_res, ops1) = p1.add(&p3);
            ops = ops.add(&ops1);
            (b_res, ops1) = a_res.add(&p6);
            ops = ops.add(&ops1);
            let (c22, ops1) = b_res.sub(&p2);
            ops = ops.add(&ops1);

            let mut c = Matrix::zero_filled(n);
            for i in 0..new_n {
                for j in 0..new_n {
                    c.set_el(i, j, c11.get_el(i, j));
                    c.set_el(i, j + new_n, c12.get_el(i, j));
                    c.set_el(i + new_n, j, c21.get_el(i, j));
                    c.set_el(i + new_n, j + new_n, c22.get_el(i, j));
                }
            }

            return (c, ops)
        }
    }

    fn multiply_strassen (&self, m: &Matrix) -> (Matrix, OperationCount) {
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
        let (c_prep, ops) = a_prep.strassen_recursion(&b_prep);
        let mut c = Matrix::zero_filled(n);
        for i in 0..n {
            for j in 0..n { 
                c.set_el(i, j, c_prep.get_el(i, j));
            }
        }
        
        (c, ops)

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

    let mut n = 1;
    while n < 1500 {
        let m1 = Matrix::new(n);
        let m2 = Matrix::new(n);


        let (res1, ops1) = m1.multiply(&m2);
        println!("loops;{n};{};{};{};{}", ops1.muls, ops1.muls0, ops1.sums, ops1.sums0);

        let (res2, ops2) = m1.multiply_strassen(&m2);
        println!("strassen;{n};{};{};{};{}", ops2.muls, ops2.muls0, ops2.sums, ops2.sums0);

        assert!(res1.is_equal(&res2));


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
