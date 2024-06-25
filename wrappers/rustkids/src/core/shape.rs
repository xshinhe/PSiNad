pub trait Size {
    fn size(&self) -> usize;
}

#[derive(Debug)]
struct Shape {
    dim: Vec<usize>
}

impl Shape {
    pub fn new(size: usize) -> Self {
        Shape {
            dim: vec![1; size], // Corrected field name to `dim`
        }
    }

    pub fn new_with_dim(dim_in: Vec<usize>) -> Self {
        Shape {
            dim: dim_in,
        }
    }
}

impl Size for Shape {
    fn size(&self) -> usize {
        self.dim.len()
    }
}

#[allow(unused_assignments)]

fn main() {
    let mut n:usize = 0;
    n = 100;
    let shape = Shape::new(n);
    n = 200;
    let shape2 = Shape::new_with_dim(vec![n,2,3,4]);
    println!("Shape size: {}", shape.size());
    println!("Shape size: {}", shape.size());
    println!("Shape size: {}", shape2.size());
    println!("Shape size: {:?}", shape2);
}