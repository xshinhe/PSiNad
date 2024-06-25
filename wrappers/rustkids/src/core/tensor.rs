pub trait Size {
    fn size()->i32;
}

#[derive(Debug)]
struct Shape {
    dim: Vec<i32>
}

impl Shape {
    pub fn new(size: usize) -> Self {
        Shape {
            data: vec![0; size],
        }
    }
}
