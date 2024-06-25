use std::fmt;
use std::fmt::Write;
use std::sync::Arc;
use std::sync::Mutex;

type SizeType = usize;

pub trait Node {
    fn repr(&self) -> String;
    fn help(&self, name: &str) -> String;
}

#[derive(Clone)]
pub struct Shape {
    dimensions: Vec<SizeType>,
}

impl Shape {
    pub fn size(&self) -> SizeType {
        self.dimensions.iter().product()
    }
}

pub struct Tensor<T> {
    _size: SizeType,
    _shape: Shape,
    _doc_info: String,
    _data: Arc<Mutex<Vec<T>>>,
}

impl<T: Default + fmt::Display + Clone + 'static> Tensor<T> {
    pub fn new(shape: Shape, info: &str) -> Self {
        let size = shape.size();
        let data = vec![T::default(); size];
        Self {
            _size: size,
            _shape: shape,
            _doc_info: info.to_string(),
            _data: Arc::new(Mutex::new(data)),
        }
    }

    pub fn data(&self) -> Arc<Mutex<Vec<T>>> {
        Arc::clone(&self._data)
    }

    pub fn size(&self) -> SizeType {
        self._size
    }

    pub fn shape(&self) -> &Shape {
        &self._shape
    }
}

impl<T: Default + fmt::Display + Clone + 'static> Node for Tensor<T> {
    fn repr(&self) -> String {
        let data = self._data.lock().unwrap();
        let mut os = String::new();
        write!(os, "{}", std::any::type_name::<T>()).unwrap();
        write!(os, " Size: {}", self._size).unwrap();
        write!(os, "\n").unwrap();
        for value in data.iter() {
            write!(os, "        {}", value).unwrap();
        }
        os
    }

    fn help(&self, _name: &str) -> String {
        self._doc_info.clone()
    }
}

// Example usage

fn main() {
    let shape = Shape { dimensions: vec![2, 3] };
    let tensor: Tensor<i32> = Tensor::new(shape, "This is a tensor of integers");
    
    println!("{}", tensor.repr());
    println!("{}", tensor.help("Tensor"));
}