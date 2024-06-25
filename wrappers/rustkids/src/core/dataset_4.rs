use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

trait Node {
    fn size(&self) -> usize;
    // 这里可以定义其他方法
}

struct Tensor<T> {
    data: Vec<T>,
}

impl<T> Node for Tensor<T> {
    fn size(&self) -> usize {
        self.data.len()
    }
    // 可以实现其他方法
}

struct DataSet {
    data: Rc<RefCell<HashMap<String, NodeRef>>>,
}

type NodeRef = Rc<RefCell<dyn Node>>;

impl Node for DataSet {
    fn size(&self) -> usize {
        self.data.borrow().len()
    }
    // 这里可以实现其他方法
}

impl DataSet {
    fn def<T>(&mut self, key: String, size: usize) -> NodeRef {
        let tensor: Tensor<T> = Tensor {
            data: vec![Default::default(); size],
        };
        let node_ref: NodeRef = Rc::new(RefCell::new(tensor));
        self.data.borrow_mut().insert(key, node_ref.clone());
        node_ref
    }
}

fn main() {
    let mut data_set = DataSet {
        data: Rc::new(RefCell::new(HashMap::new())),
    };
    let tensor_ref = data_set.def::<f64>("example_tensor".to_string(), 10);
    println!("Size of tensor: {}", tensor_ref.borrow().size());
}
