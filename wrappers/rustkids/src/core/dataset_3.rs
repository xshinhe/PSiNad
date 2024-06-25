use std::collections::HashMap;
use std::cell::RefCell;
use std::rc::Rc;

// 定义 Tensor trait
trait Tensor<T> {
    fn data(&self) -> Vec<T>;
}

// 实现 Tensor trait for Vec<T>
impl<T> Tensor<T> for Vec<T> {
    fn data(&self) -> Vec<T> {
        self.clone()
    }
}

// 定义 Node trait
trait Node {
    fn tensor<T>(&self) -> Rc<RefCell<dyn Tensor<T>>>;
}

// 实现 Node trait for Tensor<T>
impl<T> Node for Rc<RefCell<Vec<T>>> {
    fn tensor<U>(&self) -> Rc<RefCell<dyn Tensor<U>>> {
        self.clone() as Rc<RefCell<dyn Tensor<U>>>
    }
}

struct DataSet {
    data: Rc<RefCell<HashMap<String, NodeRef>>>,
}

// 定义 NodeRef 类型
type NodeRef = Rc<RefCell<dyn Node>>;

impl DataSet {
    pub fn new() -> Self {
        DataSet {
            data: Rc::new(RefCell::new(HashMap::new())),
        }
    }

    pub fn def<T: Tensor<T> + 'static>(&self, key: String, size: usize, info: String) -> Rc<RefCell<Vec<T>>> {
        // 创建一个大小为 size 的 Tensor<T>
        let tensor = vec![T::default(); size];
        let node = Rc::new(RefCell::new(tensor)); 
        self.data.borrow_mut().insert(key, node.clone());
        node
    }

    pub fn node(&self, key: &str) -> Option<NodeRef> {
        self.data.borrow().get(key).cloned()
    }

    pub fn at(&self, key: &str) -> Option<Rc<RefCell<DataSet>>> {
        self.node(key).and_then(|node| {
            if let Ok(data_set) = node.try_borrow().downcast_ref::<DataSet>() {
                Some(Rc::new(RefCell::new(data_set.clone())))
            } else {
                None
            }
        })
    }
}

fn main() {
    // 创建一个 DataSet
    let data_set = DataSet::new();

    // 定义一个大小为 5 的 Tensor<i32>，并将其添加到 DataSet 中
    let tensor_i32 = data_set.def("tensor_i32".to_string(), 5, "Info about tensor_i32".to_string());

    // 获取 tensor_i32 的数据并打印
    let data_i32 = tensor_i32.borrow().data();
    println!("Data of tensor_i32: {:?}", data_i32);

    // 定义一个大小为 3 的 Tensor<f64>，并将其添加到 DataSet 中
    let tensor_f64 = data_set.def("tensor_f64".to_string(), 3, "Info about tensor_f64".to_string());

    // 获取 tensor_f64 的数据并打印
    let data_f64 = tensor_f64.borrow().data();
    println!("Data of tensor_f64: {:?}", data_f64);
}
