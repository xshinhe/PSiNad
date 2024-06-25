use std::collections::HashMap;
use std::cell::RefCell;
use std::rc::Rc;

trait Node {
    // 定义 Node 的行为
}

struct DataSet {
    data: Rc<RefCell<HashMap<String, NodeRef>>>,
}

type NodeRef = Rc<RefCell<dyn Node>>;

impl DataSet {
    // 构造函数
    pub fn new() -> Self {
        DataSet {
            data: Rc::new(RefCell::new(HashMap::new())),
        }
    }

    // 定义一个 T 类型的变量
    pub fn def<T: Node + 'static>(&self, key: String, size: usize, info: String) -> Rc<RefCell<T>> {
        let data = vec![T::default(); size];
        let node = Rc::new(RefCell::new(data)); 
        self.data.borrow_mut().insert(key, node.clone());
        node
    }

    // 获取对应 key 的节点
    pub fn node(&self, key: &str) -> Option<NodeRef> {
        self.data.borrow().get(key).cloned()
    }

    // 获取对应 key 的数据集
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
    // 在这里编写测试代码
}
