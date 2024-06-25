use std::rc::Rc;
use std::cell::RefCell;

pub struct KernelNode {
    node: Rc<RefCell<dyn KernelTrait>>,
    children: Vec<Rc<RefCell<KernelNode>>>,
}

#[derive(Debug)]
pub struct Status {
    pub data: u32,
}

#[derive(Debug)]
pub struct ParamDict {
    pub data: u32,
}

pub trait KernelTrait {
    fn initialize(&mut self, param: &ParamDict);
    fn execute(&self, stat: &Status);
    fn finalize(&self);
}

impl KernelNode {
    pub fn new(node: Rc<RefCell<dyn KernelTrait>>) -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(KernelNode {
            node,
            children: Vec::new(),
        }))
    }

    pub fn add_child(&mut self, child: Rc<RefCell<KernelNode>>) {
        self.children.push(child);
    }

    // Assuming implementation of execute method
    pub fn execute(&self, stat: &Status) {
        // Execute code here
        self.node.borrow_mut().execute(stat);
        for child in &self.children {
            child.borrow_mut().execute(stat);
        }
    }
}
