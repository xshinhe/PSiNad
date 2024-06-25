// pub mod core;
// pub mod kernels;

// fn main() {
//     let stat = core::kernel::Status {data: 1};

//     let root = core::kernel::KernelNode::new(Box::new(
//         kernels::hello::Kernel::new())
//     );
//     let mut child1 = core::kernel::KernelNode::new(Box::new(
//         kernels::hello::Kernel::new())
//     );
//     let child2 = core::kernel::KernelNode::new(Box::new(
//         kernels::hello::Kernel::new())
//     );

//     child1.add_child(child2);
//     let mut root = root;
//     root.add_child(child1);
//     root.execute(&stat);
// }
// 

mod core {
    pub mod kernel;
}
mod kernels {
    pub mod hello;
}
mod solvers {
    pub mod nads;
}

use core::kernel::{Status, KernelNode};
use kernels::hello;
use std::rc::Rc;
use std::cell::RefCell;

fn main() {
    let stat = Status { data: 1 };
    let mut model = KernelNode::new(Rc::new(RefCell::new(hello::Kernel::new())));
    let solver = solvers::nads::builder(model);
    solver.borrow_mut().execute(&stat);
}

// fn main() {
//     let stat = Status { data: 1 };

//     let root = KernelNode::new(Box::new(hello::Kernel::new()));
//     let mut child1 = KernelNode::new(Box::new(hello::Kernel::new()));
//     let child2 = KernelNode::new(Box::new(hello::Kernel::new()));

//     child1.add_child(child2);
//     let mut root = root;
//     root.add_child(child1);
//     root.execute(&stat);
// }