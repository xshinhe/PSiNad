use crate::core::kernel::Status;
use crate::core::kernel::KernelNode;
use crate::kernels;
use crate::kernels::update::*;

use std::rc::Rc;
use std::cell::RefCell;

macro_rules! new_rc_refcell {
    ($x:expr) => { Rc::new(RefCell::new($x)) };
}

macro_rules! new_kernel_node_rc_refcell {
    ($x:expr) => {
        KernelNode::new(new_rc_refcell!($x))
    };
}

pub fn builder (model: &KernelNode)->KernelNode
{
    let k_root = new_kernel_node_rc_refcell!(hello::Kernel::new());
    let mut k_iteration = new_kernel_node_rc_refcell!(update::iteratioin::Kernel::new());
    let mut k_update_x = new_kernel_node_rc_refcell!(update::for_x::Kernel::new());
    let mut k_update_p = new_kernel_node_rc_refcell!(update::for_p::Kernel::new());
    let mut k_update_u = new_kernel_node_rc_refcell!(update::for_u::Kernel::new());
    let mut k_update_f = new_kernel_node_rc_refcell!(update::for_f::Kernel::new());
    // let mut k_record = new_kernel_node_rc_refcell!(record::Kernel::new());

    // k_iteration.borrow_mut().add_child(k_record.clone());
    k_iteration.borrow_mut().add_child(k_update_p.clone());
    k_iteration.borrow_mut().add_child(k_update_x.clone());
    // k_iteration.borrow_mut().add_child(k_update_s.clone());
    k_iteration.borrow_mut().add_child(k_update_x.clone());
    k_iteration.borrow_mut().add_child(model.clone());
    k_iteration.borrow_mut().add_child(k_update_u.clone());
    k_iteration.borrow_mut().add_child(k_update_f.clone());
    k_iteration.borrow_mut().add_child(k_update_p.clone());

    // k_root.borrow_mut().add_child(k_load.clone()); // should set dim here?
    // k_root.borrow_mut().add_child(k_set_dim.clone());
    k_root.borrow_mut().add_child(k_iteration.clone());
    // k_root.borrow_mut().add_child(k_dump.clone());
    k_root
}
