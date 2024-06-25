use crate::core::kernel;
use crate::core::kernel::*;

pub struct Kernel {
    private_data: u32,
}

impl Kernel {
    pub fn new() -> Self {
        Kernel {
            private_data: 0,
        }
    }

    pub fn get_private_data(&self) -> u32 {
        self.private_data
    }

    pub fn set_private_data(&mut self, new_data: u32) {
        self.private_data = new_data;
    }
}

impl KernelTrait for Kernel {
    fn initialize(&mut self, param: &ParamDict) {
        self.private_data = 10; // setting values
        println!("[init] Hello, world!");
    }

    fn execute(&self, stat: &Status) {
        println!("Hello, world!");
        println!("doing with stat: {:?}", stat);
    }

    fn finalize(&self) {
        println!("[final] Hello, world!");
    }
}
