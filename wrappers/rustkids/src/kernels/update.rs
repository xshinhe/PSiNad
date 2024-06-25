use crate::core::kernel;
use crate::core::kernel::*;

pub mod iteration {
    pub struct Kernel {
        niter: i32,
    }

    impl Kernel {
        pub fn new() -> Self {
            Kernel {
                niter: 10,
            }
        }
    }

    impl KernelTrait for Kernel {
        fn initialize(&mut self, param: &ParamDict) {
            //
        }

        fn execute(&self, stat: &Status) {
            println!("!!!!");
        }

        fn finalize(&self) {
            println!("[final] Hello, world!");
        }
    }   
}



pub mod for_x {
    pub struct Kernel {
        x: f32,
        p: f32,
        m: f32,
    }

    impl Kernel {
        pub fn new() -> Self {
            Kernel {
                x: 0,
                p: 0,
                f: 0,
            }
        }
    }

    impl KernelTrait for Kernel {
        fn initialize(&mut self, param: &ParamDict) {
            //
        }

        fn execute(&self, stat: &Status) {
            self.x += self.p / self.m;
            println!("x!");
        }

        fn finalize(&self) {
            println!("[final] Hello, world!");
        }
    }   
}

pub mod for_p {
    pub struct Kernel {
        x: f32,
        p: f32,
        m: f32,
    }

    impl Kernel {
        pub fn new() -> Self {
            Kernel {
                x: 0,
                p: 0,
                f: 0,
            }
        }
    }

    impl KernelTrait for Kernel {
        fn initialize(&mut self, param: &ParamDict) {
            //
        }

        fn execute(&self, stat: &Status) {
            self.x += self.p / self.m;
            println!("p!");
        }

        fn finalize(&self) {
            println!("[final] Hello, world!");
        }
    }   
}

pub mod for_u {
    pub struct Kernel {
        x: f32,
        p: f32,
        m: f32,
    }

    impl Kernel {
        pub fn new() -> Self {
            Kernel {
                x: 0,
                p: 0,
                f: 0,
            }
        }
    }

    impl KernelTrait for Kernel {
        fn initialize(&mut self, param: &ParamDict) {
            //
        }

        fn execute(&self, stat: &Status) {
            self.x += self.p / self.m;
            println!("u!");
        }

        fn finalize(&self) {
            println!("[final] Hello, world!");
        }
    }   
}

pub mod for_f {
    pub struct Kernel {
        x: f32,
        p: f32,
        m: f32,
    }

    impl Kernel {
        pub fn new() -> Self {
            Kernel {
                x: 0,
                p: 0,
                f: 0,
            }
        }
    }

    impl KernelTrait for Kernel {
        fn initialize(&mut self, param: &ParamDict) {
            //
        }

        fn execute(&self, stat: &Status) {
            self.x += self.p / self.m;
            println!("f!");
        }

        fn finalize(&self) {
            println!("[final] Hello, world!");
        }
    }   
}

