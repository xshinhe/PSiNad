use std::collections::HashMap;
use std::sync::{Arc, RwLock};

// Placeholder for Node trait
trait NodeTrait: Send + Sync {
    fn as_any(&self) -> &dyn std::any::Any;
}

// Define the Shape struct with a size attribute
struct Shape {
    size: usize,
}

// Define psnd_int, psnd_real, and psnd_complex structs
struct psnd_int {
    data: Vec<i32>,
}

struct psnd_real {
    data: Vec<f64>,
}

struct psnd_complex {
    data: Vec<(f64, f64)>, // Assuming complex numbers are represented as (real, imag) tuples
}

// Implement Default for psnd_int, psnd_real, and psnd_complex
impl Default for psnd_int {
    fn default() -> Self {
        psnd_int { data: Vec::new() }
    }
}

impl Default for psnd_real {
    fn default() -> Self {
        psnd_real { data: Vec::new() }
    }
}

impl Default for psnd_complex {
    fn default() -> Self {
        psnd_complex { data: Vec::new() }
    }
}

// Implement NodeTrait for psnd_int, psnd_real, and psnd_complex
impl NodeTrait for psnd_int {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

impl NodeTrait for psnd_real {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

impl NodeTrait for psnd_complex {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

// Define DataSet
pub struct DataSet {
    data: Arc<RwLock<HashMap<String, Arc<dyn NodeTrait>>>>,
}

impl DataSet {
    pub fn new() -> Self {
        DataSet {
            data: Arc::new(RwLock::new(HashMap::new())),
        }
    }

    pub fn def_int(&self, key: &str, shape: Shape, info: &str) -> Arc<psnd_int> {
        let var = Arc::new(psnd_int {
            data: vec![0; shape.size],
        });
        self.data.write().unwrap().insert(key.to_string(), var.clone());
        var
    }

    pub fn def_real(&self, key: &str, shape: Shape, info: &str) -> Arc<psnd_real> {
        let var = Arc::new(psnd_real {
            data: vec![0.0; shape.size],
        });
        self.data.write().unwrap().insert(key.to_string(), var.clone());
        var
    }

    pub fn def_complex(&self, key: &str, shape: Shape, info: &str) -> Arc<psnd_complex> {
        let var = Arc::new(psnd_complex {
            data: vec![(0.0, 0.0); shape.size],
        });
        self.data.write().unwrap().insert(key.to_string(), var.clone());
        var
    }

    pub fn node(&self, key: &str) -> Option<Arc<dyn NodeTrait>> {
        self.data.read().unwrap().get(key).cloned()
    }

    pub fn at(&self, key: &str) -> Option<Arc<DataSet>> {
        self.node(key).and_then(|node| node.as_any().downcast_ref::<DataSet>().cloned())
    }

    pub fn _def_int(&self, key: &str, shape: Shape, info: &str) -> &Self {
        self.def_int(key, shape, info);
        self
    }

    pub fn _def_real(&self, key: &str, shape: Shape, info: &str) -> &Self {
        self.def_real(key, shape, info);
        self
    }

    pub fn _def_complex(&self, key: &str, shape: Shape, info: &str) -> &Self {
        self.def_complex(key, shape, info);
        self
    }
}

impl Default for DataSet {
    fn default() -> Self {
        DataSet::new()
    }
}

fn main() {
    let dataset = DataSet::new();
    let shape = Shape { size: 10 };

    let int_var = dataset.def_int("int_var", shape, "Some info");
    let real_var = dataset.def_real("real_var", Shape { size: 5 }, "Some info");
    let complex_var = dataset.def_complex("complex_var", Shape { size: 20 }, "Some info");

    // Further operations...
}
