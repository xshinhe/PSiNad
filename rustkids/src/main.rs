use std::path::Path;

gflags::define! {
    --dump = ""
}

gflags::define! {
    --load = ""
}

gflags::define! {
    //show some thing
    -h, --help = false
}

fn main() {
    let patterns = gflags::parse();

    if HELP.flag {
        gflags::print_help_and_exit(0);
    }

    {
        println!("searching for patterns given on command line: {:?}", patterns);
    }
}
