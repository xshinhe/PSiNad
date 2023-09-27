#ifndef Kernel_H
#define Kernel_H

class Kernel {
    Kernel();
    int run();

   protected:
    std::vector<Kernel> _sub_kernerls;
};

#endif  // Kernel_H
