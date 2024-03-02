#ifndef MODEL_HELLO
#define MODEL_HELLO

#include "../core/Kernel.h"

namespace kids {

class Model_Hello final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_Hello"; }

   private:
    int N, F;
    double *x, *V, *dV, *ddV;
    double *m, *w;

    virtual void read_param_impl(Param* PM){};

    virtual void init_data_impl(DataSet* DS){};

    virtual int exec_kernel_impl(int stat = -1) {
        std::cout << "Hello\n";
        return 0;
    }
};

};  // namespace kids

#endif  // MODEL_HELLO
