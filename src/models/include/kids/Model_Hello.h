#ifndef MODEL_HELLO
#define MODEL_HELLO

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Model_Hello final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_Hello"; }

   private:
    int     N, F;
    double *x, *V, *dV, *ddV;
    double *m, *w;

    virtual void setInputParam_impl(std::shared_ptr<Param>& PM){};

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS){};

    virtual Status& executeKernel_impl(Status& stat) {
        std::cout << "Hello\n";
        return 0;
    }
};

};  // namespace PROJECT_NS

#endif  // MODEL_HELLO
