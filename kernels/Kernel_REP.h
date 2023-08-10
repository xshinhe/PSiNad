#ifndef Kernel_REP_H
#define Kernel_REP_H

DEFINE_POLICY(EOMPOCILY,
              DDDD,  // == Diabatic, Initial, Elec-EOMs, Nucl-EOMs, TCF
              DAAD,  // == Diabatic, Initial, Elec-EOMs, Nucl-EOMs, TCF
              AAAA,  // == Adiabatic
              DDAA   // ==
);

namespace PROJECT_NS {
class Kernel_REP final : public Kernel {
   public:
   private:
    void init_calc_impl(int stat = -1){

    };
};
};  // namespace PROJECT_NS


#endif  // Kernel_REP_H