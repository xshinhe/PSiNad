#ifndef Func_Macro_H
#define Func_Macro_H
#include "types.h"

inline double REAL_OF(const double& val) { return val; }

inline double IMAG_OF(const double& val) { return 0.0f; }

inline double CONJ_OF(const double& val) { return val; }

inline double NORM_OF(const double& val) { return val * val; }

inline double ABS_OF(const double& val) { return std::abs(val); }

inline double PHAS_OF(const double& val) { return 1.0f; }

inline double REAL_OF(const num_complex& val) { return std::real(val); }

inline double IMAG_OF(const num_complex& val) { return std::imag(val); }

inline num_complex CONJ_OF(const num_complex& val) { return std::conj(val); }

inline double NORM_OF(const num_complex& val) { return std::norm(val); }

inline double ABS_OF(const num_complex& val) { return std::abs(val); }

inline num_complex PHAS_OF(const num_complex& val) { return val / std::abs(val); }

#endif  // Func_Macro_H
