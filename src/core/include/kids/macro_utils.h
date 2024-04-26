#ifndef KIDS_MACRO_UTILS_H
#define KIDS_MACRO_UTILS_H

#if defined(_MSC_VER)
#define FUNCTION_NAME __FUNCSIG__
#elif defined(__INTEL_COMPILER)
#define FUNCTION_NAME __FUNCTION_SIGNATURE__
#else
#define FUNCTION_NAME __PRETTY_FUNCTION__
#endif

#endif  // KIDS_MACRO_UTILS_H
