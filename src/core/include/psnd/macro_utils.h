#ifndef PSND_MACRO_UTILS_H
#define PSND_MACRO_UTILS_H

#if defined(_MSC_VER)
#define FUNCTION_NAME __FUNCSIG__
#elif defined(__INTEL_COMPILER)
#define FUNCTION_NAME __FUNCTION_SIGNATURE__
#else
#define FUNCTION_NAME __PRETTY_FUNCTION__
#endif

#define NEW(X)

#endif  // PSND_MACRO_UTILS_H
