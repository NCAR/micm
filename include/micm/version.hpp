// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  const char* getmicmVersion()
  {
    return "3.0.0";
  }
  unsigned getmicmVersionMajor()
  {
    return 3;
  }
  unsigned getmicmVersionMinor()
  {
    return 0+0;
  }
  unsigned getmicmVersionPatch()
  {
    return 0+0;
  }
  unsigned getmicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
