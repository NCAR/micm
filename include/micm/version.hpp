// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  const char* getMicmVersion()
  {
    return "3.3.1";
  }
  unsigned getMicmVersionMajor()
  {
    return 3;
  }
  unsigned getMicmVersionMinor()
  {
    return 3+0;
  }
  unsigned getMicmVersionPatch()
  {
    return 1+0;
  }
  unsigned getMicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
