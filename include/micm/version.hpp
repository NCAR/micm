// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  const char* getMicmVersion()
  {
    return "3.4.0";
  }
  unsigned getMicmVersionMajor()
  {
    return 3;
  }
  unsigned getMicmVersionMinor()
  {
    return 4+0;
  }
  unsigned getMicmVersionPatch()
  {
    return 0+0;
  }
  unsigned getMicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
