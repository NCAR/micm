// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  const char* GetMicmVersion()
  {
    return "3.5.0";
  }
  unsigned GetMicmVersionMajor()
  {
    return 3;
  }
  unsigned GetMicmVersionMinor()
  {
    return 5+0;
  }
  unsigned GetMicmVersionPatch()
  {
    return 0+0;
  }
  unsigned GetMicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
