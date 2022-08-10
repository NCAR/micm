// Copyright (C) 2022 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
//
#include <micm/micm.h>
#include <micm/system/system.h>
#include <iostream>

#if __cplusplus
extern "C" {
#endif

void create_solver(char* file_path) {
  std::cout << "Hi" << std::endl;
  System system;
}

#if __cplusplus
}
#endif
