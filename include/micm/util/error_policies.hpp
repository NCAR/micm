// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>

namespace micm
{
  // JSON-Parsing error policies
  class ThrowPolicy
  {
    class Exception : public std::exception
    {
     public:
      const char* msg_;

     public:
      Exception(const char* msg)
          : msg_(msg)
      {
      }

      virtual const char* what()
      {
        return msg_;
      }
    };

   public:
    void OnError(std::string message)
    {
      throw Exception(message.c_str());
    }
  };

  class LogToStandardErrorPolicy
  {
   public:
    void OnError(std::string message)
    {
      std::cerr << message << std::endl;
    }
  };

  class IgnoreErrorsPolicy
  {
   public:
    void OnError(std::string message)
    {
    }
  };

}