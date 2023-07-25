// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>

namespace micm
{
  // JSON-Parsing error policies
  template<class Object>
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
    Object OnError(std::string message)
    {
      throw Exception(message.c_str());
    }
  };

  template<class Object>
  class NoThrowPolicy
  {
   public:
    Object OnError(std::string message)
    {
      std::cerr << message << std::endl;
      return Object();
    }
  };

}