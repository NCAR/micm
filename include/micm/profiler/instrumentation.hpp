/* Copyright (C) 2019-2022 Hazel Engine
 * Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

namespace micm
{

  using FloatingPointMicroseconds = std::chrono::duration<double, std::micro>;

  struct ProfileResult
  {
    std::string name_;
    FloatingPointMicroseconds start_;
    std::chrono::microseconds elapsed_time_;
    std::thread::id thread_id_;
  };

  struct InstrumentationSession
  {
    std::string name_;
  };

  class Instrumentor
  {
   public:
    Instrumentor(const Instrumentor&) = delete;
    Instrumentor(Instrumentor&&) = delete;

    void BeginSession(const std::string& name, const std::string& filepath = "results.json")
    {
      std::lock_guard lock(mutex_);
      if (current_session_)
      {
        std::cerr << "Instrumentor::BeginSession('" << name << "') when session '" << current_session_->name_
                  << "' already open." << std::endl;
        InternalEndSession();
      }

      output_stream_.open(filepath);

      if (output_stream_.is_open())
      {
        current_session_ = new InstrumentationSession({ name });
        WriteHeader();
      }
      else
      {
        std::cerr << "Instrumentor could not open results file '" << filepath << "'." << std::endl;
      }
    }

    void EndSession()
    {
      std::lock_guard lock(mutex_);
      InternalEndSession();
    }

    void WriteProfile(const ProfileResult& result)
    {
      std::stringstream json;

      json << std::setprecision(3) << std::fixed;
      json << ",{";
      json << "\"cat\":\"function\",";
      json << "\"dur\":" << (result.elapsed_time_.count()) << ',';
      json << "\"name\":\"" << result.name_ << "\",";
      json << "\"ph\":\"X\",";
      json << "\"pid\":0,";
      json << "\"tid\":\"" << result.thread_id_ << "\",";
      json << "\"ts\":" << result.start_.count();
      json << "}";

      std::lock_guard lock(mutex_);
      if (current_session_)
      {
        output_stream_ << json.str();
        output_stream_.flush();
      }
    }

    static Instrumentor& Get()
    {
      static Instrumentor instance;
      return instance;
    }

   private:
    Instrumentor()
        : current_session_(nullptr)
    {
    }

    ~Instrumentor()
    {
      EndSession();
    }

    void WriteHeader()
    {
      output_stream_ << "{\"otherData\": {},\"traceEvents\":[{}";
      output_stream_.flush();
    }

    void WriteFooter()
    {
      output_stream_ << "]}";
      output_stream_.flush();
    }

    void InternalEndSession()
    {
      if (current_session_)
      {
        WriteFooter();
        output_stream_.close();
        delete current_session_;
        current_session_ = nullptr;
      }
    }

    std::mutex mutex_;
    InstrumentationSession* current_session_;
    std::ofstream output_stream_;
  };

  class InstrumentationTimer
  {
   public:
    InstrumentationTimer(const char* name)
        : name_(name),
          stopped_(false)
    {
      start_time_point_ = std::chrono::steady_clock::now();
    }

    ~InstrumentationTimer()
    {
      if (!stopped_)
      {
        Stop();
      }
    }

    void Stop()
    {
      auto end_time_point = std::chrono::steady_clock::now();
      auto start = FloatingPointMicroseconds{ start_time_point_.time_since_epoch() };
      auto elapsed_time = std::chrono::time_point_cast<std::chrono::microseconds>(end_time_point).time_since_epoch() -
                          std::chrono::time_point_cast<std::chrono::microseconds>(start_time_point_).time_since_epoch();

      Instrumentor::Get().WriteProfile({ name_, start, elapsed_time, std::this_thread::get_id() });

      stopped_ = true;
    }

   private:
    const char* name_;
    std::chrono::time_point<std::chrono::steady_clock> start_time_point_;
    bool stopped_;
  };

  namespace InstrumentorUtils
  {
    template<size_t N>
    struct ChangeResult
    {
      char data_[N];
    };

    template<size_t N, size_t K>
    constexpr auto CleanupOutputString(const char (&expr)[N], const char (&remove)[K])
    {
      ChangeResult<N> result = {};

      size_t srcIndex = 0;
      size_t dstIndex = 0;
      while (srcIndex < N)
      {
        size_t matchIndex = 0;
        while (matchIndex < K - 1 && srcIndex + matchIndex < N - 1 && expr[srcIndex + matchIndex] == remove[matchIndex])
        {
          matchIndex++;
        }
        if (matchIndex == K - 1)
        {
          srcIndex += matchIndex;
        }
        result.data_[dstIndex++] = expr[srcIndex] == '"' ? '\'' : expr[srcIndex];
        srcIndex++;
      }

      return result;
    }
  }  // namespace InstrumentorUtils

}  // namespace micm

#if MICM_PROFILE
#  if defined(__GNUC__) || defined(__ICC)
#    define MICM_FUNC_SIG __func__
   // #define MICM_FUNC_SIG __PRETTY_FUNCTION__
#  elif (defined(__FUNCSIG__) || (_MSC_VER))
#    define MICM_FUNC_SIG __func__
   // #define MICM_FUNC_SIG __FUNCSIG__
#  else
#    define MICM_FUNC_SIG "MICM_FUNC_SIG unknown!"
#  endif

#  define MICM_PROFILE_BEGIN_SESSION(name, filepath) ::micm::Instrumentor::Get().BeginSession(name, filepath)
#  define MICM_PROFILE_END_SESSION()                 ::micm::Instrumentor::Get().EndSession()
#  define MICM_PROFILE_SCOPE_LINE2(name, line)                                                         \
    constexpr auto fixedName##line = ::micm::InstrumentorUtils::CleanupOutputString(name, "__cdecl "); \
    ::micm::InstrumentationTimer timer##line(fixedName##line.data_)
#  define MICM_PROFILE_SCOPE_LINE(name, line) MICM_PROFILE_SCOPE_LINE2(name, line)
#  define MICM_PROFILE_SCOPE(name)            MICM_PROFILE_SCOPE_LINE(name, __LINE__)
#  define MICM_PROFILE_FUNCTION()             MICM_PROFILE_SCOPE(MICM_FUNC_SIG)
#else
#  define MICM_PROFILE_BEGIN_SESSION(name, filepath)
#  define MICM_PROFILE_END_SESSION()
#  define MICM_PROFILE_SCOPE(name)
#  define MICM_PROFILE_FUNCTION()
#endif