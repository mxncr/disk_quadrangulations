#pragma once

#include <iostream>
#include <cstdint>
#include <typeinfo>

#include "fmt/format.h"
#include "fmt/ostream.h"
#include "third_party/prettyprint.hpp"

#define FCT_TRACE() trace("{}() ...",std::string(__func__))
#define FCT_TRACE_ARGS(...) trace_fct_call(__func__, #__VA_ARGS__, __VA_ARGS__)
#define DBG(...) fprintf(stdout, "(DBG) %s:%i: ", __FILE__,__LINE__); Logging::show_name_value(std::cout, #__VA_ARGS__, __VA_ARGS__); printf("\n"); fflush(stdout)

namespace Logging {
  /* Logging usage:
   *   - formatting based on fmtlib (i.e. std::format in C++20)
   *   - example: 
   *        info("values are {:.2f} and {}", 0.3, "foo bar");
   *   - for the trace() function:
   *      - disabled at compile time if NDEBUG is defined
   *      - can be enabled/disable at run-time via the global variable LOG_ENABLE_TRACE
   *   - FCT_TRACE() macro can be used to trace function calls, and function arguments values 
   *     by passing them to the macro. See trace() to enable/disable.
   */

  void vreport_error(const char* format, fmt::format_args args);
  void vreport_warn(const char* format, fmt::format_args args);
  void vreport_info(const char* format, fmt::format_args args);
  void vreport_trace(const char* format, fmt::format_args args);
  extern bool LOG_ENABLE_TRACE;

  template <typename... Args>
    void error(const char* format, const Args & ... args) {
      vreport_error(format, fmt::make_format_args(args...));
      fmt::print("\n");
      fflush(stdout);
    }

  template <typename... Args>
    void warn(const char* format, const Args & ... args) {
      vreport_warn(format, fmt::make_format_args(args...));
      fmt::print("\n");
      fflush(stdout);
    }

  template <typename... Args>
    void info(const char* format, const Args & ... args) {
      vreport_info(format, fmt::make_format_args(args...));
      fmt::print("\n");
    }

  template <typename... Args>
    void trace(const char* format, const Args & ... args) {
// #if !defined(NDEBUG)
      if (LOG_ENABLE_TRACE) {
        vreport_trace(format, fmt::make_format_args(args...));
        fmt::print("\n");
      }
// #endif
    }

  // void trace_fct_call(const char* functionName) {
  //   vreport_trace("{} ...\n", functionName);
  // }

  template<typename H1>
    std::ostream& show_name_value(std::ostream& out, const char* label, H1&& value) {
      return out << label << "=" << std::forward<H1>(value);
    }

  template<typename H1, typename ...T>
    std::ostream& show_name_value(std::ostream& out, const char* label, H1&& value, T&&... rest) {
      const char* pcomma = strchr(label, ',');
      return show_name_value(out.write(label, pcomma - label) << "="
          << std::forward<H1>(value)
          << ',',
          pcomma + 1,
          std::forward<T>(rest)...);
    }

  template <typename... Args>
    void trace_fct_call(const char* functionName, const char* format, const Args & ... args) {
      vreport_trace("{}(", fmt::make_format_args(std::string(functionName)));
      show_name_value(std::cout,format,args...);
      // fmt::vprint(format, fmt::make_format_args(args...));
      fmt::print(") ...\n");
    }

  // std::ostream& show_name_value(std::ostream& out) {
  //   return out;
  // }
}

