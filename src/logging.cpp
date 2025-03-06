#include "src/logging.h"
#include <chrono>

namespace Logging {
    auto start = std::chrono::system_clock::now();

#if !defined(NDEBUG)
    bool LOG_ENABLE_TRACE = true;
#else
    bool LOG_ENABLE_TRACE = false;
#endif

    /* Colors used in spdlog */
    // Formatting codes
    const std::string reset = "\033[m";
    const std::string bold = "\033[1m";
    const std::string dark = "\033[2m";
    const std::string underline = "\033[4m";
    const std::string blink = "\033[5m";
    const std::string reverse = "\033[7m";
    const std::string concealed = "\033[8m";
    const std::string clear_line = "\033[K";

    // Foreground colors
    const std::string black = "\033[30m";
    const std::string red = "\033[31m";
    const std::string green = "\033[32m";
    const std::string yellow = "\033[33m";
    const std::string blue = "\033[34m";
    const std::string magenta = "\033[35m";
    const std::string cyan = "\033[36m";
    const std::string white = "\033[37m";

    // Background colors
    const std::string on_black = "\033[40m";
    const std::string on_red = "\033[41m";
    const std::string on_green = "\033[42m";
    const std::string on_yellow = "\033[43m";
    const std::string on_blue = "\033[44m";
    const std::string on_magenta = "\033[45m";
    const std::string on_cyan = "\033[46m";
    const std::string on_white = "\033[47m";

    // Bold colors
    const std::string yellow_bold = "\033[33m\033[1m";
    const std::string red_bold = "\033[31m\033[1m";
    const std::string bold_on_red = "\033[1m\033[41m";

    const std::string color_trace = white;
    const std::string color_debug = cyan;
    const std::string color_info = green;
    const std::string color_warn = yellow_bold;
    const std::string color_err = red_bold;
    const std::string color_critical = bold_on_red;
    const std::string color_off = reset;


    void vreport_error(const char* format, fmt::format_args args) {
        const std::string prefix = "[" + color_err + "error" + color_off + "] ";
        auto end = std::chrono::system_clock::now();
        double et = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.;
        fmt::print("[{:>8.1f}] {}", et, prefix);
        fmt::vprint(format, args);
    }

    void vreport_warn(const char* format, fmt::format_args args) {
        const std::string prefix = "[" + color_warn + "warning" + color_off + "] ";
        auto end = std::chrono::system_clock::now();
        double et = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.;
        fmt::print("[{:>8.1f}] {}", et, prefix);
        fmt::vprint(format, args);
    }

    void vreport_info(const char* format, fmt::format_args args) {
        const std::string prefix = "[" + color_info + "info" + color_off + "] ";
        auto end = std::chrono::system_clock::now();
        double et = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.;
        fmt::print("[{:>8.1f}] {}", et, prefix);
        fmt::vprint(format, args);
    }

    void vreport_trace(const char* format, fmt::format_args args) {
        const std::string prefix = "[" + color_trace + "trace" + color_off + "] ";
        auto end = std::chrono::system_clock::now();
        double et = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.;
        fmt::print("[{:>8.1f}] {}", et, prefix);
        fmt::vprint(format, args);
    }
}
