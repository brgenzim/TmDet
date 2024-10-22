#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <format>
#include <utility>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <source_location>

#define LOG_LEVEL_TRACE 1
#define LOG_LEVEL_DEBUG 2
#define LOG_LEVEL_INFO 3
#define LOG_LEVEL_WARN 4
#define LOG_LEVEL_ERROR 5
#define LOG_LEVEL_CRITICAL 6
#define LOG_LEVEL_OFF 7

#define LOG_LEVEL_TRACE_STR trace
#define LOG_LEVEL_DEBUG_STR debug
#define LOG_LEVEL_INFO_STR info
#define LOG_LEVEL_WARN_STR warn
#define LOG_LEVEL_ERROR_STR err
#define LOG_LEVEL_CRITICAL_STR critical
#define LOG_LEVEL_OFF_STR off

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_TRACE_STR
#define TRACE_LOG(...) logger.trace(__VA_ARGS__)
#else
#define TRACE_LOG(...)
#endif

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_DEBUG_STR
#define DEBUG_LOG(...) logger.debug(__VA_ARGS__)
#else
#define DEBUG_LOG(...)
#endif

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_INFO_STR
#define INFO_LOG(...) logger.info(__VA_ARGS__)
#else
#define INFO_LOG(...)
#endif

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_WARN_STR
#define WARN_LOG(...) logger.warn(__VA_ARGS__)
#else
#define WARN_LOG(...)
#endif

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_ERROR_STR
#define ERROR_LOG(...) logger.error(__VA_ARGS__)
#else
#define ERROR_LOG(...)
#endif

#if defined TMDET_LOG_LEVEL && TMDET_LOG_LEVEL == LOG_LEVEL_CRITICAL_STR
#define CRITICAL_LOG(...) logger.critical(__VA_ARGS__)
#else
#define CRITICAL_LOG(...)
#endif

namespace Tmdet::System {

    enum level : int {
        LOG_LEVEL_TRACE_STR = LOG_LEVEL_TRACE,
        LOG_LEVEL_DEBUG_STR = LOG_LEVEL_DEBUG,
        LOG_LEVEL_INFO_STR = LOG_LEVEL_INFO,
        LOG_LEVEL_WARN_STR = LOG_LEVEL_WARN,
        LOG_LEVEL_ERROR_STR = LOG_LEVEL_ERROR,
        LOG_LEVEL_CRITICAL_STR = LOG_LEVEL_CRITICAL,
        LOG_LEVEL_OFF_STR = LOG_LEVEL_OFF,
        n_levels
    };

    class Logger {
        private:
            std::vector<std::ostream*> logStreams;
            level logLevel = level::err;
            const std::vector<std::string> logLevels = {
                "nop", "trace", "debug", "info", "warning", "error", "critical", "off"
            };

            bool shouldLog(level lvl) const {
                return lvl >= logLevel;
            }

            template <typename... Args>
            void logIt(level lvl, std::format_string<Args...> fmt, Args &&... args) {
                std::string str = std::format(fmt,std::forward<Args>(args)...);
                std::string dt = getCurrentDateTime();
                for( std::ostream* os : logStreams) {
                    (*os) << dt << "[" << logLevels[lvl] << "] " << str << std::endl;
                }
            }

            template <typename... Args>
            void log(level lvl, std::format_string<Args...> fmt, Args &&... args) {
                if (shouldLog(lvl)) {
                    logIt(lvl, fmt, std::forward<Args>(args)...);
                }
            }

            std::string getCurrentDateTime() {
                auto now = std::chrono::system_clock::now();
                std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
                std::tm localTime = *std::localtime(&currentTime);
                auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()) % 1000000;
                std::ostringstream oss;
                oss << std::put_time(&localTime, "[%Y-%m-%d %H:%M:%S.");
                oss << std::setfill('0') << std::setw(6) << microseconds.count() << "] ";
                return oss.str();
            }

        public:
            void addStream(std::ostream& os) {
                logStreams.push_back(&os);
            }

            void setLevel(level log_level) {
                logLevel = log_level;
            }

            template <typename... Args>
            void trace(std::format_string<Args...> fmt, Args &&...args) {
                log(level::trace, fmt, std::forward<Args>(args)...);
            }

            template <typename... Args>
            void debug(std::format_string<Args...> fmt, Args &&...args) {
                log(level::debug, fmt, std::forward<Args>(args)...);
            }

            template <typename... Args>
            void info(std::format_string<Args...> fmt, Args &&...args) {
                log(level::info, fmt, std::forward<Args>(args)...);
            }

            template <typename... Args>
            void warn(std::format_string<Args...> fmt, Args &&...args) {
                log(level::warn, fmt, std::forward<Args>(args)...);
            }

            template <typename... Args>
            void error(std::format_string<Args...> fmt, Args &&...args) {
                log(level::err, fmt, std::forward<Args>(args)...);
            }

            template <typename... Args>
            void critical(std::format_string<Args...> fmt, Args &&...args) {
                log(level::critical, fmt, std::forward<Args>(args)...);
            }
    };
}
