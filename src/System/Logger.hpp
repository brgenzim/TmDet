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

#define LOG_LEVEL_TRACE 1
#define LOG_LEVEL_DEBUG 2
#define LOG_LEVEL_INFO 3
#define LOG_LEVEL_WARN 4
#define LOG_LEVEL_ERROR 5
#define LOG_LEVEL_CRITICAL 6
#define LOG_LEVEL_OFF 7

namespace Tmdet::System {

    enum level : int {
        trace = LOG_LEVEL_TRACE,
        debug = LOG_LEVEL_DEBUG,
        info = LOG_LEVEL_INFO,
        warn = LOG_LEVEL_WARN,
        err = LOG_LEVEL_ERROR,
        critical = LOG_LEVEL_CRITICAL,
        off = LOG_LEVEL_OFF,
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
            void logIt(std::format_string<Args...> fmt, Args &&... args) {
                std::string str = std::format(fmt,std::forward<Args>(args)...);
                std::string dt = getCurrentDateTime();
                for( std::ostream* os : logStreams) {
                    (*os) << dt << "[" << logLevels[logLevel] << "] " << str << std::endl;
                }
            }

            template <typename... Args>
            void log(level lvl, std::format_string<Args...> fmt, Args &&... args) {
                if (shouldLog(lvl)) {
                    logIt(fmt, std::forward<Args>(args)...);
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