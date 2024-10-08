#pragma once

#include <chrono>
#include <string>
#include <sstream>

namespace Tmdet::System {

    /**
     * @brief Helper for getting date
     * 
     */
    class Date {

        public:
            
            /**
             * @brief get the curent date
             * 
             * @return std::string 
             */
            static std::string get() {
                const std::chrono::time_point now{std::chrono::system_clock::now()};
                const std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)};
                std::stringstream ss;
                ss << ymd;
                return ss.str();
            }
    };
}
