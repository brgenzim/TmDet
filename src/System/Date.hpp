// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <chrono>
#include <string>
#include <sstream>

/**
 * @brief namespace for tmdet system
 *
 * @namespace Tmdet
 * @namespace System
 */
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
