// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <format>

#define TMDET_VERSION_MAJOR 4
#define TMDET_VERSION_MINOR 1
#define TMDET_VERSION_PATCH 4

/**
 * @namespace Tmdet
 */
namespace Tmdet {

    /**
    * @brief semver version of tmdet
    * 
    * @return std::string 
    */
    static std::string version() {
        return std::format("{}.{}.{}",
                    TMDET_VERSION_MAJOR, TMDET_VERSION_MINOR, TMDET_VERSION_PATCH);
    }
}