#pragma once

#include <string>
#include <format>

#define TMDET_VERSION_MAJOR 4
#define TMDET_VERSION_MINOR 0
#define TMDET_VERSION_PATCH 0

#define TMDET_VERSION_DEVELOP_MAJOR 0
#define TMDET_VERSION_DEVELOP_MINOR 6
#define TMDET_VERSION_DEVELOP_PATCH 0

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
        return std::format("{}.{}.{}-{}.{}.{}",
                    TMDET_VERSION_MAJOR, TMDET_VERSION_MINOR, TMDET_VERSION_PATCH,
                    TMDET_VERSION_DEVELOP_MAJOR, TMDET_VERSION_DEVELOP_MINOR, TMDET_VERSION_DEVELOP_PATCH);
    }
}