#pragma once

#include <string>
#include <VOs/Region.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    struct Region {

        /**
         * @brief print value object content of the region
         */
        static std::string toString(const Tmdet::VOs::Region& region);
    };
}
