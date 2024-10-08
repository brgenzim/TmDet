#pragma once

#include <string>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief description of modifications
     */
    struct Modification {

        /**
         * @brief date of the modification
         */
        std::string date;

        /**
         * @brief short description of changes
         */
        std::string descr;
    };
}
