#pragma once

#include <string>
#include <vector>
#include <ValueObjects/Matrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief the assembly in the membrane
     */
    struct BioMatrix {

        /**
         * @brief transformation and chain ids
         */
        std::vector<Matrix> matrices;

        /**
         * @brief deleted chains (e.g. antibodies, nanobodies, biologically 
         *        not relevant assemblies, etc)
         */
        std::vector<std::string> deletedChainIds;
    };
}
