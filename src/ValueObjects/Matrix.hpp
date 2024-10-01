#pragma once

#include <string>
#include <ValueObjects/TMatrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief transformation matrix and source chain identifiers
     *        for generating Bio Matrix
     */
    struct Matrix {

        /**
         * @brief identifier of the matrix
         */
        int id;

        /**
         * @brief chain identifier of the chain copied
         */
        std::string sourceChainId;

        /**
         * @brief chain identifier of the new chain
         */
        std::string newChainId;

        /**
         * @brief transformation matrix for generating coordinates of
         *        new chain from the source chain
         */
        TMatrix tmatrix;
    };
}
