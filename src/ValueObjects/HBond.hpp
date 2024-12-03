#pragma once

namespace Tmdet::ValueObjects {

    /**
     * @brief simple representation of a hydrogen bond
     */
    struct HBond {

        /**
         * @brief calculated energy of the hydrogen bond
         */
        double energy = 1e30;

        /**
         * @brief chain index of the akceptor residue
         */
        int toChainIdx = -1;

        /**
         * @brief residue index of the akceptor
         */
        int toResIdx = -1;

        int fromResIdx = -1;
    };
}
