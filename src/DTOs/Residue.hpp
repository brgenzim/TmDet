#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for data transfer objects
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief reading and writing residue data and parsing
     *        it to Residue Value Object
     */
    struct Residue {

        /**
         * @brief get data from gemmi residue
         * 
         * @param residue 
         * @param chainIdx 
         * @param residueIdx 
         * @return Tmdet::VOs::Residue 
         */
        static Tmdet::VOs::Residue get(gemmi::Residue& residue, int chainIdx, int residueIdx);

        /**
         * @brief string representation of the residue
         * 
         * @param residueVO 
         */
        static std::string toString(const Tmdet::VOs::Residue& residueVO);
    };
}
