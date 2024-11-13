#pragma once

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <gemmi/model.hpp>

/**
 * @brief namespace for value objects
 * @namespace Tmdet
 * @namespace ValueObjects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief description of an atom
     */
    struct Atom {
        /**
         * @brief gemmi Atom
         */
        gemmi::Atom& gemmi;

        /**
         * @brief water accessible surface of the atom in the structure
         */
        double surface = 0.0;

        /**
         * @brief atom index in gemmi residue atoms
         */
        int idx;

        /**
         * @brief chain index in gemmi
         */
        int chainIdx;

        /**
         * @brief residue index in gemmi
         */
        int residueIdx;

        /**
         * @brief temporary container for claculating various
         *        properties for the atom
         */
        std::unordered_map<std::string, std::any> temp;

        /**
         * @brief Construct a new Atom object
         * 
         * @param _gemmi 
         * @param residueVO 
         * @param chainVO 
         */
        explicit Atom(gemmi::Atom& _gemmi) :
            gemmi(_gemmi) {}

    };
}
