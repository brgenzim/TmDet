// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <any>
#include <unordered_map>
#include <gemmi/model.hpp>
#include <VOs/TMatrix.hpp>

/**
 * @brief namespace for value objects
 * @namespace Tmdet
 * @namespace VOs
 */
namespace Tmdet::VOs {

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
         * @brief outside water accessible surface of the atom
         */
        double outSurface = 0.0;

        /**
         * @brief atom index in gemmi residue atoms
         */
        int idx = 0;

        /**
         * @brief chain index in gemmi
         */
        int chainIdx = 0;

        /**
         * @brief residue index in gemmi
         */
        int residueIdx = 0;

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

        /**
         * @brief transform atom coordinates using transformation matrix tmatrix
         * 
         * @param tmatrix 
         */
        void transform(Tmdet::VOs::TMatrix& tmatrix) {
            tmatrix.transform(gemmi.pos);
        }
    };
}
