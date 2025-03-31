// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <format>
#include <VOs/SecStrVec.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    struct SecStrVec {

        /**
         * @brief convert value object content into string
         */
        static std::string toString(const Tmdet::VOs::Protein& protein, const Tmdet::VOs::SecStrVec& vec) {
            return std::format(R"(
    SecStrVec type: {} sheet: {} barrel: {}
        begin: [{}, {}, {}] end: [{}, {}, {}]
        chainIdx: {} begResIdx: {} endResIdx: {}
)",
                vec.type.name, vec.sheetIdx, vec.barrelIdx,
                vec.begin.x,vec.begin.y,vec.begin.z,
                vec.end.x, vec.end.y, vec.end.z,
                protein.chains[vec.chainIdx].id,
                protein.chains[vec.chainIdx].residues[vec.begResIdx].authId,
                protein.chains[vec.chainIdx].residues[vec.endResIdx].authId);
        }
    };
}