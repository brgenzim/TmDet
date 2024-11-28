#pragma once

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Region.hpp>
#include <ValueObjects/SecStrVec.hpp>
#include <gemmi/model.hpp>

/**
 * @brief namespace for value objects
 * @namespace Tmdet
 * @namespace ValueObjects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief chain in the protein
     */
    struct Chain {
        /**
         * @brief chain identifier (auth_sym_id in cif)
         */
        std::string id;

        /**
         * @brief chain label
         */
        std::string labelId;

        /**
         * @brief entity identifier in cif file (same as label_entity_id in ATOM lines)
         */
        std::string entityId;

        /**
         * @brief entity index in gemmi structrue entities vector
         */
        int entityIdx;

        /**
         * @brief flag for selection
         */
        bool selected = true;

        /**
         * @brief number of transmembrane segments in the chain
         */
        int numtm = 0;

        /**
         * @brief amino acid sequence of the chain in one letter code
         */
        std::string seq;

        /**
         * @brief list of Tmdet ValueObjects Residues in the chain
         */
        std::vector<Residue> residues;

        /**
         * @brief chain index in Tmdet ValueObject protein.chains
         */
        int idx = 0;

        /**
         * @brief length of the chain in residues
         */
        int length = 0;

        /**
         * @brief annotated regions (i.e. side1, side2, membrane, loop etc)
         */
        std::vector<Region> regions;

        /**
         * @brief chain type (e.g. alpha, beta, globular etc)
         */
        Tmdet::Types::Chain type = Tmdet::Types::ChainType::UNK;

        /**
         * @brief transform atom coordinates in the chain
         * 
         * @param tmatrix 
         */
        void transform(Tmdet::ValueObjects::TMatrix& tmatrix);

        template<typename T>
        void eachResidue(T func) {
            for(auto& residue: residues) {
                func(residue);
            }
        }

        template<typename T>
        void eachSelectedResidue(T func) {
            if (selected) {
                for(auto& residue: residues) {
                    if (residue.selected) {
                        func(residue);
                    }
                }
            }
        }

    };
}
