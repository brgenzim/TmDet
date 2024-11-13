#pragma once

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Region.hpp>
#include <gemmi/model.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief forward definition of Residue
     */
    struct Residue;
    
    /**
     * @brief chain in the protein
     */
    struct Chain {
        /**
         * @brief chain identifier (auth_sym_id in cif)
         */
        std::string id;

        /**
         * @brief chain identifier (label_sym_id in cif)
         */
        std::string labId;

        /**
         * @brief entity identifier in cif file
         */
        std::string entityId;

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
         * @brief gemmi chain
         */
        gemmi::Chain gemmi;

        /**
         * @brief list of Tmdet ValueObjects Residues in the chain
         */
        std::vector<Residue> residues;

        /**
         * @brief chain index in gemmi Model
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

        void addStructure(const gemmi::Chain& _gemmi);

        gemmi::Vec3 centre();
    };
}
