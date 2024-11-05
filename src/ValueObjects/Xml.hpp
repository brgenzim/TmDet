#pragma once

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <Types/Protein.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/Region.hpp>
#include <ValueObjects/TMatrix.hpp>

namespace Tmdet::ValueObjects {

    struct XmlChain {

        /**
         * @brief chain identifier (auth_sym_id in cif)
         */
        std::string id;

        /**
         * @brief chain identifier (label_sym_id in cif)
         */
        std::string labId;

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
         * @brief annotated regions (i.e. side1, side2, membrane, loop etc)
         */
        std::vector<Region> regions;

        /**
         * @brief chain type (e.g. alpha, beta, globular etc)
         */
        Tmdet::Types::Chain type = Tmdet::Types::ChainType::UNK;
    };

    struct Xml {
        /**
         * @brief PDB code of the protein
         */
        std::string code = "";

        /**
         * @brief list of XmlChains
         */
        std::vector<XmlChain> chains;

        /**
         * @brief flag if the protein is transmembrane or not
         */
        bool tmp = false;

        /**
         * @brief date of creation of the entry
         */
        std::string date = "DD-MM-YYYY";

        /**
         * @brief tmdet version that was used for creating/modifying
         *        the entry
         */
        std::string version = "";

        /**
         * @brief description of modifications
         */
        std::vector<Modification> modifications;

        /**
         * @brief the final Q value that measure the fitness of the
         *        structure - membrane orientation
         * 
         */
        double qValue = 0.0;

        /**
         * @brief protein type (alpha, beta, mixed etc.)
         * 
         */
        Tmdet::Types::Protein type = Tmdet::Types::ProteinType::SOLUBLE;

        /**
         * @brief facts about the transmembrane type in the 
         *        corresponding swiss-prot (UniProt) entry
         *        it is not used from tmdet version 2.0
         */
        std::string spres = "";

        /**
         * @brief facts about the transmembrane type in the 
         *        corresponding PDB entry it is not used
         *        from tmdet version 2.0
         */
        std::string pdbkwres = "";

        /**
         * @brief biomatrix for generating the biological
         *        ensembly of chains
         */
        BioMatrix bioMatrix;

        /**
         * @brief description of membrane(s)
         */
        std::vector<Membrane> membranes;

        /**
         * @brief transformation matrix for the protein
         *        to be the membrane plane in the xy plane
         *        and the membrane normal parallel to the z axes
         */
        TMatrix tmatrix;

        void notTransmembrane();

    };
}