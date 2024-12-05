#pragma once

#include <string>
#include <vector>
#include <any>
#include <iostream>
#include <unordered_map>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Atom.hpp>
#include <gemmi/model.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief description of a Residue structure value objects
     */
    struct Residue {

        /**
         * @brief auth_id in cif file
         */
        int authId;

        /**
         * @brief label_id in cif file
         */
        int labelId;

        /**
         * @brief residue index in Tmdet ValueObject protein.chains[].residues
         */
        int idx;

        /**
         * @brief residue index in gemmi structure (protein.gemmi.chains[].residues)
         */
        int gemmiIdx = -1;

        /**
         * @brief insertion code
         */
        char authIcode = ' ';

        gemmi::Residue& gemmi;

        /**
         * @brief list of Atom structure value objects in the residue
         */
        std::vector<Atom> atoms;

        /**
         * @brief flag for selection
         */
        bool selected = true;

        /**
         * @brief solvent accessible surface of the residue
         */
        double surface = 0.0;

                /**
         * @brief outside water accessible surface of the residue
         */
        double outSurface = 0.0;

        /**
         * @brief type of the residue
         */
        Tmdet::Types::Residue type = Tmdet::Types::ResidueType::UNK;

        /**
         * @brief secondary structure of the residue
         */
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructType::U;

        /**
         * @brief number of backbone atoms
         */
        int nba = 0;

        /**
         * @brief number of side chain atoms
         */
        int nsa = 0;

        /**
         * @brief chain index of gemmi chain object
         * 
         */
        int chainIdx;

        /**
         * @brief secondary structure vector index
         */
        int secStrVecIdx = -1;

        /**
         * @brief temporary container for claculating various
         *        properties for the residue
         */
        std::unordered_map<std::string,std::any> temp;

        explicit Residue(gemmi::Residue& residue) :
            gemmi(residue) {}

        /**
         * @brief check if residue has all side chain atoms
         * 
         * @return true 
         * @return false 
         */
        bool hasAllSideChainAtoms() const;

        /**
         * @brief check if residue has all atoms
         * 
         * @return true 
         * @return false 
         */
        bool hasAllAtoms() const;

        /**
         * @brief Set the Number Of sidechain and backbone atoms 
         *        according to the residue type
         * 
         */
        void setProperties(std::string name);

        /**
         * @brief check if the residue has only backbone atoms
         * 
         * @return true 
         * @return false 
         */
        bool hasOnlyBackBoneAtoms() const;

        /**
         * @brief Get the C alpha or C beta or N or C atom
         * 
         * @return const gemmi::Atom* 
         */
        const gemmi::Atom* getCa() const;

        /**
         * @brief check if residue is a gap (i.e. has no any atoms in the 
         *         gemmi structure)
         */
        bool isGap() const;

        void transform(Tmdet::ValueObjects::TMatrix& tmatrix);

        bool isInside() const;
    };
}
