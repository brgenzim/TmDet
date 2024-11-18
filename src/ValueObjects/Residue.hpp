#pragma once

#include <string>
#include <vector>
#include <any>
#include <iostream>
#include <unordered_map>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Atom.hpp>
#include <ValueObjects/HBond.hpp>
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
         * @brief the gemmi residue object
         */
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
         * @brief type of the residue
         */
        Tmdet::Types::Residue type = Tmdet::Types::ResidueType::UNK;

        /**
         * @brief secondary structure of the residue
         */
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructType::U;

        /**
         * @brief lowest energy back chain hydrogen bond of the residue
         */
        HBond hbond1;

        /**
         * @brief second lowest energy back chain hydrogen bond of the residue
         */
        HBond hbond2;
        
        /**
         * @brief order of residues after alignment
         *        Use this to measure residues distance in the sequence
         */
        int idx;

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

        /**
         * @brief Construct a new Residue object
         * 
         * @param _gemmi 
         * @param chainVO 
         */
        explicit Residue(gemmi::Residue& _gemmi) : 
            gemmi(_gemmi) {
                if (Tmdet::Types::Residues.contains(_gemmi.name)) {
                    type = Tmdet::Types::Residues.at(_gemmi.name);
                }
        }

        /**
         * @brief pdb sequence number of the residue
         * 
         * @return int 
         */
        int resn() const;
        
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

        void setNumberOfAtoms();

        bool hasOnlyBackBoneAtoms() const;

        gemmi::Vec3 centre() const;

        const gemmi::Atom* getCa() const;
    };
}
