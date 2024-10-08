#pragma once

#include <string>
#include <vector>
#include <any>
#include <iostream>
#include <unordered_map>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Atom.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/HBond.hpp>
#include <gemmi/model.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief forward declaration of Chain
     */
    struct Chain;

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
         * @brief solvent accessible surface of the residue
         */
        double surface;

        /**
         * @brief type of the residue
         */
        Tmdet::Types::Residue type;

        /**
         * @brief secondary structure of the residue
         */
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructs.at('-');

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
         * @brief chain index of gemmi chain object
         * 
         */
        int chainIdx;

        /**
         * @brief Chain structure value object the the residue belongs to
         */
        Chain& parentChain;

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
        explicit Residue(gemmi::Residue& _gemmi, Chain& chainVO) : 
            gemmi(_gemmi),
            parentChain(chainVO) {
                if (Tmdet::Types::Residues.contains(_gemmi.name)) {
                    type = Tmdet::Types::Residues.at(_gemmi.name);
                }
        }

        /**
         * @brief Destroy the Residue object
         */
        ~Residue()=default;

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
    };
}
