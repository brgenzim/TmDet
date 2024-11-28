#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Atom.hpp>

/**
 * @brief namespace for data transfer objects
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief reading and writing protein structure and parsing
     *        it to Protein Value Object
     */
    struct Protein {

        /**
         * @brief write transformed structure into file in cif format
         * 
         * @param protein 
         * @param path 
         */
        static void writeCif(const Tmdet::ValueObjects::Protein& protein, const std::string& path);

        /**
         * @brief get pdb structure and parse it into Protein Value Object
         * 
         * @param inputPath
         */
        static Tmdet::ValueObjects::Protein get(const std::string& inputPath);

        /**
         * @brief Unselect antibody chains
         *
         * @param protein
         */
        static void unselectAntiBodyChains(Tmdet::ValueObjects::Protein& protein);

        /**
         * @brief unselect chains
         * 
         * @param chainList
         * @param protein
         */
        static void unselectChains(const std::string& chainList, Tmdet::ValueObjects::Protein& protein);

        /**
         * @brief string representation of the protein
         * 
         * @param protein 
         * @return std::string 
         */
        static std::string toString(const Tmdet::ValueObjects::Protein& protein);
    };
}
