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
         * @brief print the TmdetStruct Value Object into an out stream
         * 
         * @param os 
         * @param protein 
         */
        static void print(std::ostream& os, const Tmdet::ValueObjects::Protein& protein);
        
        /**
         * @brief print the Tmdet Chain Value Object into an out stream
         * 
         * @param os 
         * @param chainVO 
         */
        static void printChain(std::ostream& os, const Tmdet::ValueObjects::Chain& chainVO);
        
        /**
         * @brief print the Tmdet Residue Value Object into an out stream
         * 
         * @param os 
         * @param residueVO 
         */
        static void printResidue(std::ostream& os, const Tmdet::ValueObjects::Residue& residueVO);
        
        /**
         * @brief print the Tmdet Atom Value Object into an out stream
         * 
         * @param os 
         * @param atomVO 
         */
        static void printAtom(std::ostream& os, const Tmdet::ValueObjects::Atom& atomVO);

        /**
         * @brief Get the amino acid sequence of a chain
         * 
         * @param protein 
         * @param chainVO 
         * @return std::vector<std::string> 
         */
        static std::vector<std::string> getChainSequence(const Tmdet::ValueObjects::Protein& protein,
            const gemmi::Chain& chainVO);

        static void transform(Tmdet::ValueObjects::Protein& protein);

        /**
         * @brief Unselect polymer chains based on their name and the given string set in TMDET_POLYMER_FILTER_FILE.
         */
        static void unselectPolymers(Tmdet::ValueObjects::Protein& protein);
    };
}
