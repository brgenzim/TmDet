#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Atom.hpp>

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
        static void writeCif(Tmdet::VOs::Protein& protein, const std::string& path);

        /**
         * @brief get pdb structure and parse it into Protein Value Object
         * 
         * @param inputPath
         */
        static Tmdet::VOs::Protein get(const std::string& inputPath);

        /**
         * @brief Unselect antibody chains
         *
         * @param protein
         */
        static void unselectAntiBodyChains(Tmdet::VOs::Protein& protein);

        /**
         * @brief unselect chains
         * 
         * @param chainList
         * @param protein
         */
        static void unselectChains(const std::string& chainList, Tmdet::VOs::Protein& protein);

        /**
         * @brief string representation of the protein
         * 
         * @param protein 
         * @return std::string 
         */
        static std::string toString(const Tmdet::VOs::Protein& protein);

        /**
         * @brief add silver atoms to the boundary of membrane
         */
        static std::vector<gemmi::Vec3> addMembraneAtoms(Tmdet::VOs::Protein& protein);

        static void addPlaneMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>&ret);
        static void addBlendedMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret);
    };
}
