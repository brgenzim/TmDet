#pragma once

#include <string>
#include <vector>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Membrane.hpp>


/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief Protein contains both gemmi structure 
     *        description and all information about
     *        protein ensembly with orientation to the membrane
     */
    struct Protein {

        /**
         * @brief PDB code of the protein
         */
        std::string code = "";

        /**
         * @brief the gemmi structure
         */
        gemmi::Structure gemmi;

        /**
         * @brief cif document description of the protein
         */
        gemmi::cif::Document document;

        /**
         * @brief gemmi neighbor search object for fast determination
         *        nearby atoms
         */
        gemmi::NeighborSearch neighbors;

        /**
         * @brief list of Tmdet ValueObjects Chains in the protein
         */
        std::vector<Chain> chains;

        /**
         * @brief flag if the protein is transmembrane or not
         */
        bool tmp;

        /**
         * @brief date of creation of the entry
         */
        std::string date;

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
        double qValue;

        /**
         * @brief protein type (alpha, beta, mixed etc.)
         * 
         */
        Tmdet::Types::Protein type;

        /**
         * @brief facts about the transmembrane type in the 
         *        corresponding swiss-prot (UniProt) entry
         *        it is not used from tmdet version 2.0
         */
        std::string spres;

        /**
         * @brief facts about the transmembrane type in the 
         *        corresponding PDB entry it is not used
         *        from tmdet version 2.0
         */
        std::string pdbkwres;

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

        /**
         * @brief set transmembrane to no and clear data
         */
        void notTransmembrane();

        /**
         * @brief clear/update the content of the tmdetVO
         */
        void clear();
    
        /**
        * @brief helper for fetching and parsing a new protein object from a cif file
        *        and stroring the gemmi structure as well as the cif document
        * 
        * @param inputPath 
        */
        void getStructure(const std::string& inputPath);

        /**
         * @brief Create a hash for structure that unique for each structure
         */
        std::string hash() const;

        /**
         * @brief unselect all chains
         */
        void unSelectAll();

        /**
         * @brief search chain by id (auth_asym_id)
         */
        int searchChainById(const std::string& id) const;

        /**
         * @brief search chain by id (label_asym_id)
         */
        int searchChainByLabId(const std::string& id) const;

        gemmi::Vec3 centre();
    };
}
