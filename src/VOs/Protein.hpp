// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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
#include <VOs/Chain.hpp>
#include <VOs/Modification.hpp>
#include <VOs/BioMatrix.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/SecStrVec.hpp>

#define EACH_SELECTED_CHAIN(protein) for( auto& chain: protein.chains) if (chain.selected)
/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {
    
    /**
     * @brief Protein contains both gemmi structure 
     *        description and all information about
     *        protein assembly with orientation to the membrane
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
         * @brief list of Tmdet VOs Chains in the protein
         */
        std::vector<Chain> chains;

        /**
         * @brief flag if the protein is transmembrane or not
         */
        bool tmp = false;

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

        /**
         * @brief contains names of polymer entities (_entity.pdbx_description);
         *        keys are entity names (_entity.id)
         */
        std::map<std::string, std::string> polymerNames;

        /**
         * @brief vectors constructed from secondary structrure elements
         */
        std::vector<Tmdet::VOs::SecStrVec> secStrVecs;

        /**
         * @brief force single membrane (in case of curved membrane and
         *        fragment analysis)
         */
        bool forceSingleMembrane = false;

        /**
         * @brief number of beta barrels in the protein
         */
        int numBarrels = 0;

        std::string inputFile;

        bool hasIdenticalChains = false;

        int modelIndex = 0;

        gemmi::Vec3 firstTranslation = gemmi::Vec3(0, 0, 0);

        /**
         * @brief set transmembrane to no and clear data
         */
        void notTransmembrane();

        /**
         * @brief clear/update the content of the tmdetVO
         */
        void clear();

        /**
        * @brief helper for fetching and parsing a new protein object from a pdb file
        *        and stroring the gemmi structure
        * 
        * @param inputPath 
        */
        void getStructure(const std::string& inputPath);

        /**
        * @brief helper for fetching and parsing a new protein object from a cif file
        *        and stroring the gemmi structure as well as the cif document
        * 
        * @param inputPath 
        */
        void getCifStructure(const std::string& inputPath);

        /**
        * @brief helper for fetching and parsing a new protein object from an ent file
        *        and stroring the gemmi structure as well as create a cif document
        * 
        * @param inputPath 
        */
        void getEntStructure(const std::string& inputPath);

        void setupPolymerNames();

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

        /**
         * @brief calculate the mass centre of the protein
         * 
         * @return gemmi::Vec3 
         */
        gemmi::Vec3 centre();

        /**
         * @brief transform protein's atom coordinates by tmatrix
         */
        void transform();

        template<typename T>
        void eachChain(T func) {
            for(auto& chain: chains) {
                func(chain);
            }
        }

        template<typename T>
        void eachSelectedChain(T func) {
            for(auto& chain: chains) {
                if (chain.selected) {
                    func(chain);
                }
            }
        }

        template<typename T>
        void eachResidue(T func) {
            for(auto& chain: chains) {
                for(auto& residue: chain.residues) {
                    func(residue);
                }
            }
        }

        template<typename T>
        void eachSelectedResidue(T func) {
            for(auto& chain: chains) {
                if (chain.selected) {
                    for(auto& residue: chain.residues) {
                        if (residue.selected) {
                            func(residue);
                        }
                    }
                }
            }
        }

        /**
         * @brief number of selected chains
         */
        int numberOfSelectedChains() {
            int ret = 0;
            for(const auto& chain: chains) {
                ret += (chain.selected?1:0);
            }
            return ret;
        }

        /**
         * @brief number of selected residues
         */
        int numberOfSelectedResidues() {
            int ret = 0;
            eachSelectedResidue(
                [&](Tmdet::VOs::Residue& residue) {
                    ret++;
                }
            );
            return ret;
        }
    };
}
