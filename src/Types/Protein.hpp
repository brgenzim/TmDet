// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of tmdet types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    /**
     * @brief definition of protein type
     */
    struct Protein {
        /**
         * @brief name of the protein type
         */
        std::string name;

        /**
         * @brief description of protein type
         */
        std::string description;
        
        /**
         * @brief check if two protein types are equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator == (const Protein& other) const {
            return name == other.name;
        }

        /**
         * @brief check if protein type is alpha
         * 
         * @return true 
         * @return false 
         */
        bool isAlpha() const {
            return name == "Tm_Alpha" || name == "Tm_Mixed";
        }

        /**
         * @brief check if protein type is beta
         * 
         * @return true 
         * @return false 
         */
        bool isBeta() const {
            return name == "Tm_Beta" || name == "Tm_Mixed";
        }

        /**
         * @brief check if protein type is mixed
         * 
         * @return true 
         * @return false 
         */
        bool isMixed() const {
            return name == "Tm_Mixed";
        }
    };

    namespace ProteinType {
        const Protein TM_ALPHA = { 
            "Tm_Alpha",
            "Alpha helical transmembrane protein"
        };
        const Protein TM_BETA = { 
            "Tm_Beta",
            "Beta barrel transmembrane protein"
        };
        const Protein TM_MIXED = { 
            "Tm_Mixed",
            "Transmembrane protein containing both alpha helical membrane regions and beta barrel"
        };
        const Protein SOLUBLE = {
            "Soluble",
            "Globular, water soluble, not transmembrane protein"
        };
        const Protein CA_TM = { 
            "Ca_Tm",
            "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms"
        };
        const Protein CA_GLOBULAR = { 
            "Ca_Globular",
            "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms"
        };
        const Protein NOPROTEIN = { 
            "Noprotein",
            "PDB entry does not contain any amino acid polymer"
        };
    }
    
    const std::unordered_map<std::string, const Protein> Proteins = {
        { "Tm_Alpha", ProteinType::TM_ALPHA },
        { "Tm_Beta", ProteinType::TM_BETA },
        { "Tm_Mixed", ProteinType::TM_MIXED},
        { "Soluble", ProteinType::SOLUBLE},
        { "Ca_Tm", ProteinType::CA_TM},
        { "Ca_Globular", ProteinType::CA_GLOBULAR},
        { "Noprotein", ProteinType::NOPROTEIN}
    };
}
