#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    struct Protein {
        std::string name;
        std::string description;
        bool globular;
        bool membrane;
        bool alpha;

        bool operator == (const Protein& other) const {
            return name == other.name;
        }

        bool isAlpha() const {
            return name == "Tm_Alpha" || name == "Tm_Mixed";
        }

        bool isBeta() const {
            return name == "Tm_Beta" || name == "Tm_Mixed";
        }

        bool isMixed() const {
            return name == "Tm_Mixed";
        }
    };

    namespace ProteinType {
        const Protein TM_ALPHA = { 
            "Tm_Alpha",
            "Alpha helical transmembrane protein",
            false, true, true
        };
        const Protein TM_BETA = { 
            "Tm_Beta",
            "Beta barrel transmembrane protein",
            false, true, false
        };
        const Protein TM_MIXED = { 
            "Tm_Mixed",
            "Transmembrane protein containing both alpha helical membrane regions and beta barrel",
            false, true, true 
        };
        const Protein SOLUBLE = {
            "Soluble",
            "Globular, water soluble, not transmembrane protein",
            true, false, false
        };
        const Protein CA_TM = { 
            "Ca_Tm",
            "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
            false, true, true 
        };
        const Protein CA_GLOBULAR = { 
            "Ca_Globular",
            "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
            true, false, false 
        };
        const Protein NOPROTEIN = { 
            "Noprotein",
            "PDB entry does not contain any amino acid polymer",
            false, false, false 
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
