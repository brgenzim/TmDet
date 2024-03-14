#ifndef __UNITMP_PDBLIB_TYPES_PDB_CHAIN_TYPES__
#define __UNITMP_PDBLIB_TYPES_PDB_CHAIN_TYPES__

#include <unordered_map>

namespace UniTmp::PdbLib::Types {

    struct PdbChainType {
        const char *name;
        const char *description;
    };

    namespace PdbChainTypes {
        const PdbChainType PRO {"Protein", "Chain belonging to a protein in a PDB file"};
        const PdbChainType NUC {"Nucleic", "Chain belonging to a dna or rna in a PDB file"};
        const PdbChainType CAO {"Calpha only", "Chain belonging to a protein in a PDB file containing only C alpha atoms"};
        const PdbChainType BBO {"Backbone only", "Chain belonging to a protein in a PDB file containing only backbone atoms"};
        const PdbChainType HET {"Hetero atoms only", "Chain in a PDB file containing only hetero atoms"};
        const PdbChainType OTH {"Other", "Other, unspecified chain type"};
        const PdbChainType UNK {"Unknown", "Unknown chain type"};

        const std::unordered_map<const char*, PdbChainType> all {
            {"PRO",PRO},{"NUC",NUC},{"CAO",CAO},{"BBO",BBO},{"HET",HET},{"OTH",OTH},{"UNK",UNK}
        };
    }

}

#endif
