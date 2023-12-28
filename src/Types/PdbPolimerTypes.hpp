#ifndef __UNITMP_PDBLIB_TYPES_PDB_POLIMER_TYPES__
#define __UNITMP_PDBLIB_TYPES_PDB_POLIMER_TYPES__

#include <unordered_map>

namespace UniTmp::PdbLib::Types {

    struct PdbPolimerType {
        const char *code;
        const char *name;
    };

    namespace PdbPolimerTypes {
        const PdbPolimerType CyclicPeptide {"CyclicPeptide","cyclic-pseudo-peptide"};
        const PdbPolimerType Other {"Other","other"};
        const PdbPolimerType PeptideNA {"PeptideNA","peptide nucleic acid"};
        const PdbPolimerType DNA {"DNA", "polydeoxyribonucleotide"};
        const PdbPolimerType DRNA {"DRNA", "polydeoxyribonucleotide/polyribonucleotide hybrid"};
        const PdbPolimerType DProtein {"DProtein", "polypeptide(D)"};
        const PdbPolimerType Protein {"Protein", "polypeptide(L)"};
        const PdbPolimerType RNA {"RNA", "polyribonucleotide"};
        
        const std::unordered_map<const char*, PdbPolimerType> all {
            {"CyclicPeptide",CyclicPeptide},
            {"Other", Other},
            {"PeptideNA",PeptideNA},
            {"DNA",DNA},
            {"DRNA", DRNA},
            {"DProtein", DProtein},
            {"Protein", Protein},
            {"RNA", RNA},
        };

        const std::unordered_map<const char*, PdbPolimerType> fromName {
            {"cyclic-pseudo-peptide",CyclicPeptide},
            {"other", Other},
            {"peptide nucleic acid",PeptideNA},
            {"polydeoxyribonucleotide",DNA},
            {"polydeoxyribonucleotide/polyribonucleotide hybrid", DRNA},
            {"polypeptide(D)", DProtein},
            {"polypeptide(L)", Protein},
            {"polyribonucleotide", RNA},
        };
    }

}

#endif
