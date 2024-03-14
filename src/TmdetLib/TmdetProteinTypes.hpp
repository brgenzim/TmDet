#ifndef __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__
#define __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__

#include <unordered_map>
#include <iostream>

using namespace std;

namespace UniTmp::TmdetLib {

    struct TmdetProteinType {
        char *name;
        char *description;
        bool globular;
        bool membrane;
        bool alpha;
        TmdetProteinType(TmdetProteinType& other) : 
            name(other.name),
            description(other.description),
            globular(other.globular),
            membrane(other.membrane),
            alpha(other.alpha) {};
    };

    namespace TmdetProteinTypes {
            TmdetProteinType ALPHA = {
                "Tm_Alpha\0", 
                "Alpha helical transmembrane protein",
                false, true, true
            };
            TmdetProteinType BETA = {
                "Tm_Beta", 
                "Beta barrel transmembrane protein",
                false, true, false
            };
            TmdetProteinType GLOB = {
                "Soluble",
                "Globular, water soluble, not transmembrane protein",
                true, false, false
            };
            TmdetProteinType CA_TM = {
                "Ca_Tm", 
                "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                false, true, true 
            };
            TmdetProteinType CA_GLOB = {
                "Ca_Globular", 
                "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                true, false, false 
            };
            TmdetProteinType NUCL = {
                "Nucleotide", 
                "Nucleotide polimer chain",
                false, false, false 
            };
            TmdetProteinType MIXED = {
                "Tm_Mixed", 
                "Transmembrane protein containing both alpha helical membrane regions and beta barrel",
                false, true, true 
            };

            std::unordered_map<const char*, TmdetProteinType> all = {
                {"Tm_Alpha", ALPHA},
                {"Tm_Beta", BETA},
                {"Soluble", GLOB},
                {"Ca_Tm", CA_TM},
                {"Ca_Globular", CA_GLOB},
                {"Nucleotide", NUCL},
                {"Tm_Mixed", MIXED}
            };

           /*
            static TmdetProteinType get(const char* name) {
                TmdetProteinTypes t;
                std::unordered_map<const char*, TmdetProteinType> ta = t.all;
                TmdetProteinType type = ta[name];
                cerr << "+++" << type.name << "+++" << endl;
                return type;
            }*/
    };

}

#endif
