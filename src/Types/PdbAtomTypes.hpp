#ifndef __UNITMP_PDBLIB_TYPES_PDB_ATOM_TYPES__
#define __UNITMP_PDBLIB_TYPES_PDB_ATOM_TYPES__

#include <unordered_map>

namespace UniTmp::PdbLib::Types {

    struct PdbAtomType {
        const char *symbol;
        const char *name;
        double vdw;
    };

    namespace PdbAtomTypes {
        const PdbAtomType AG {"AG","Silver",1.72};
        const PdbAtomType AL {"AL","Aluminium",0.675};
        const PdbAtomType AR {"AR","Argon",1.88};
        const PdbAtomType AS {"AS","Arsenic",1.85};
        const PdbAtomType AU {"AU","Gold",1.66};
        const PdbAtomType BR {"BR","Bromine",1.85};
        const PdbAtomType C  {"C","Carbon",1.76};
        const PdbAtomType C_ALI {"C_ALI","Alifase carbon",1.87};
        const PdbAtomType C_CAR {"C_CAR","Aromatic carbon",1.76};
        const PdbAtomType C_NUC{"C_NUC","Carbon in nucleic acids",1.80};
        const PdbAtomType CA{"CA","Calcium in of HETATM line",1.26};
        const PdbAtomType CD{"CD","Cadmium in of HETATM line",1.58};
        const PdbAtomType CL{"CL","Chlorine",1.75};
        const PdbAtomType CU{"CU","Copper",1.40};
        const PdbAtomType F{"F","Fluorine",1.47};
        const PdbAtomType FE{"FE","Iron",1.47};
        const PdbAtomType GA{"GA","Gallium",1.87};
        const PdbAtomType H{"H","Hydrogen",1.2};
        const PdbAtomType HE{"HE","Helium",1.40};
        const PdbAtomType HG{"HG","Mercury",1.55};
        const PdbAtomType I{"I","Iodine",1.98};
        const PdbAtomType IN{"IN","Indium",1.93};
        const PdbAtomType K{"K","Potassium",2.75};
        const PdbAtomType KR{"KR","Krypton",2.02};
        const PdbAtomType LI{"LI","Lithium",1.82};
        const PdbAtomType MG{"MG","Magnesium",1.73};
        const PdbAtomType MN{"MN","Manganese",0.81};
        const PdbAtomType N{"N","Nitrogen",1.65};
        const PdbAtomType N_AMN{"N_AMN","Amin nitrogen",1.50};
        const PdbAtomType N_AMD{"N_AMD","Amid nitrogen",1.65};
        const PdbAtomType N_NUC{"N_NUC","Nitrogen in nucleic acids",1.60};
        const PdbAtomType NA{"NA","Sodium",2.27};
        const PdbAtomType NI{"NI","Nickel",1.63};
        const PdbAtomType O{"O","Oxygen",1.40};
        const PdbAtomType O_CAR{"O_CAR","Delocalised carbonyl oxygen",1.50};
        const PdbAtomType P{"P","Phosphorus",1.90};
        const PdbAtomType PB{"PB","Lead",2.02};
        const PdbAtomType PD{"PD","Palladium",1.63};
        const PdbAtomType PT{"PT","Platinium",1.72};
        const PdbAtomType S{"S","Sulfur",1.85};
        const PdbAtomType SE{"SE","Selenium",1.90};
        const PdbAtomType SI{"SI","Silicon",2.10};
        const PdbAtomType SN{"SN","Tin",2.17};
        const PdbAtomType TE{"TE","Tellurium",2.06};
        const PdbAtomType TL{"TL","Thallium",1.96};
        const PdbAtomType U{"U","Uranium",1.86};
        const PdbAtomType XE{"XE","Xenon",2.16};
        const PdbAtomType ZN{"ZN","Zincz",1.39};
        const PdbAtomType UNK{"UNK","Unknown",1.80};

        const std::unordered_map<const char*, PdbAtomType> all {
            {"AG",AG},
            {"AL",AL},
            {"AR",AR},
            {"AS",AS},
            {"AU",AU},
            {"BR",BR},
            {"C",C},
            {"C_ALI",C_ALI},
            {"C_CAR",C_CAR},
            {"C_NUC",C_NUC},
            {"CA",CA},
            {"CD",CD},
            {"CL",CL},
            {"CU",CU},
            {"F",F},
            {"FE",FE},
            {"GA",GA},
            {"H",H},
            {"HE",HE},
            {"HG",HG},
            {"I",I},
            {"IN",IN},
            {"K",K},
            {"KR",KR},
            {"LI",LI},
            {"MG",MG},
            {"MN",MN},
            {"N",N},
            {"N_AMN",N_AMN},
            {"N_AMD",N_AMD},
            {"N_NUC",N_NUC},
            {"NA",NA},
            {"NI",NI},
            {"O",O},
            {"O_CAR",O_CAR},
            {"P",P},
            {"PB",PB},
            {"PD",PD},
            {"PT",PT},
            {"S",S},
            {"SE",SE},
            {"SI",SI},
            {"SN",SN},
            {"TE",TE},
            {"TL",TL},
            {"U",U},
            {"XE",XE},
            {"ZN",ZN},
            {"UNK",UNK},
        };
    }
}

#endif
