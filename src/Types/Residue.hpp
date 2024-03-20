#ifndef __TMDET_TYPES_RESIDUE__
#define __TMDET_TYPES_RESIDUE__

#include <unordered_map>
#include <string>
#include <Types/Atom.hpp>

using namespace std;

namespace Tmdet::Types {

    struct Residue {
        string name;
        char a1;
        unordered_map<string,Atom> atoms;
    };

    namespace ResidueType {
        const Residue ALA = {
            "ALA", 'A',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue CYS = {
            "CYS", 'C',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"SG", AtomType::S},
                {"OXT", AtomType::O}
            }
        };
        const Residue ASP = {
            "ASP", 'D',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"OD1", AtomType::O},
                {"OD2", AtomType::O},
                {"OXT", AtomType::O}
            }
        };
        const Residue GLU = {
            "GLU", 'E',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD", AtomType::C_CAR},
                {"OE1", AtomType::O},
                {"OE2", AtomType::O},
                {"OXT", AtomType::O}
            }
        };
        const Residue PHE = {
            "PHE", 'F',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"CD1", AtomType::C_CAR},
                {"CD2", AtomType::C_CAR},
                {"CE1", AtomType::C_CAR},
                {"CE2", AtomType::C_CAR},
                {"CZ", AtomType::C_CAR},
                {"OXT", AtomType::O}
            }
        };
        const Residue GLY = {
            "GLY", 'G',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"OXT", AtomType::O}
            }
        };
        const Residue HIS = {
            "HIS", 'H',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"ND1", AtomType::N_AMD},
                {"CD2", AtomType::C_CAR},
                {"CE1", AtomType::C_CAR},
                {"NE2", AtomType::N_AMD},
                {"OXT", AtomType::O}
            }
        };
        const Residue ILE = {
            "ILE", 'I',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG1", AtomType::C_ALI},
                {"CG2", AtomType::C_ALI},
                {"CD1", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue LYS = {
            "LYS", 'K',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD", AtomType::C_ALI},
                {"CE", AtomType::C_ALI},
                {"NZ", AtomType::N_AMN},
                {"OXT", AtomType::O}
            }
        };
        const Residue LEU = {
            "LEU", 'L',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD1", AtomType::C_ALI},
                {"CD2", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue MET = {
            "MET", 'M',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"SD", AtomType::S},
                {"CE", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue ASN = {
            "ASN", 'N',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"OD1", AtomType::O},
                {"ND2", AtomType::N_AMD},
                {"OXT", AtomType::O}
            }
        };
        const Residue PRO = {
            "PRO", 'P',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue GLN = {
            "GLN", 'Q',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD", AtomType::C_CAR},
                {"OE1", AtomType::O},
                {"NE2", AtomType::N_AMD},
                {"OXT", AtomType::O}
            }
        };
        const Residue ARG = {
            "ARG", 'R',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_ALI},
                {"CD", AtomType::C_ALI},
                {"NE", AtomType::N_AMD},
                {"CZ", AtomType::C_CAR},
                {"NH1", AtomType::N_AMD},
                {"NH2", AtomType::N_AMD},
                {"OXT", AtomType::O}
            }
        };
        const Residue SER = {
            "SER", 'S',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"OG", AtomType::O},
                {"OXT", AtomType::O}
            }
        };
        const Residue THR = {
            "THR", 'T',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"OG1", AtomType::O},
                {"CG2", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue VAL = {
            "VAL", 'V',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG1", AtomType::C_ALI},
                {"CG2", AtomType::C_ALI},
                {"OXT", AtomType::O}
            }
        };
        const Residue TRP = {
            "TRP", 'W',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"CD1", AtomType::C_CAR},
                {"CD2", AtomType::C_CAR},
                {"NE1", AtomType::N_AMD},
                {"CE2", AtomType::C_CAR},
                {"CE3", AtomType::C_CAR},
                {"CZ2", AtomType::C_CAR},
                {"CZ3", AtomType::C_CAR},
                {"CH2", AtomType::C_CAR},
                {"OXT", AtomType::O}
            }
        };
        const Residue TYR = {
            "TYR", 'Y',
            {
                {"N", AtomType::N},
                {"CA", AtomType::C_ALI},
                {"C", AtomType::C_CAR},
                {"O", AtomType::O},
                {"CB", AtomType::C_ALI},
                {"CG", AtomType::C_CAR},
                {"CD1", AtomType::C_CAR},
                {"CD2", AtomType::C_CAR},
                {"CE1", AtomType::C_CAR},
                {"CE2", AtomType::C_CAR},
                {"CZ", AtomType::C_CAR},
                {"OH", AtomType::O},
                {"OXT", AtomType::O}
            }
        };
    };

    const unordered_map<string, const Residue> Residues = {
        {"ALA", ResidueType::ALA},
        {"CYS", ResidueType::CYS},
        {"ASP", ResidueType::ASP},
        {"GLU", ResidueType::GLU},
        {"PHE", ResidueType::PHE},
        {"GLY", ResidueType::GLY},
        {"HIS", ResidueType::HIS},
        {"ILE", ResidueType::ILE},
        {"LYS", ResidueType::LYS},
        {"LEU", ResidueType::LEU},
        {"MET", ResidueType::MET},
        {"ASN", ResidueType::ASN},
        {"PRO", ResidueType::PRO},
        {"GLN", ResidueType::GLN},
        {"ARG", ResidueType::ARG},
        {"SER", ResidueType::SER},
        {"THR", ResidueType::THR},
        {"VAL", ResidueType::VAL},
        {"TRP", ResidueType::TRP},
        {"TYR", ResidueType::TYR}
    };
    
}

#endif
