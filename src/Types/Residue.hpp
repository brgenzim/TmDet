#ifndef __TMDET_TYPES_RESIDUE__
#define __TMDET_TYPES_RESIDUE__

#include <unordered_map>
#include <map>
#include <string>
#include <functional>
#include <Types/Atom.hpp>

namespace Tmdet::Types {

    const double voronotaMeanMin = -0.567174;
    const double voronotaMeanMax = -0.227304;

    /**
     * mean and sds are taken from voronota 
     * (voromqa_v1_energy_means_and_sds)
    */
    struct AtomData {
        Atom atom;
        bool bb = false;
        double mean;
        double sds;
    };

    struct Residue {
        std::string name;
        char a1;
        float apol = 0.0;
        float hsc = 0.0;
        int nsa = 0;
        std::unordered_map<std::string,AtomData> atoms;
    };

    namespace ResidueType {
        const Residue ALA = {
            "ALA", 'A', 0, 1.6, 1,
            {
                {"N", {AtomType::N, true, -0.366099, 0.260687}},
                {"CA", {AtomType::C_ALI, true, -0.345553, 0.267743}},
                {"C", {AtomType::C_CAR, true, -0.361355, 0.271633}},
                {"O", {AtomType::O, true, -0.345946, 0.258291}},
                {"CB", {AtomType::C_ALI, false, -0.3448, 0.243304}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue CYS = {
            "CYS", 'C', 0, 2.0, 2,
            {
                {"N", {AtomType::N, true, -0.452772, 0.25979}},
                {"CA", {AtomType::C_ALI, true, -0.452816, 0.264808}},
                {"C", {AtomType::C_CAR, true, -0.433833, 0.263039}},
                {"O", {AtomType::O, true, -0.415418, 0.249034}},
                {"CB", {AtomType::C_ALI, false, -0.445575, 0.241479}},
                {"SG", {AtomType::S, false, -0.453821, 0.233319}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue ASP = {
            "ASP", 'D', 0,  -9.2, 4,
            {
                {"N", {AtomType::N, true, -0.284671, 0.219865}},
                {"CA", {AtomType::C_ALI,  true, -0.227304, 0.211296}},
                {"C", {AtomType::C_CAR, true, -0.276127, 0.230521}},
                {"O", {AtomType::O, true, -0.275534, 0.213786}},
                {"CB", {AtomType::C_ALI, false, -0.234205, 0.187007}},
                {"CG", {AtomType::C_CAR, false, -0.23162, 0.164127}},
                {"OD1", {AtomType::O, false, -0.240172, 0.170026}},
                {"OD2", {AtomType::O, false, -0.240172, 0.170026}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue GLU = {
            "GLU", 'E', 0, -8.2, 5,
            {
                {"N", {AtomType::N, true, -0.280941, 0.224805}},
                {"CA", {AtomType::C_ALI, true, -0.228705, 0.219997}},
                {"C", {AtomType::C_CAR, true, -0.274567, 0.234897}},
                {"O", {AtomType::O, true, -0.271903, 0.219933}},
                {"CB", {AtomType::C_ALI, false, -0.228811, 0.189035}},
                {"CG", {AtomType::C_ALI, false, -0.246744, 0.177932}},
                {"CD", {AtomType::C_CAR, false, -0.280044, 0.165997}},
                {"OE1", {AtomType::O, false, -0.289284, 0.184738}},
                {"OE2", {AtomType::O, false, -0.289284, 0.184738}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue PHE = {
            "PHE", 'F', 1, 3.7, 7,
            {
                {"N", {AtomType::N, true, -0.414952, 0.27603}},
                {"CA", {AtomType::C_ALI, true, -0.4315, 0.282044}},
                {"C", {AtomType::C_CAR, true, -0.402941, 0.27988}},
                {"O", {AtomType::O, true, -0.410731, 0.267482}},
                {"CB", {AtomType::C_ALI, false, -0.4593, 0.272564}},
                {"CG", {AtomType::C_CAR, false, -0.567174, 0.402211}},
                {"CD1", {AtomType::C_CAR, false, -0.504217, 0.279886}},
                {"CD2", {AtomType::C_CAR, false, -0.504217, 0.279886}},
                {"CE1", {AtomType::C_CAR, false, -0.517273, 0.287964}},
                {"CE2", {AtomType::C_CAR, false, -0.517273, 0.287964}},
                {"CZ", {AtomType::C_CAR, false, -0.525165, 0.297149}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue GLY = {
            "GLY", 'G', 1, 1.0, 0,
            {
                {"N", {AtomType::N, true, -0.271878, 0.243837}},
                {"CA", {AtomType::C_ALI, true, -0.255495, 0.224329}},
                {"C", {AtomType::C_CAR, true, -0.284029, 0.238359}},
                {"O", {AtomType::O, true, -0.29364, 0.234042}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue HIS = {
            "HIS", 'H', 0, -3.0, 6,
            {
                {"N", {AtomType::N, true, -0.33991, 0.253808}},
                {"CA", {AtomType::C_ALI, true, -0.304916, 0.242045}},
                {"C", {AtomType::C_CAR, true, -0.333108, 0.26356}},
                {"O", {AtomType::O, true, -0.331277, 0.246533}},
                {"CB", {AtomType::C_ALI, false, -0.313616, 0.232557}},
                {"CG", {AtomType::C_CAR, false, -0.321819, 0.23152}},
                {"ND1", {AtomType::N_AMD, false, -0.297591, 0.198435}},
                {"CD2", {AtomType::C_CAR, false, -0.315305, 0.202641}},
                {"CE1", {AtomType::C_CAR, false, -0.292742, 0.17448}},
                {"NE2", {AtomType::N_AMD, false, -0.310276, 0.183333}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue ILE = {
            "ILE", 'I', 1, 3.1, 4,
            {
                {"N", {AtomType::N, true, -0.461542, 0.287267}},
                {"CA", {AtomType::C_ALI, true, -0.479819, 0.309649}},
                {"C", {AtomType::C_CAR, true, -0.450304, 0.293395}},
                {"O", {AtomType::O, true, -0.449107, 0.276286}},
                {"CB", {AtomType::C_ALI, false, -0.518727, 0.337553}},
                {"CG1", {AtomType::C_ALI, false, -0.518525, 0.291153}},
                {"CG2", {AtomType::C_ALI, false, -0.506426, 0.274073}},
                {"CD1", {AtomType::C_ALI, false, -0.514598, 0.274205}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue LYS = {
            "LYS", 'K', 0, -8.8, 5,
            {
                {"N", {AtomType::N, true, -0.288103, 0.231916}},
                {"CA", {AtomType::C_ALI, true, -0.231758, 0.227465}},
                {"C", {AtomType::C_CAR, true, -0.271316, 0.243753}},
                {"O", {AtomType::O, true, -0.268066, 0.229215}},
                {"CB", {AtomType::C_ALI, false, -0.236498, 0.193188}},
                {"CG", {AtomType::C_ALI, false, -0.246717, 0.177021}},
                {"CD", {AtomType::C_ALI, false, -0.26921, 0.162298}},
                {"CE", {AtomType::C_ALI, false, -0.293129, 0.175876}},
                {"NZ", {AtomType::N_AMN, false, -0.334014, 0.223843}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue LEU = {
            "LEU", 'L', 1, 2.8, 4,
            {
                {"N", {AtomType::N, true, -0.430719, 0.273607}},
                {"CA", {AtomType::C_ALI, true, -0.433988, 0.288229}},
                {"C", {AtomType::C_CAR, true, -0.428931, 0.282663}},
                {"O", {AtomType::O, true, -0.412602, 0.269545}},
                {"CB", {AtomType::C_ALI, false, -0.462253, 0.280208}},
                {"CG", {AtomType::C_ALI, false, -0.519465, 0.308807}},
                {"CD1", {AtomType::C_ALI, false, -0.494805, 0.276094}},
                {"CD2", {AtomType::C_ALI, false, -0.484707, 0.275624}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue MET = {
            "MET", 'M', 1, 3.4, 4,
            {
                {"N", {AtomType::N, true, -0.391078, 0.284179}},
                {"CA", {AtomType::C_ALI, true, -0.396421, 0.286734}},
                {"C", {AtomType::C_CAR, true, -0.398077, 0.27971}},
                {"O", {AtomType::O, true, -0.387303, 0.263566}},
                {"CB", {AtomType::C_ALI, false, -0.411255, 0.273216}},
                {"CG", {AtomType::C_ALI, false, -0.439682, 0.275508}},
                {"SD", {AtomType::S, false, -0.458938, 0.273885}},
                {"CE", {AtomType::C_ALI, false, -0.445468, 0.264173}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue ASN = {
            "ASN", 'N', 0, -4.8, 4,
            {
                {"N", {AtomType::N, true, -0.285685, 0.227598}},
                {"CA", {AtomType::C_ALI, true, -0.239685, 0.216018}},
                {"C", {AtomType::C_CAR, true, -0.277629, 0.234746}},
                {"O", {AtomType::O, true, -0.275947, 0.221397}},
                {"CB", {AtomType::C_ALI, false, -0.240536, 0.192117}},
                {"CG", {AtomType::C_CAR, false, -0.243489, 0.177516}},
                {"OD1", {AtomType::O, false, -0.240778, 0.178939}},
                {"ND2", {AtomType::N_AMD, false, -0.248695, 0.169529}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue PRO = {
            "PRO", 'P', 0, -0.2, 3,
            {
                {"N", {AtomType::N, true, -0.265365, 0.239899}},
                {"CA", {AtomType::C_ALI, true, -0.235571, 0.22395}},
                {"C", {AtomType::C_CAR, true, -0.251864, 0.227495}},
                {"O", {AtomType::O, true, -0.274112, 0.216166}},
                {"CB", {AtomType::C_ALI, false, -0.253491, 0.210246}},
                {"CG", {AtomType::C_ALI, false, -0.285889, 0.202525}},
                {"CD", {AtomType::C_ALI, false, -0.27945, 0.208038}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue GLN = {
            "GLN", 'Q', 0, -4.1, 5,
            {
                {"N", {AtomType::N, true, -0.294536, 0.233243}},
                {"CA", {AtomType::C_ALI, true, -0.249327, 0.227324}},
                {"C", {AtomType::C_CAR, true, -0.288729, 0.24214}},
                {"O", {AtomType::O, true, -0.283414, 0.229182}},
                {"CB", {AtomType::C_ALI, false, -0.247521, 0.200981}},
                {"CG", {AtomType::C_ALI, false, -0.244308, 0.189226}},
                {"CD", {AtomType::C_CAR, false, -0.268025, 0.170638}},
                {"OE1", {AtomType::O, false, -0.269828, 0.173859}},
                {"NE2", {AtomType::N_AMD, false, -0.263844, 0.168541}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue ARG = {
            "ARG", 'R', 0, -12.3, 7,
            {
                {"N", {AtomType::N, true, -0.326924, 0.237998}},
                {"CA", {AtomType::C_ALI, true, -0.283062, 0.235648}},
                {"C", {AtomType::C_CAR, true, -0.314033, 0.250329}},
                {"O", {AtomType::O, true, -0.306149, 0.235503}},
                {"CB", {AtomType::C_ALI, false, -0.282617, 0.196033}},
                {"CG", {AtomType::C_ALI, false, -0.286845, 0.183697}},
                {"CD", {AtomType::C_ALI, false, -0.300179, 0.166284}},
                {"NE", {AtomType::N_AMD, false, -0.309583, 0.161413}},
                {"CZ", {AtomType::C_CAR, false, -0.311343, 0.154622}},
                {"NH1", {AtomType::N_AMD, false, -0.299401, 0.14972}},
                {"NH2", {AtomType::N_AMD, false, -0.299401, 0.14972}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue SER = {
            "SER", 'S', 0, 0.6, 2,
            {
                {"N", {AtomType::N, true, -0.305418, 0.242504}},
                {"CA", {AtomType::C_ALI, true, -0.265312, 0.240888}},
                {"C", {AtomType::C_CAR, true, -0.291562, 0.250821}},
                {"O", {AtomType::O, true, -0.292216, 0.235508}},
                {"CB", {AtomType::C_ALI, false, -0.260352, 0.203396}},
                {"OG", {AtomType::O, false, -0.252789, 0.197016}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue THR = {
            "THR", 'T', 0, 1.2, 3,
            {
                {"N", {AtomType::N, true, -0.336243, 0.250223}},
                {"CA", {AtomType::C_ALI, true, -0.30935, 0.25289}},
                {"C", {AtomType::C_CAR, true, -0.340493, 0.264137}},
                {"O", {AtomType::O, true, -0.335836, 0.244642}},
                {"CB", {AtomType::C_ALI, false, -0.301043, 0.220135}},
                {"OG1", {AtomType::O, false, -0.280829, 0.205848}},
                {"CG2", {AtomType::C_ALI, false, -0.318906, 0.220701}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue VAL = {
            "VAL", 'V', 1, 2.6, 3,
            {
                {"N", {AtomType::N, true, -0.446777, 0.287806}},
                {"CA", {AtomType::C_ALI, true, -0.47288, 0.3119}},
                {"C", {AtomType::C_CAR, true, -0.452268, 0.294451}},
                {"O", {AtomType::O, true, -0.441403, 0.275012}},
                {"CB", {AtomType::C_ALI, false, -0.486325, 0.312156}},
                {"CG1", {AtomType::C_ALI, false, -0.476699, 0.272099}},
                {"CG2", {AtomType::C_ALI, false, -0.471999, 0.267144}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue TRP = {
            "TRP", 'W', 1, 1.9, 10,
            {
                {"N", {AtomType::N, true, -0.388018, 0.255962}},
                {"CA", {AtomType::C_ALI, true, -0.407684, 0.260605}},
                {"C", {AtomType::C_CAR, true, -0.390675, 0.268288}},
                {"O", {AtomType::O, true, -0.390417, 0.252122}},
                {"CB", {AtomType::C_ALI, false, -0.427691, 0.260025}},
                {"CG", {AtomType::C_CAR, false, -0.490856, 0.349568}},
                {"CD1", {AtomType::C_CAR, false, -0.396679, 0.216544}},
                {"CD2", {AtomType::C_CAR, false, -0.514994, 0.27336}},
                {"NE1", {AtomType::N_AMD, false, -0.418557, 0.220445}},
                {"CE2", {AtomType::C_CAR, false, -0.495801, 0.273607}},
                {"CE3", {AtomType::C_CAR, false, -0.5046, 0.26119}},
                {"CZ2", {AtomType::C_CAR, false, -0.451098, 0.231715}},
                {"CZ3", {AtomType::C_CAR, false, -0.503666, 0.264548}},
                {"CH2", {AtomType::C_CAR, false, -0.479774, 0.256487}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue TYR = {
            "TYR", 'Y', 1, -0.7, 8,
            {
                {"N", {AtomType::N, true, -0.395839, 0.270795}},
                {"CA", {AtomType::C_ALI, true, -0.397171, 0.266279}},
                {"C", {AtomType::C_CAR, true, -0.394996, 0.277939}},
                {"O", {AtomType::O, true, -0.392777, 0.260632}},
                {"CB", {AtomType::C_ALI, false, -0.424335, 0.262804}},
                {"CG", {AtomType::C_CAR, false, -0.462158, 0.352833}},
                {"CD1", {AtomType::C_CAR, false, -0.418049, 0.235199}},
                {"CD2", {AtomType::C_CAR, false, -0.418049, 0.235199}},
                {"CE1", {AtomType::C_CAR, false, -0.406279, 0.217514}},
                {"CE2", {AtomType::C_CAR, false, -0.406279, 0.217514}},
                {"CZ", {AtomType::C_CAR, false, -0.420874, 0.236456}},
                {"OH", {AtomType::O, false, -0.382486, 0.195056}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        const Residue UNK = {
            "UNK", 'X', 0, 0, 5,
            {
                {"N", {AtomType::N, true, -0.395839, 0.270795}},
                {"CA", {AtomType::C_ALI, true, -0.397171, 0.266279}},
                {"C", {AtomType::C_CAR, true, -0.394996, 0.277939}},
                {"O", {AtomType::O, true, -0.392777, 0.260632}},
                {"CB", {AtomType::C_ALI, false, -0.424335, 0.262804}},
                {"OXT", {AtomType::O, true, 0, 0}}
            }
        };
        extern std::map<std::string, Residue, std::less<>> ChemicalCompoundDictionary;

        extern Residue getResidue(const std::string& threeLetterCode);
    };

    const std::unordered_map<std::string, const Residue> Residues = {
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
        {"TYR", ResidueType::TYR},
        {"UNK", ResidueType::UNK}
    };
    
}

#endif
