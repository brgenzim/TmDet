#ifndef __TMDET_TYPES_ATOM__
#define __TMDET_TYPES_ATOM__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Atom {
        string name;
        double vdw;
    };

    namespace AtomType {
        const double DEFAULT_VDW = 1.8;

        const Atom AG = {"AG", 1.72};
        const Atom AL = {"AL", 0.675};
        const Atom AR = {"AR", 1.88};
        const Atom AS = {"AS", 1.85};
        const Atom AU = {"AU", 1.66};
        const Atom BR = {"BR", 1.85};
        const Atom C = {"C", 1.76};
        const Atom C_ALI = {"C_ALI", 1.87};
        const Atom C_CAR = {"C_CAR", 1.76};
        const Atom C_NUC = {"C_NUC", 1.80};
        const Atom CA = {"CA", 1.26};
        const Atom CD = {"CD", 1.58};
        const Atom CL = {"CL", 1.75};
        const Atom CU = {"CU", 1.40};
        const Atom F = {"F", 1.47};
        const Atom FE = {"FE", 1.47};
        const Atom GA = {"GA", 1.87};
        const Atom H = {"H", 1.20};
        const Atom HE = {"HE", 1.40};
        const Atom HG = {"HG", 1.55};
        const Atom I = {"I", 1.98};
        const Atom IN = {"IN", 1.93};
        const Atom K = {"K", 2.75};
        const Atom KR = {"KR", 2.02};
        const Atom LI = {"LI", 1.82};
        const Atom MG = {"MG", 1.73};
        const Atom MN = {"MN", 0.81};
        const Atom N = {"N", 1.65};
        const Atom N_AMN = {"N_AMN", 1.50};
        const Atom N_AMD = {"N_AMD", 1.65};
        const Atom N_NUC = {"N_NUC", 1.60};
        const Atom NA = {"NA", 2.27};
        const Atom NI = {"NI", 1.63};
        const Atom O = {"O", 1.40};
        const Atom O_CAR = {"O_CAR", 1.50};
        const Atom P = {"P", 1.90};
        const Atom PB = {"PB", 2.02};
        const Atom PD = {"PD", 1.63};
        const Atom PT = {"PT", 1.72};
        const Atom S = {"S", 1.85};
        const Atom SE = {"SE", 1.90};
        const Atom SI = {"SI", 2.10};
        const Atom SN = {"SN", 2.17};
        const Atom TE = {"TE", 2.06};
        const Atom TL = {"TL", 1.96};
        const Atom U = {"U", 1.86};
        const Atom V = {"V", 2.42};
        const Atom XE = {"XE", 2.16};
        const Atom ZN = {"ZN", 1.39};
        const Atom UNK = {"UNK", DEFAULT_VDW};
    };

    const unordered_map<string, const Atom> Atoms = {
        {"AG", AtomType::AG},
        {"AL", AtomType::AL},
        {"AR", AtomType::AR},
        {"AS", AtomType::AS},
        {"AU", AtomType::AU},
        {"BR", AtomType::BR},
        {"C", AtomType::C},
        {"C_ALI", AtomType::C_ALI},
        {"C_CAR", AtomType::C_CAR},
        {"C_NUC", AtomType::C_NUC},
        {"CA", AtomType::CA},
        {"CD", AtomType::CD},
        {"CL", AtomType::CL},
        {"CU", AtomType::CU},
        {"F", AtomType::F},
        {"FE", AtomType::FE},
        {"GA", AtomType::GA},
        {"H", AtomType::H},
        {"HE", AtomType::HE},
        {"HG", AtomType::HG},
        {"I", AtomType::I},
        {"IN", AtomType::IN},
        {"K", AtomType::K},
        {"KR", AtomType::KR},
        {"LI", AtomType::LI},
        {"MG", AtomType::MG},
        {"MN", AtomType::MN},
        {"N", AtomType::N},
        {"N_AMN", AtomType::N_AMN},
        {"N_AMD", AtomType::N_AMD},
        {"N_NUC", AtomType::N_NUC},
        {"NA", AtomType::NA},
        {"NI", AtomType::NI},
        {"O", AtomType::O},
        {"O_CAR", AtomType::O_CAR},
        {"P", AtomType::P},
        {"PB", AtomType::PB},
        {"PD", AtomType::PD},
        {"PT", AtomType::PT},
        {"S", AtomType::S},
        {"SE", AtomType::SE},
        {"SI", AtomType::SI},
        {"SN", AtomType::SN},
        {"TE", AtomType::TE},
        {"TL", AtomType::TL},
        {"U", AtomType::U},
        {"V", AtomType::V},
        {"XE", AtomType::XE},
        {"ZN", AtomType::ZN},
        {"UNK", AtomType::UNK}
    };

}

#endif
