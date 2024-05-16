#ifndef __TMDET_VALUE_OBJECTS_TMDETSTRUCT__
#define __TMDET_VALUE_OBJECTS_TMDETSTRUCT__

#include <string>
#include <vector>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Chain.hpp>
#include <Types/Protein.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct TmdetStruct {
        string code;
        bool tmp;
        string date;
        vector<Modification> modifications;
        double qValue;
        Tmdet::Types::Protein type;
        string spres;
        string pdbkwres;
        BioMatrix bioMatrix;
        vector<Membrane> membranes;
        vector<Chain> chains;
        gemmi::Structure& gemmi;
        gemmi::cif::Document& document;
        gemmi::NeighborSearch neighbors;
        TmdetStruct(gemmi::Structure& _gemmi, gemmi::cif::Document& _document) : gemmi(_gemmi), document(_document) {
        }
        ~TmdetStruct() {}
    };
}

#endif