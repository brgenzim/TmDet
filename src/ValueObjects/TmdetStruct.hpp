#ifndef __TMDET_VALUE_OBJECTS_TMDETSTRUCT__
#define __TMDET_VALUE_OBJECTS_TMDETSTRUCT__

#include <string>
#include <vector>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Chain.hpp>
#include <Types/Protein.hpp>

namespace Tmdet::ValueObjects {

    struct TmdetStruct {
        std::string code;
        bool tmp;
        std::string date;
        std::vector<Modification> modifications;
        double qValue;
        Tmdet::Types::Protein type;
        std::string spres;
        std::string pdbkwres;
        BioMatrix bioMatrix;
        std::vector<Membrane> membranes;
        std::vector<Chain> chains;
        // NOTE: gemmi::Structure& reference is short-term and can be invalid.
        //       It may reference freed memory.
        gemmi::Structure gemmi;
        gemmi::cif::Document& document;
        gemmi::NeighborSearch neighbors;

        explicit TmdetStruct(gemmi::Structure& _gemmi, gemmi::cif::Document& _document) :
            code(_gemmi.name),
            gemmi(_gemmi), 
            document(_document) {
        }
        ~TmdetStruct()=default;
    };

    TmdetStruct get(const std::string &inputPath, gemmi::Structure &pdb, gemmi::cif::Document &document);

}

#endif
