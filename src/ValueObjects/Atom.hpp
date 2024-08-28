#ifndef __TMDET_VALUE_OBJECTS_ATOM__
#define __TMDET_VALUE_OBJECTS_ATOM__

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <gemmi/model.hpp>

namespace Tmdet::ValueObjects {

    struct Atom {
        gemmi::Atom& gemmi;
        double surface = 0.0;
        int idx;
        int chainIdx;
        int residueIdx;
        std::unordered_map<std::string, std::any> temp;

        explicit Atom(gemmi::Atom& _gemmi) :
            gemmi(_gemmi) {
            }

        ~Atom()=default;
    };
}

#endif
