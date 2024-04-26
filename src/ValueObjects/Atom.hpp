#ifndef __TMDET_VALUE_OBJECTS_ATOM__
#define __TMDET_VALUE_OBJECTS_ATOM__

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <gemmi/model.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Atom {
        gemmi::Atom& gemmi;
        double surface;
        int idx;
        int chainIdx;
        int residueIdx;
        unordered_map<string,any> temp;

        Atom(gemmi::Atom& _gemmi) : 
            gemmi(_gemmi) {}
        ~Atom() {}
    };
}

#endif