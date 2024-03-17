#ifndef __TMDET_VALUE_OBJECTS_CHAIN__
#define __TMDET_VALUE_OBJECTS_CHAIN__

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Region.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Chain {
        string id;
        bool selected;
        int numtm;
        string seq;
        vector<Tmdet::ValueObjects::Region> regions;
        Tmdet::Types::Chain type;
    };
}

#endif