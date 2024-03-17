#ifndef __TMDET_VALUE_OBJECTS_REGION__
#define __TMDET_VALUE_OBJECTS_REGION__

#include <string>
#include <vector>
#include <Types/Region.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Region {
        int beg;
        string begi;
        int end;
        string endi;
        int rbeg;
        int rend;
        Tmdet::Types::Region type;
    };
}

#endif