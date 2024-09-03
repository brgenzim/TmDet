#ifndef __TMDET_VALUE_OBJECTS_REGION__
#define __TMDET_VALUE_OBJECTS_REGION__

#include <string>
#include <Types/Region.hpp>

namespace Tmdet::ValueObjects {

    struct Region {
        int beg;
        std::string begi;
        int end;
        std::string endi;
        int rbeg;
        int rend;
        Tmdet::Types::Region type;
    };
}

#endif