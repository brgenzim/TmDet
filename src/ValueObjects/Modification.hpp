#ifndef __TMDET_VALUE_OBJECTS_MODIFICATION__
#define __TMDET_VALUE_OBJECTS_MODIFICATION__

#include <string>

namespace Tmdet::ValueObjects {

    struct Modification {
        std::string date;
        std::string descr;
    };
}

#endif