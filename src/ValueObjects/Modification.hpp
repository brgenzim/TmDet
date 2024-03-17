#ifndef __TMDET_VALUE_OBJECTS_MODIFICATION__
#define __TMDET_VALUE_OBJECTS_MODIFICATION__

#include <string>
#include <vector>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Modification {
        string date;
        string descr;
    };
}

#endif