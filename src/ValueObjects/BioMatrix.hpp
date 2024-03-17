#ifndef __TMDET_VALUE_OBJECTS_BIOMATRIX__
#define __TMDET_VALUE_OBJECTS_BIOMATRIX__

#include <string>
#include <vector>
#include <ValueObjects/Matrix.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct BioMatrix {
        vector<Matrix> matrices;
        vector<string> deletedChainIds;
    };
}

#endif