#ifndef __TMDET_VALUE_OBJECTS_BIOMATRIX__
#define __TMDET_VALUE_OBJECTS_BIOMATRIX__

#include <string>
#include <vector>
#include <ValueObjects/Matrix.hpp>

namespace Tmdet::ValueObjects {

    struct BioMatrix {
        std::vector<Matrix> matrices;
        std::vector<std::string> deletedChainIds;
    };
}

#endif