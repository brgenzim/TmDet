#ifndef __TMDET_VALUE_OBJECTS_MATRIX__
#define __TMDET_VALUE_OBJECTS_MATRIX__

#include <string>
#include <ValueObjects/TMatrix.hpp>

namespace Tmdet::ValueObjects {

    struct Matrix {
        int id;
        std::string sourceChainId;
        std::string newChainId;
        TMatrix tmatrix;
    };
}

#endif