#ifndef __TMDET_VALUE_OBJECTS_MATRIX__
#define __TMDET_VALUE_OBJECTS_MATRIX__

#include <string>
#include <vector>
#include <ValueObjects/TMatrix.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Matrix {
        int id;
        string sourceChainId;
        string newChainId;
        TMatrix tmatrix;
    };
}

#endif