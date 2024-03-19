#ifndef __TMDET_VALUE_OBJECTS_HBOND__
#define __TMDET_VALUE_OBJECTS_HBOND__

namespace Tmdet::ValueObjects {

    struct HBond {
        double energy = 1e30;
        int toChainIdx = -1;
        int toResIdx = -1;
    };
}

#endif