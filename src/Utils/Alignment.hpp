#ifndef __TMDET_UTILS_ALIGNMENT__
#define __TMDET_UTILS_ALIGNMENT__

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

// #define __ALIGNMENT_DBG 1 // debug flag for util functions and verbose outputs


using namespace std;

namespace Tmdet::Utils::Alignment {

    extern void alignResidues(const Tmdet::ValueObjects::TmdetStruct& tmdetVO);

}

#endif
