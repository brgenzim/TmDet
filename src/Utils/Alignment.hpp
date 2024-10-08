#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>

// #define __ALIGNMENT_DBG 1 // debug flag for util functions and verbose outputs


using namespace std;

namespace Tmdet::Utils::Alignment {

    extern void alignResidues(Tmdet::ValueObjects::Protein& proteinVO);

}
