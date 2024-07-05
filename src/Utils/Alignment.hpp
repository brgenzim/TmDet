#ifndef __TMDET_UTILS_ALIGNMENT__
#define __TMDET_UTILS_ALIGNMENT__

#include <string>
#include <vector>

using namespace std;

namespace Tmdet::Utils::Alignment {

    extern vector<string> alignSequences(vector<string> &query, vector<string> &target);

}

#endif
