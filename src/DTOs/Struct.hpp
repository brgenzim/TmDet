#ifndef __TMDET_DTOS_STRUCT__
#define __TMDET_DTOS_STRUCT__

#include <string>
#include <vector>
#include <ValueObjects/Struct.hpp>

using namespace std;

namespace Tmdet::DTOS {

    struct Struct {
        static void read(Tmdet::ValueObjects::Struct& tmdetVO, string path);
        static void write(Tmdet::ValueObjects::Struct& tmdetVO, string path);
    };
}

#endif